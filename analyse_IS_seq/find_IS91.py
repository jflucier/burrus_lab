#!/usr/bin/env python3
import sys, csv, os, subprocess, tempfile
from collections import defaultdict

COLS = ["barcode","contig","contig_len","contig_assignation","is_type","exp_target","on_target","event_type",
        "array_id","copy_index","n_copies","strand","ins_start","ins_end",
        "ter_start","ter_end","ori_start","ori_end","ori_to_ter_len","inner_len",
        "arr_start","arr_end","upstream_10nt","downstream_10nt","target_4nt",
        "ref_left","ref_right","ins_strand_on_ref","inner_seq"]

if len(sys.argv) > 1 and sys.argv[1] == "--header":
    print("\t".join(COLS)); sys.exit()

bc, fasta = sys.argv[1], sys.argv[2]
ORI, TER = "sequences/ori.fa", "sequences/ter.fa"
REFDB = "analysis/db/refs"
MIN_ID, MIN_COV = 90.0, 0.80
MAX_SPAN = 20000     # ter->ori max pour former une paire
TANDEM_GAP = 50      # écart max entre copies adjacentes
FLANK, SHORT = 150, 10
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}

meta = {}
if os.path.exists("analysis/contig_assignment.tsv"):
    for r in csv.DictReader(open("analysis/contig_assignment.tsv"), delimiter="\t"):
        meta[(r["barcode"], r["contig"])] = (r["assignation"],
            r.get("Type_IS", "NA"), r.get("Cible_IS", "NA"))

seqs = read_fasta(fasta)
ori = list(read_fasta(ORI).values())[0]
ter = list(read_fasta(TER).values())[0]
tmp = tempfile.mkdtemp()
qf = f"{tmp}/q.fa"
open(qf, "w").write(f">ORI\n{ori}\n>TER\n{ter}\n")

FMT = "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore"
out = subprocess.run(["blastn", "-task", "blastn-short", "-query", qf, "-subject", fasta,
                      "-outfmt", FMT, "-evalue", "1e-5", "-dust", "no"],
                     capture_output=True, text=True, check=True).stdout

cand = []
for l in out.splitlines():
    f = l.split("\t")
    qlen, pid, aln = int(f[1]), float(f[4]), int(f[5])
    s, e = int(f[8]), int(f[9])
    if pid < MIN_ID or aln / qlen < MIN_COV: continue
    cand.append(dict(kind=f[0], contig=f[2], s=min(s, e), e=max(s, e),
                     strand="+" if s <= e else "-", bit=float(f[11])))
cand.sort(key=lambda h: -h["bit"])
kept = []   # déduplication (les palindromes d'ori peuvent donner un hit sur l'autre brin)
for h in cand:
    if any(k["kind"] == h["kind"] and k["contig"] == h["contig"] and
           min(k["e"], h["e"]) - max(k["s"], h["s"]) + 1 > 0.5 * (h["e"] - h["s"] + 1)
           for k in kept): continue
    kept.append(h)

by = defaultdict(list)
for h in kept: by[h["contig"]].append(h)

# --- appariement ter -> ori dans le sens de l'IS ---
units = []
for contig, hs in by.items():
    for strand in "+-":
        ev = [h for h in hs if h["strand"] == strand]
        ev.sort(key=(lambda h: h["s"]) if strand == "+" else (lambda h: -h["e"]))
        def orphan(h):
            units.append(dict(contig=contig, strand=strand, type=h["kind"].lower() + "_only",
                              ter=h if h["kind"] == "TER" else None,
                              ori=h if h["kind"] == "ORI" else None, s=h["s"], e=h["e"]))
        pend = None
        for h in ev:
            if h["kind"] == "TER":
                if pend: orphan(pend)
                pend = h
            else:
                if pend and max(pend["e"], h["e"]) - min(pend["s"], h["s"]) + 1 <= MAX_SPAN:
                    units.append(dict(contig=contig, strand=strand, type="complete",
                                      ter=pend, ori=h,
                                      s=min(pend["s"], h["s"]), e=max(pend["e"], h["e"])))
                else:
                    if pend: orphan(pend)
                    orphan(h)
                pend = None
        if pend: orphan(pend)

# --- regroupement en tandems ---
groups = []
for contig in by:
    for strand in "+-":
        us = [u for u in units if u["contig"] == contig and u["strand"] == strand
              and u["type"] == "complete"]
        us.sort(key=(lambda u: u["s"]) if strand == "+" else (lambda u: -u["e"]))
        cur = []
        for u in us:
            if cur:
                p = cur[-1]
                gap = u["s"] - p["e"] - 1 if strand == "+" else p["s"] - u["e"] - 1
                if -10 <= gap <= TANDEM_GAP:
                    cur.append(u); continue
                groups.append(cur)
            cur = [u]
        if cur: groups.append(cur)
for u in units:
    if u["type"] != "complete": groups.append([u])

# --- cartographie des flancs sur les références ---
fq = open(f"{tmp}/flanks.fa", "w")
for gi, g in enumerate(groups):
    sq = seqs[g[0]["contig"]]
    a_s, a_e = min(u["s"] for u in g), max(u["e"] for u in g)
    L, R = sq[max(0, a_s - 1 - FLANK):a_s - 1], sq[a_e:a_e + FLANK]
    if len(L) >= 30: fq.write(f">{gi}_L\n{L}\n")
    if len(R) >= 30: fq.write(f">{gi}_R\n{R}\n")
fq.close()
best = {}
if os.path.getsize(f"{tmp}/flanks.fa") > 0:
    o = subprocess.run(["blastn", "-query", f"{tmp}/flanks.fa", "-db", REFDB, "-evalue", "1e-10",
                        "-outfmt", "6 qseqid qlen sseqid pident length qstart qend sstart send bitscore"],
                       capture_output=True, text=True, check=True).stdout
    for l in o.splitlines():
        f = l.split("\t"); q = f[0]; ql = int(f[1]); pid = float(f[3])
        qs, qe, ss, se, bit = int(f[5]), int(f[6]), int(f[7]), int(f[8]), float(f[9])
        if pid < MIN_ID: continue
        side = q.split("_")[1]
        if side == "L" and qe < ql - 5: continue   # doit toucher l'insertion
        if side == "R" and qs > 5: continue
        pos = se if side == "L" else ss
        if q not in best or bit > best[q][0]:
            best[q] = (bit, f[2], pos, "+" if ss <= se else "-")

# --- sorties ---
os.makedirs("analysis/is91", exist_ok=True)
fa = open(f"analysis/is91/{bc}.fa", "w")
rows = []
for gi, g in enumerate(groups):
    contig = g[0]["contig"]; sq = seqs[contig]
    a_s, a_e = min(u["s"] for u in g), max(u["e"] for u in g)
    strand = g[0]["strand"]
    lb, rb = best.get(f"{gi}_L"), best.get(f"{gi}_R")
    ref_left = f"{lb[1]}:{lb[2]}" if lb else "NA"
    ref_right = f"{rb[1]}:{rb[2]}" if rb else "NA"
    hs = (lb or rb)
    ref_strand = "NA" if not hs else (strand if hs[3] == "+" else ("-" if strand == "+" else "+"))
    ordered = sorted(g, key=(lambda u: u["s"]) if strand == "+" else (lambda u: -u["e"]))
    for ci, u in enumerate(ordered, 1):
        s, e = u["s"], u["e"]
        t, o_ = u["ter"], u["ori"]
        if u["type"] == "complete":
            inner = sq[t["e"]:o_["s"] - 1] if strand == "+" else rc(sq[o_["e"]:t["s"] - 1])
            insert = sq[s - 1:e] if strand == "+" else rc(sq[s - 1:e])
            span, inl = e - s + 1, len(inner)
            fa.write(f">{bc}|{contig}|{s}-{e}|{strand}|{span}bp\n{insert}\n")
        else:
            inner, span, inl = "NA", "NA", "NA"
        if strand == "+":
            up, down = sq[max(0, s - 1 - SHORT):s - 1], sq[e:e + SHORT]
        else:
            up, down = rc(sq[e:e + SHORT]), rc(sq[max(0, s - 1 - SHORT):s - 1])

        SN = {"+": "plus", "-": "minus"}
        strand_lbl = SN[strand]
        ref_strand_lbl = SN[ref_strand]
        asg, ist, tgt = meta.get((bc, contig), ("NA", "NA", "NA"))
        on = "yes" if (("pXO38" in tgt and asg == "pOX38Km") or
                       ("Chromosome" in tgt and asg == "chromosome")) else "no"
        rows.append([bc, contig, len(sq), asg, ist, tgt, on, u["type"],
                     gi + 1, ci, len(g), strand_lbl, s, e,
                     t["s"] if t else "NA", t["e"] if t else "NA",
                     o_["s"] if o_ else "NA", o_["e"] if o_ else "NA",
                     span, inl, a_s, a_e, up, down, down[:4],
                     ref_left, ref_right, ref_strand_lbl, inner])
for r in sorted(rows, key=lambda r: (r[1], r[12])):
    print("\t".join(map(str, r)))