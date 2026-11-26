#!/usr/bin/env python3
"""find_IS91.py v2 -- localise les insertions IS91 (ter ... ori) dans les contigs d'un barcode.
Python seul (pas de BLAST). Usage:
    find_IS91.py --header
    find_IS91.py BARCODE contigs.fa  > analysis/is91/BARCODE.tsv
Lit analysis/samples.tsv et analysis/contig_assignment_meta.tsv (variable IS91_ANALYSIS pour changer le dossier).
Conventions: tout est exprimé dans le sens de l'IS (ter -> ori). strand 'plus' = ter->ori vers les
coordonnées croissantes du contig. Coordonnées 1-based incluses."""
import sys, os, csv
from collections import defaultdict

SEQDIR = os.environ.get("IS91_SEQDIR", "sequences")
ANA = os.environ.get("IS91_ANALYSIS", "analysis")
REF_FILES = [("CHR", "GCF_056501015.1_ASM5650101v1_genomic.fna"), ("POX38", "pOX38Km.fa"),
             ("pCP54", "pCP54.fa"), ("pCP55", "pCP55.fa"), ("pCP7", "pCP7.fa")]
TARGETS, DONORS = {"CHR", "POX38"}, {"pCP54", "pCP55", "pCP7"}
EXPECTED = {"pCP7": 945, "pCP54": 1984, "pCP55": 1984}   # ter(début) -> ori(fin) attendu
TARGET_MAP = {"pxo38km": "pOX38Km", "pox38km": "pOX38Km", "ecoli_chr": "chromosome"}
SEED, SEEDF = 20, 24
MIN_COV, MIN_ID = 0.80, 0.88
MAX_SPAN, TANDEM_GAP, JUNC_GAP, FLANK = 20000, 40, 20, 150
SN = {"+": "plus", "-": "minus"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

COLS = ["barcode", "is_type", "exp_target", "passage", "bio_rep", "clone", "contig", "contig_len",
        "assignation", "on_target", "event_class", "oet", "array_id", "copy_index", "n_copies", "strand",
        "ins_start", "ins_end", "ter_start", "ter_end", "ori_start", "ori_end", "ter_to_ori_len",
        "expected_len", "len_vs_expected", "inner_len", "arr_start", "arr_end", "target_4nt", "ori_tail4",
        "left_class", "left_map", "right_class", "right_map", "ref_gap", "ins_strand_on_ref",
        "upstream_10nt", "downstream_10nt", "inner_seq"]
if len(sys.argv) > 1 and sys.argv[1] == "--header":
    print("\t".join(COLS)); sys.exit()

def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; d[n] = []
        elif n is not None: d[n].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}

def sw(q, t, match=2, mis=-3, gap=-4):
    n, m = len(q), len(t)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best = (0, 0, 0)
    for i in range(1, n + 1):
        qi, Hi, Hp = q[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            v = max(0, Hp[j - 1] + (match if qi == t[j - 1] else mis), Hp[j] + gap, Hi[j - 1] + gap)
            Hi[j] = v
            if v > best[0]: best = (v, i, j)
    sc, i, j = best
    qe, te, mt, al = i, j, 0, 0
    while i > 0 and j > 0 and H[i][j] > 0:
        v = H[i][j]
        if v == H[i - 1][j - 1] + (match if q[i - 1] == t[j - 1] else mis):
            mt += q[i - 1] == t[j - 1]; i -= 1; j -= 1
        elif v == H[i - 1][j] + gap: i -= 1
        else: j -= 1
        al += 1
    return sc, i, qe, j, te, mt, al      # qs, qe (0-based, qe exclusif), ts, te (exclusif)

def find_hits(seq, q, kind):
    L, out = len(q), []
    for strand in "+-":
        qs = q if strand == "+" else rc(q)
        diags = set()
        for o in range(0, L - SEED + 1, 8):
            sd = qs[o:o + SEED]; i = seq.find(sd); c = 0
            while i != -1 and c < 200:
                diags.add(i - o); c += 1; i = seq.find(sd, i + 1)
        reps = []
        for d in sorted(diags):
            if not reps or d - reps[-1] > 8: reps.append(d)
        for d in reps:
            ws, we = max(0, d - 20), min(len(seq), d + L + 20)
            sc, q0, q1, t0, t1, mt, al = sw(qs, seq[ws:we])
            if al == 0 or (q1 - q0) / L < MIN_COV or mt / al < MIN_ID: continue
            ts, te = ws + t0, ws + t1                 # 0-based, te exclusif
            if strand == "+": lo, hi = ts - q0, te + (L - q1)    # étendue de la requête entière
            else:             lo, hi = ts - (L - q1), te + q0
            out.append(dict(kind=kind, strand=strand, s=ts + 1, e=te, score=sc, lo=max(0, lo),
                            hi=min(len(seq), hi), ident=mt / al, L=L, q0=q0, q1=q1, ts=ts, te=te))
    out.sort(key=lambda h: -h["score"]); kept = []
    for h in out:
        if any(k["kind"] == h["kind"] and min(k["e"], h["e"]) - max(k["s"], h["s"]) + 1 > 0.5 * (h["e"] - h["s"] + 1)
               for k in kept): continue
        kept.append(h)
    return kept

# ---------- références ----------
REFSEQ = {}
for lab, fn in REF_FILES:
    p = os.path.join(SEQDIR, fn)
    if os.path.exists(p):
        for name, s in read_fasta(p).items(): REFSEQ[f"{lab}|{name}"] = s

def map_flank(seq, side):
    n = len(seq)
    if n < 30: return None
    groups = defaultdict(list)
    for strand in "+-":
        qs = seq if strand == "+" else rc(seq)
        for o in range(0, n - SEEDF + 1, 6):
            sd = qs[o:o + SEEDF]
            for key, fw in REFSEQ.items():
                i = fw.find(sd); c = 0
                while i != -1 and c < 20:
                    orig = o if strand == "+" else n - o - SEEDF
                    groups[(key, strand)].append((i - o, orig)); c += 1; i = fw.find(sd, i + 1)
    best = None
    for (key, strand), lst in groups.items():
        lst.sort(); cl = [lst[0:1]]
        for x in lst[1:]:
            if x[0] - cl[-1][-1][0] <= 6: cl[-1].append(x)
            else: cl.append([x])
        for c in cl:
            if best is None or len(c) > len(best[2]): best = (key, strand, c)
    if best is None or len(best[2]) < 2: return None
    key, strand, c = best
    diag, _ = (max(c, key=lambda x: x[1]) if side == "L" else min(c, key=lambda x: x[1]))
    if strand == "+": pos = n + diag if side == "L" else diag + 1
    else:             pos = diag + 1 if side == "L" else n + diag
    lab, name = key.split("|", 1)
    return dict(label=lab, ref=name, strand=strand, pos=pos, nseeds=len(c))

def fmt_map(m): return "NA" if not m else f"{m['label']}:{m['ref']}:{m['pos']}({m['strand']},{m['nseeds']}s)"
def fl_class(m, flank):
    if len(flank) < 30: return "edge"
    if not m: return "unmapped"
    return "target" if m["label"] in TARGETS else ("donor" if m["label"] in DONORS else "other")

# ---------- métadonnées ----------
bc, fasta = sys.argv[1], sys.argv[2]
sp = os.path.join(ANA, "samples.tsv"); meta = {}
if os.path.exists(sp):
    for r in csv.DictReader(open(sp, encoding="utf-8-sig"), delimiter="\t"):
        if r["sample"].strip() == bc: meta = r
is_type = meta.get("type_IS", "NA"); exp_t = TARGET_MAP.get(meta.get("cible_IS", "").lower(), meta.get("cible_IS", "NA"))
passage, brep, clone = meta.get("passage", "NA"), meta.get("N_biol", "NA"), meta.get("clone", "NA")
assign = {}
ap = os.path.join(ANA, "contig_assignment_meta.tsv")
if os.path.exists(ap):
    for r in csv.DictReader(open(ap, encoding="utf-8-sig"), delimiter="\t"):
        if r["barcode"] == bc: assign[r["contig"]] = r["assignation"]
expected = EXPECTED.get(is_type, "NA")

ORI = list(read_fasta(os.path.join(SEQDIR, "ori.fa")).values())[0]
TER = list(read_fasta(os.path.join(SEQDIR, "ter.fa")).values())[0]
contigs = read_fasta(fasta)
rows = []

for cname, seq in contigs.items():
    hits = find_hits(seq, ORI, "ORI") + find_hits(seq, TER, "TER")
    if not hits: continue
    units = []   # dict(type, strand, ter, ori)
    for strand in "+-":
        ev = sorted([h for h in hits if h["strand"] == strand],
                    key=(lambda h: h["s"]) if strand == "+" else (lambda h: -h["e"]))
        pend, orphans = None, []
        for h in ev:
            if h["kind"] == "TER":
                if pend: orphans.append(pend)
                pend = h
            elif pend and max(pend["e"], h["e"]) - min(pend["s"], h["s"]) + 1 <= MAX_SPAN:
                units.append(dict(type="complete", strand=strand, ter=pend, ori=h)); pend = None
            else:
                if pend: orphans.append(pend); pend = None
                orphans.append(h)
        if pend: orphans.append(pend)
        orphans.sort(key=(lambda h: h["s"]) if strand == "+" else (lambda h: -h["e"]))
        used = set()
        for a, b in zip(orphans, orphans[1:]):
            if a["kind"] == "ORI" and b["kind"] == "TER":
                gap = (b["s"] - a["e"] - 1) if strand == "+" else (a["s"] - b["e"] - 1)
                if -10 <= gap <= JUNC_GAP:
                    units.append(dict(type="junction_donor", strand=strand, ori=a, ter=b)); used |= {id(a), id(b)}
        for h in orphans:
            if id(h) not in used:
                units.append(dict(type=h["kind"].lower() + "_only", strand=strand,
                                  ori=h if h["kind"] == "ORI" else None, ter=h if h["kind"] == "TER" else None))
    # --- tandems ---
    arrays = []
    for strand in "+-":
        us = [u for u in units if u["type"] == "complete" and u["strand"] == strand]
        key = (lambda u: u["ter"]["lo"]) if strand == "+" else (lambda u: -u["ter"]["hi"])
        us.sort(key=key); cur = []
        for u in us:
            if cur:
                p = cur[-1]
                gap = (u["ter"]["lo"] - p["ori"]["hi"]) if strand == "+" else (p["ori"]["lo"] - u["ter"]["hi"])
                if -10 <= gap <= TANDEM_GAP: cur.append(u); continue
                arrays.append(cur)
            cur = [u]
        if cur: arrays.append(cur)
    asg = assign.get(cname, "NA")
    on = "NA" if exp_t == "NA" or asg == "NA" else ("yes" if exp_t == asg else "no")
    base = [bc, is_type, exp_t, passage, brep, clone, cname, len(seq), asg, on]

    for ai, arr in enumerate(arrays, 1):
        st = arr[0]["strand"]
        if st == "+": a_lo, a_hi = arr[0]["ter"]["lo"], arr[-1]["ori"]["hi"]
        else:         a_lo, a_hi = arr[-1]["ori"]["lo"], arr[0]["ter"]["hi"]
        if st == "+": Lf, Rf = seq[max(0, a_lo - FLANK):a_lo], seq[a_hi:a_hi + FLANK]
        else:         Lf, Rf = rc(seq[a_hi:a_hi + FLANK]), rc(seq[max(0, a_lo - FLANK):a_lo])
        lm, rm = map_flank(Lf, "L"), map_flank(Rf, "R")
        lc, rc_ = fl_class(lm, Lf), fl_class(rm, Rf)
        span = a_hi - a_lo; n = len(arr)
        ref_gap, ins_ref = "NA", "NA"
        if lm and rm and lm["ref"] == rm["ref"] and lm["strand"] == rm["strand"]:
            ref_gap = (rm["pos"] - lm["pos"] - 1) if lm["strand"] == "+" else (lm["pos"] - rm["pos"] - 1)
        hm = lm or rm
        if hm: ins_ref = SN[hm["strand"]]
        extra = span - n * expected if isinstance(expected, int) else None
        if lc == "donor": oet = "yes"
        elif lc in ("unmapped", "edge") and extra is not None and extra > 100: oet = "maybe"
        else: oet = "no"
        if n > 1: ecl = "tandem"
        elif oet == "yes": ecl = "OET"
        elif oet == "maybe": ecl = "OET_candidate"
        elif lc == "target" and rc_ == "target":
            ecl = "simple" if ref_gap != "NA" and abs(ref_gap) <= 6 else "rearranged_junction"
        else: ecl = "other_junction"
        ordered = arr
        for ci, u in enumerate(ordered, 1):
            t, o = u["ter"], u["ori"]
            if st == "+":
                s_, e_ = t["lo"] + 1, o["hi"]; inner = seq[t["e"]:o["s"] - 1]; tail = seq[o["hi"] - 4:o["hi"]]
                up, dn = seq[max(0, s_ - 11):s_ - 1], seq[e_:e_ + 10]
            else:
                s_, e_ = o["lo"] + 1, t["hi"]; inner = rc(seq[o["e"]:t["s"] - 1]); tail = rc(seq[o["lo"]:o["lo"] + 4])
                up, dn = rc(seq[e_:e_ + 10]), rc(seq[max(0, s_ - 11):s_ - 1])
            ln = e_ - s_ + 1
            rows.append(base + [ecl, oet, ai, ci, n, SN[st], s_, e_, t["s"], t["e"], o["s"], o["e"], ln,
                                expected, (ln - expected) if isinstance(expected, int) else "NA", len(inner),
                                a_lo + 1, a_hi, rc(tail), tail, lc, fmt_map(lm), rc_, fmt_map(rm), ref_gap,
                                ins_ref, up, dn, inner])
    for u in units:
        if u["type"] == "complete": continue
        h = u["ter"] or u["ori"]
        t, o = u["ter"], u["ori"]
        rows.append(base + [u["type"], "NA", "NA", 1, 1, SN[u["strand"]],
                            min(x["s"] for x in (t, o) if x), max(x["e"] for x in (t, o) if x),
                            t["s"] if t else "NA", t["e"] if t else "NA", o["s"] if o else "NA", o["e"] if o else "NA",
                            "NA", expected, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                            "NA", "NA", "NA", "NA"])

def sk(r):
    v = r[16]; return (r[6], int(v) if str(v).isdigit() else 0)
for r in sorted(rows, key=sk):
    print("\t".join(map(str, r)))
