#!/usr/bin/env python3
"""find_IS91.py v3 -- localise les insertions IS91 (ter ... ori) et les copies du phage Mu dans les contigs d'un barcode.
Python seul (pas de BLAST). Usage:
    find_IS91.py --header | --header-mu
    find_IS91.py BARCODE contigs.fa  > analysis/is91/BARCODE.tsv
Lit analysis/samples.tsv et analysis/contig_assignment_meta.tsv (variable IS91_ANALYSIS pour changer le dossier).
Mu: si sequences/Mu.fa existe, les copies de Mu sont détectées (écrites dans analysis/mu/BARCODE.tsv) et
les colonnes mu_in_contig / mu_dist sont ajoutées à la table des insertions IS91.
Conventions: tout est exprimé dans le sens de l'IS (ter -> ori). strand 'plus' = ter->ori vers les
coordonnées croissantes du contig. Coordonnées 1-based incluses."""
import sys, os, csv, re
from collections import defaultdict

SEQDIR = os.environ.get("IS91_SEQDIR", "sequences")
ANA = os.environ.get("IS91_ANALYSIS", "analysis")
REF_FILES = [("CHR", "GCF_056501015.1_ASM5650101v1_genomic.fna"), ("POX38", "pOX38Km.fa"),
             ("pCP54", "pCP54.fa"), ("pCP55", "pCP55.fa"), ("pCP7", "pCP7.fa"), ("MU", "Mu.fa")]
TARGETS, DONORS = {"CHR", "POX38"}, {"pCP54", "pCP55", "pCP7"}
EXPECTED = {"pCP7": 945, "pCP54": 1988, "pCP55": 1988}   # ter(début) -> ori(fin) attendu (ori recouvre ter de 4 nt dans pCP54/55)
TARGET_MAP = {"pxo38km": "pOX38Km", "pox38km": "pOX38Km", "ecoli_chr": "chromosome"}
SEED, SEEDF = 20, 24
MIN_COV, MIN_ID = 0.80, 0.88
MAX_SPAN, TANDEM_GAP, JUNC_GAP = 20000, 40, 20
FLANK = int(os.environ.get("IS91_FLANK", "1000"))   # flancs plus longs = meilleure résolution des répétitions
SN = {"+": "plus", "-": "minus"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

COLS = ["barcode", "is_type", "exp_target", "passage", "bio_rep", "clone", "contig", "contig_len",
        "assignation", "on_target", "event_class", "oet", "array_id", "copy_index", "n_copies", "strand",
        "ins_start", "ins_end", "ter_start", "ter_end", "ori_start", "ori_end", "ter_to_ori_len",
        "expected_len", "len_vs_expected", "inner_len", "arr_start", "arr_end", "target_4nt", "ori_tail4",
        "left_class", "left_map", "right_class", "right_map", "ref_gap", "ref_gap_seq", "ins_strand_on_ref", "array_open",
        "mu_in_contig", "mu_dist", "upstream_10nt", "downstream_10nt", "inner_seq"]
MU_COLS = ["barcode", "is_type", "exp_target", "passage", "bio_rep", "clone", "contig", "contig_len",
           "mu_start", "mu_end", "mu_len", "mu_strand", "mu_cov", "inverted_bp", "state", "left_map", "right_map",
           "tsd_len", "tsd_seq", "site_ref", "site_start", "site_end", "n_mu_in_contig", "snp_profile", "snp_diff"]
if len(sys.argv) > 1 and sys.argv[1] == "--header":
    print("\t".join(COLS)); sys.exit()
if len(sys.argv) > 1 and sys.argv[1] == "--header-mu":
    print("\t".join(MU_COLS)); sys.exit()

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
    allc = []
    for (key, strand), lst in groups.items():
        lst.sort(); cl = [lst[0:1]]
        for x in lst[1:]:
            if x[0] - cl[-1][-1][0] <= 6: cl[-1].append(x)
            else: cl.append([x])
        for c in cl:
            bonus = 3 if key.split("|")[0] == EXPLAB else 0      # préférence pour le réplicon cible attendu
            allc.append((len(c) + bonus, len(c), key, strand, c))
    if not allc: return None
    allc.sort(key=lambda x: -x[0])
    _, nbest, key, strand, c = allc[0]
    if nbest < 2: return None
    near = [x for x in allc if x[1] >= 0.8 * nbest]
    nloci = len(near)                                              # >1 : séquence répétée, position ambiguë
    def mk(x):
        _, nb, k, st, cc = x
        diag, _o = (max(cc, key=lambda y: y[1]) if side == "L" else min(cc, key=lambda y: y[1]))
        if st == "+": ps = n + diag if side == "L" else diag + 1
        else:         ps = diag + 1 if side == "L" else n + diag
        ps = (ps - 1) % len(REFSEQ[k]) + 1                         # références circulaires
        lb, nm = k.split("|", 1)
        return dict(label=lb, ref=nm, strand=st, pos=ps, nseeds=nb, nloci=nloci)
    out = mk(allc[0]); out["alts"] = [mk(x) for x in near]
    return out

def fmt_map(m): return "NA" if not m else f"{m['label']}:{m['ref']}:{m['pos']}({m['strand']},{m['nseeds']}s,{m['nloci']}loci)"
def fl_class(m, flank):
    if len(flank) < 30: return "edge"
    if not m: return "unmapped"
    if m["label"] == "MU": return "mu"
    return "target" if m["label"] in TARGETS else ("donor" if m["label"] in DONORS else "other")

def reconcile(lm, rm):
    """Si un flanc est ambigu (répétition), choisit les loci des deux flancs qui sont consécutifs."""
    if not lm or not rm: return lm, rm
    best = None
    for a in lm["alts"]:
        for b in rm["alts"]:
            if (a["label"], a["ref"], a["strand"]) != (b["label"], b["ref"], b["strand"]): continue
            RL = len(REFSEQ[a["label"] + "|" + a["ref"]])
            g = (b["pos"] - a["pos"] - 1) if a["strand"] == "+" else (a["pos"] - b["pos"] - 1)
            g = (g + RL // 2) % RL - RL // 2
            if -10 <= g <= 10 and (best is None or a["nseeds"] + b["nseeds"] > best[0]): best = (a["nseeds"] + b["nseeds"], a, b)
    return (best[1], best[2]) if best else (lm, rm)

# ---------- phage Mu ----------
KMU = 25
MU_REF = next((v for k, v in REFSEQ.items() if k.startswith("MU|")), None)
MU_IDX = {}
if MU_REF:
    for i in range(len(MU_REF) - KMU + 1):
        k = MU_REF[i:i + KMU]
        MU_IDX.setdefault(k, (i, "+")); MU_IDX.setdefault(rc(k), (i, "-"))

# positions (coordonnées Mu.fa = NC_000929.1) des SNP du Mu des clones par rapport à NCBI (diffseq : 4 SNP)
MU_SNP_POS = [int(x) for x in re.split(r"[,\s]+", os.environ.get("IS91_MU_SNPS", "790,19272,29509,32841").strip()) if x]

def mu_snps(sub):
    """Allèle à chaque position SNP (ancrage sur 12 nt de part et d'autre, brin direct ou segment inversé).
    snp_diff = nombre d'allèles différents de Mu.fa (copie de référence extraite du barcode 07)."""
    if not MU_REF: return "NA", "NA"
    out, diff = [], 0
    for p in MU_SNP_POS:
        l, r = MU_REF[p - 13:p - 1], MU_REF[p:p + 12]
        base = "NA"
        if len(l) == 12 and len(r) == 12:
            lo = max(0, p - 1 - 300)
            i = sub.find(l, lo)
            while i != -1 and i < p + 300:
                if sub[i + 13:i + 25] == r: base = sub[i + 12]; break
                i = sub.find(l, i + 1)
            if base == "NA":                                   # segment inversé : rc(r) base rc(l)
                j = sub.find(rc(r))                           # le segment inversé est réfléchi : pas de fenêtre
                while j != -1:
                    if sub[j + 13:j + 25] == rc(l): base = rc(sub[j + 12]); break
                    j = sub.find(rc(r), j + 1)
        if base not in ("NA", MU_REF[p - 1]): diff += 1
        out.append(f"{p}:{base}")
    return ",".join(out), diff

def mu_scan(seq):
    """Copies de Mu dans un contig : intervalles denses en 25-mers de Mu (>= 5 kb)."""
    if not MU_IDX: return []
    get, hits = MU_IDX.get, []
    for i in range(len(seq) - KMU + 1):
        h = get(seq[i:i + KMU])
        if h: hits.append((i, h[0], h[1]))
    copies, cur = [], []
    def close():
        if not cur: return
        span = cur[-1][0] - cur[0][0] + KMU
        if span >= 5000 and len(cur) >= 0.5 * span:
            plus = sum(1 for x in cur if x[2] == "+"); minus = len(cur) - plus
            copies.append(dict(lo=cur[0][0], hi=cur[-1][0] + KMU, span=span,
                               strand="+" if plus >= minus else "-", minor=min(plus, minus),
                               cov=len({x[1] for x in cur}) / (len(MU_REF) - KMU + 1)))
    for h in hits:
        if cur and h[0] - cur[-1][0] > 300: close(); cur = []
        cur.append(h)
    close()
    return copies

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
EXPLAB = {"pOX38Km": "POX38", "chromosome": "CHR"}.get(exp_t)

ORI = list(read_fasta(os.path.join(SEQDIR, "ori.fa")).values())[0]
TER = list(read_fasta(os.path.join(SEQDIR, "ter.fa")).values())[0]
contigs = read_fasta(fasta)
rows = []

MU_BY_CONTIG, mu_rows = {}, []
for cname, seq in contigs.items():
    cps = mu_scan(seq)
    if not cps: continue
    MU_BY_CONTIG[cname] = cps
    for c in cps:
        clen, mlen = len(seq), c["hi"] - c["lo"]
        lm = rm = None; tsd_len, tsd_seq, site = "NA", "NA", ["NA", "NA", "NA"]
        if clen <= mlen * 1.2: state = "free"
        elif c["lo"] >= 500 and clen - c["hi"] >= 500:
            state = "integrated"
            Lf, Rf = seq[max(0, c["lo"] - FLANK):c["lo"]], seq[c["hi"]:c["hi"] + FLANK]
            lm, rm = reconcile(map_flank(Lf, "L"), map_flank(Rf, "R"))
            if lm and rm and (lm["label"], lm["ref"], lm["strand"]) == (rm["label"], rm["ref"], rm["strand"]):
                RL = len(REFSEQ[lm["label"] + "|" + lm["ref"]])
                g = (rm["pos"] - lm["pos"] - 1) if lm["strand"] == "+" else (lm["pos"] - rm["pos"] - 1)
                g = (g + RL // 2) % RL - RL // 2
                tsd_len = -g if -12 <= g < 0 else 0
                tsd_seq = Lf[-tsd_len:] if tsd_len else "NA"
                site = [f"{lm['label']}:{lm['ref']}", min(lm["pos"], rm["pos"]), max(lm["pos"], rm["pos"])]
            elif lm or rm:
                hm_ = lm or rm; site = [f"{hm_['label']}:{hm_['ref']}", hm_["pos"], hm_["pos"]]
        else: state = "integrated_edge"
        mu_rows.append([bc, is_type, exp_t, passage, brep, clone, cname, clen, c["lo"] + 1, c["hi"], mlen,
                        SN[c["strand"]], f"{c['cov']:.2f}", c["minor"] if c["minor"] >= 500 else 0, state,
                        fmt_map(lm), fmt_map(rm), tsd_len, tsd_seq] + site + [len(cps)] +
                       list(mu_snps(seq[c["lo"]:c["hi"]] if c["strand"] == "+" else rc(seq[c["lo"]:c["hi"]]))))
if MU_IDX and mu_rows:
    os.makedirs(os.path.join(ANA, "mu"), exist_ok=True)
    with open(os.path.join(ANA, "mu", f"{bc}.tsv"), "w") as fh:
        for r in mu_rows: fh.write("\t".join(map(str, r)) + "\n")

def mu_cols(cname, lo0, hi0):
    cps = MU_BY_CONTIG.get(cname, [])
    if not cps: return [0, "NA"]
    return [len(cps), min(max(0, c["lo"] - hi0, lo0 - c["hi"]) for c in cps)]

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
        lm, rm = reconcile(lm, rm)
        lc, rc_ = fl_class(lm, Lf), fl_class(rm, Rf)
        span = a_hi - a_lo; n = len(arr)
        ml = (lm or rm or {}).get("label")
        if ml in TARGETS: on_a = "yes" if ml == EXPLAB else "no"       # d'après la position cartographiée
        else: on_a = on
        b = base[:9] + [on_a]
        opens = []
        for u2 in units:
            if u2["strand"] != st or u2["type"] not in ("ter_only", "ori_only"): continue
            h2 = u2["ter"] or u2["ori"]
            if u2["type"] == "ter_only":     # copie partielle après le dernier ori
                g = (h2["lo"] - arr[-1]["ori"]["hi"]) if st == "+" else (arr[-1]["ori"]["lo"] - h2["hi"])
                if -10 <= g <= TANDEM_GAP: opens.append("ter_after_ori")
            else:                            # ori orphelin avant le premier ter
                g = (arr[0]["ter"]["lo"] - h2["hi"]) if st == "+" else (h2["lo"] - arr[0]["ter"]["hi"])
                if -10 <= g <= TANDEM_GAP: opens.append("ori_before_ter")
        aopen = ",".join(sorted(set(opens))) or "no"
        ref_gap, ins_ref, gseq = "NA", "NA", "NA"
        if lm and rm and lm["ref"] == rm["ref"] and lm["strand"] == rm["strand"]:
            fw = REFSEQ[lm["label"] + "|" + lm["ref"]]; RL = len(fw)
            ref_gap = (rm["pos"] - lm["pos"] - 1) if lm["strand"] == "+" else (lm["pos"] - rm["pos"] - 1)
            ref_gap = (ref_gap + RL // 2) % RL - RL // 2                # distance circulaire signée
            if 0 < ref_gap <= 6:
                ext = fw + fw
                gseq = ext[lm["pos"]:lm["pos"] + ref_gap] if lm["strand"] == "+" else \
                       rc(ext[RL + lm["pos"] - 1 - ref_gap:RL + lm["pos"] - 1])
        if lm and rm: amb = bool(lm["nloci"] > 1 and rm["nloci"] > 1)       # ambigu seulement si les deux le sont
        else: amb = bool((lm or rm or {}).get("nloci", 1) > 1)
        hm = lm or rm
        if hm: ins_ref = SN[hm["strand"]]
        extra = span - n * expected if isinstance(expected, int) else None
        if lc == "donor": oet = "yes"
        elif lc in ("unmapped", "edge") and extra is not None and extra > 100: oet = "maybe"
        else: oet = "no"
        if n > 1: ecl = "tandem"
        elif oet == "yes": ecl = "donor_both_sides" if rc_ == "donor" else "OET"
        elif oet == "maybe": ecl = "OET_candidate"
        elif lc == "target" and rc_ == "donor": ecl = "ori_side_donor"
        elif lc == "target" and rc_ == "target":
            if ref_gap != "NA" and abs(ref_gap) <= 6: ecl = "simple"
            else: ecl = "repeat_site" if amb else "rearranged_junction"
        elif "mu" in (lc, rc_): ecl = "mu_junction"
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
            rows.append(b + [ecl, oet, ai, ci, n, SN[st], s_, e_, t["s"], t["e"], o["s"], o["e"], ln,
                                expected, (ln - expected) if isinstance(expected, int) else "NA", len(inner),
                                a_lo + 1, a_hi, rc(tail), tail, lc, fmt_map(lm), rc_, fmt_map(rm), ref_gap,
                                gseq, ins_ref, aopen] + mu_cols(cname, a_lo, a_hi) + [up, dn, inner])
    for u in units:
        if u["type"] == "complete": continue
        h = u["ter"] or u["ori"]
        t, o = u["ter"], u["ori"]
        d = dict(zip(COLS, base)); lo0 = min(x["s"] for x in (t, o) if x) - 1; hi0 = max(x["e"] for x in (t, o) if x)
        d.update(event_class=u["type"], copy_index=1, n_copies=1, strand=SN[u["strand"]], ins_start=lo0 + 1, ins_end=hi0,
                 ter_start=t["s"] if t else "NA", ter_end=t["e"] if t else "NA",
                 ori_start=o["s"] if o else "NA", ori_end=o["e"] if o else "NA", expected_len=expected)
        d["mu_in_contig"], d["mu_dist"] = mu_cols(cname, lo0, hi0)
        rows.append([d.get(c, "NA") for c in COLS])

for r in rows: assert len(r) == len(COLS), (len(r), len(COLS))
def sk(r):
    v = r[16]; return (r[6], int(v) if str(v).isdigit() else 0)
for r in sorted(rows, key=sk):
    print("\t".join(map(str, r)))
