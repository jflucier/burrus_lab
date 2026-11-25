#!/bin/bash
#SBATCH --job-name=IS91_prep
#SBATCH --account=def-burrusvi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --partition=bio
#SBATCH --output=logs/prep_%j.out

set -euo pipefail
export LC_ALL=C
module load gcc/9.3.0 blast+/2.12.0

BASE=/net/nfs-bio/jbod2/def-burrusvi/20261120_IS91_sequencing
cd "$BASE"
mkdir -p analysis/db analysis/blast analysis/assign logs

# Références : LABEL:fichier  (l'ordre définit l'ordre des colonnes de sortie)
REFS=(
  "CHR:sequences/GCF_056501015.1_ASM5650101v1_genomic.fna"
  "POX38:sequences/pOX38Km.fa"
  "pCP54:sequences/pCP54.fa"
  "pCP55:sequences/pCP55.fa"
  "pCP7:sequences/pCP7.fa"
)

: > analysis/db/refs.fa
: > analysis/db/labels.txt

for r in "${REFS[@]}"; do
    lab=${r%%:*}
    f=${r#*:}
    [[ -s "$f" ]] || { echo "Fichier manquant ou vide: $f (exporter le .dna en FASTA ?)" >&2; exit 1; }
    # En-têtes ASCII propres : LABEL|accession (CHR) ou LABEL|LABEL_n (plasmides)
    awk -v L="$lab" 'BEGIN{n=0}
        /^>/ {n++; print ">" L "|" (L=="CHR" ? substr($1,2) : L "_" n); next}
        {print}' "$f" >> analysis/db/refs.fa
    echo "$lab" >> analysis/db/labels.txt
    echo "$lab <- $f : $(grep -c '>' "$f") record(s)"
done

# Contrôle ASCII des en-têtes
if grep '>' analysis/db/refs.fa | grep -qP '[^\x00-\x7F]'; then
    echo "Caractères non ASCII dans les en-têtes" >&2; exit 1
fi

rm -f analysis/db/refs.n*
makeblastdb -in analysis/db/refs.fa -dbtype nucl -out analysis/db/refs

cat > "${BASE}/analysis/assign_contigs.py" << 'EOF'
#!/usr/bin/env python3
import sys
from collections import defaultdict

bc, fasta = sys.argv[1], sys.argv[2]
MIN_FRAC = 0.5    # couverture min du contig pour assigner
MIN_ID   = 90.0   # identité min des alignements
TIE      = 0.10   # écart min de couverture entre 1re et 2e référence
RENAME   = {"CHR": "chromosome", "POX38": "pOX38Km"}

labels = [l.strip() for l in open("analysis/db/labels.txt") if l.strip()]

def merged_len(iv):
    iv = sorted(iv); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    if cs is not None: tot += ce - cs + 1
    return tot

# Longueurs de tous les contigs (inclut ceux sans hit)
qlen = {}; name = None
for l in open(fasta):
    if l.startswith(">"):
        name = l[1:].split()[0]; qlen[name] = 0
    else:
        qlen[name] += len(l.strip())

# hits[contig][label] = [(sid, qs, qe, ss, se, aln_len), ...]
hits = defaultdict(lambda: defaultdict(list))
for l in open(f"analysis/blast/{bc}.tsv"):
    f = l.rstrip().split("\t")
    if float(f[4]) < MIN_ID: continue
    qs, qe = sorted((int(f[6]), int(f[7])))
    ss, se = int(f[8]), int(f[9])
    sid = f[2]
    hits[f[0]][sid.split("|")[0]].append((sid, qs, qe, ss, se, int(f[5])))

def aln_info(hl):
    """Coordonnées de l'alignement sur la référence : enregistrement dominant."""
    by_sid = defaultdict(list)
    for h in hl: by_sid[h[0]].append(h)
    sid = max(by_sid, key=lambda s: sum(h[5] for h in by_sid[s]))
    hs = by_sid[sid]
    plus  = sum(h[5] for h in hs if h[3] <= h[4])
    total = sum(h[5] for h in hs)
    pf = plus / total if total else 0.0
    strand = "+" if pf >= 0.5 else "-"
    r_start = min(min(h[3], h[4]) for h in hs)
    r_end   = max(max(h[3], h[4]) for h in hs)
    c_start = min(h[1] for h in hs)
    c_end   = max(h[2] for h in hs)
    return sid.split("|", 1)[1], r_start, r_end, strand, pf, c_start, c_end, len(hs)

for q, L in qlen.items():
    iv = {lab: [(h[1], h[2]) for h in hits[q][lab]] for lab in labels}
    cov = {lab: (merged_len(iv[lab]) / L if L else 0.0) for lab in labels}
    passing = sorted(((c, lab) for lab, c in cov.items() if c >= MIN_FRAC), reverse=True)
    best = None
    if not passing:
        a = "partiel_chimere" if any(c > 0 for c in cov.values()) else "aucun_hit"
        if a == "partiel_chimere":
            best = max(cov, key=cov.get)
    elif len(passing) == 1:
        best = passing[0][1]; a = RENAME.get(best, best)
    else:
        top = [lab for c, lab in passing if passing[0][0] - c < TIE]
        best = passing[0][1]
        a = RENAME.get(top[0], top[0]) if len(top) == 1 else \
            "AMBIGU(" + "/".join(RENAME.get(t, t) for t in top) + ")"

    if best is not None and hits[q][best]:
        rid, rs, re_, st, pf, cs, ce, n = aln_info(hits[q][best])
        SN = {"+": "plus", "-": "minus"}
        extra = [RENAME.get(best, best), rid, rs, re_, SN[st], f"{pf:.2f}", cs, ce, n]
    else:
        extra = ["NA"] * 9

    covs = "\t".join(f"{cov[lab]:.2f}" for lab in labels)
    print(f"{bc}\t{q}\t{L}\t{covs}\t{a}\t" + "\t".join(map(str, extra)))
EOF