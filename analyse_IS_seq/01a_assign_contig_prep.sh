#!/bin/bash
#SBATCH --job-name=IS91_prep
#SBATCH --account=def-burrusvi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --partition=bio
#SBATCH --output=logs/prep_%j.out

set -euo pipefail
module load gcc/9.3.0 blast+/2.12.0

BASE=/net/nfs-bio/jbod2/def-burrusvi/20261120_IS91_sequencing
cd "$BASE"
mkdir -p ${BASE}/analysis/db ${BASE}/analysis/blast ${BASE}/analysis/assign logs

[[ -s sequences/pOX38Km.fa ]] || { echo "sequences/pOX38Km.fa manquant" >&2; exit 1; }

# En-têtes ASCII propres
awk '/^>/ {print ">CHR|" substr($1,2); next} {print}' \
    sequences/GCF_056501015.1_ASM5650101v1_genomic.fna > analysis/db/refs.fa
awk 'BEGIN{n=0} /^>/ {n++; print ">POX38|pOX38Km_" n; next} {print}' \
    sequences/pOX38Km.fa >> analysis/db/refs.fa

rm -f analysis/db/refs.n*
makeblastdb -in analysis/db/refs.fa -dbtype nucl -out analysis/db/refs

cat > ${BASE}/analysis/assign_contigs.py << 'EOF'
#!/usr/bin/env python3
import sys
from collections import defaultdict

bc, fasta = sys.argv[1], sys.argv[2]
MIN_FRAC, MIN_ID = 0.5, 90.0

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

iv = defaultdict(lambda: defaultdict(list))
for l in open(f"analysis/blast/{bc}.tsv"):
    f = l.rstrip().split("\t")
    if float(f[4]) < MIN_ID: continue
    qs, qe = sorted((int(f[6]), int(f[7])))
    iv[f[0]][f[2].split("|")[0]].append((qs, qe))

for q, L in qlen.items():
    cc = merged_len(iv[q]["CHR"])   / L if L else 0
    cp = merged_len(iv[q]["POX38"]) / L if L else 0
    if cc >= MIN_FRAC and cp >= MIN_FRAC: a = "AMBIGU"
    elif cp >= MIN_FRAC and cp > cc:      a = "pOX38Km"
    elif cc >= MIN_FRAC:                  a = "chromosome"
    elif cc > 0 or cp > 0:                a = "partiel_chimere"
    else:                                 a = "aucun_hit"
    print(f"{bc}\t{q}\t{L}\t{cc:.2f}\t{cp:.2f}\t{a}")
EOF