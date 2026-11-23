#!/bin/bash
#SBATCH --job-name=IS91_assign
#SBATCH --account=def-burrusvi
#SBATCH --array=1-96
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/assign_%A_%a.out
#SBATCH --error=logs/assign_%A_%a.err

set -euo pipefail

module load gcc/9.3.0 blast+/2.12.0 python/3

BASE=/net/nfs-bio/jbod2/def-burrusvi/20261120_IS91_sequencing
cd "$BASE"
mkdir -p analysis/db analysis/blast analysis/assign logs

BC=$(printf "barcode%02d" "$SLURM_ARRAY_TASK_ID")
QUERY=data/long_read/${BC}.fa

# Barcode absent ou vide : on sort proprement
if [[ ! -s "$QUERY" ]]; then
    echo "Fichier absent ou vide: $QUERY" >&2
    exit 0
fi

# BLAST
blastn -query "$QUERY" -db analysis/db/refs \
    -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore" \
    -evalue 1e-20 -num_threads "$SLURM_CPUS_PER_TASK" \
    > analysis/blast/${BC}.tsv

# Assignation
python3 analysis/assign_contigs.py "$BC" "$QUERY" > analysis/assign/${BC}.tsv