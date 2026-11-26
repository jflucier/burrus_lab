#!/bin/bash
#SBATCH --job-name=IS91_find
#SBATCH --account=def-burrusvi
#SBATCH --array=1-96
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=bio
#SBATCH --output=logs/find_%A_%a.out

set -euo pipefail
module load python/3

export IS91_FLANK=3000      # longueur des flancs cartographiés de part et d'autre de l'insertion

BASE=/net/nfs-bio/jbod2/def-burrusvi/20261120_IS91_sequencing
cd "$BASE"
mkdir -p analysis/is91 logs

BC=$(printf "barcode%02d" "$SLURM_ARRAY_TASK_ID")
QUERY=data/long_read/${BC}.fa
[[ -s "$QUERY" ]] || { echo "absent: $QUERY" >&2; exit 0; }

python3 analysis/find_IS91.py "$BC" "$QUERY" > analysis/is91/${BC}.tsv
