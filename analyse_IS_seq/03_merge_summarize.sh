#!/bin/bash
#SBATCH --job-name=IS91_summ
#SBATCH --account=def-burrusvi
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --partition=bio
#SBATCH --output=logs/summ_%j.out

set -euo pipefail
module load python/3
cd /net/nfs-bio/jbod2/def-burrusvi/20261120_IS91_sequencing

{ python3 analysis/find_IS91.py --header; cat analysis/is91/barcode*.tsv; } > analysis/IS91_insertions.tsv
python3 analysis/summarize_IS91.py analysis/IS91_insertions.tsv analysis/samples.tsv analysis/
