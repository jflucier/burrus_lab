#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 80gb
#SBATCH --time 72:00:00
#SBATCH --output=find_iso_%A.log
#SBATCH -A def-burrusvi

## to call
#sbatch find_isoIS.slurm \
#--fasta '/path/to/fa' \
#--ter_model '/path/to/model' \
#--ori_model '/path/to/model'
#--skip_blast_step 'true|false' \
#--find_orf121 'true|false'

FA=""
ORI=""
TER=""
SKIP_BLAST_STEP="false"
FIND_ORF121="false"
# FIXED: Added missing colon after z to allow arguments (z:)
TEMP=$(getopt -o f:o:t:s:z: --long fasta:,ter_model:,ori_model:,skip_blast_step:,find_orf121: -n 'find_isoIS.slurm' -- "$@")
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -f|--fasta) FA="$2" ; shift 2 ;;
        -t|--ter_model) TER="$2" ; shift 2 ;;
        -o|--ori_model) ORI="$2" ; shift 2 ;;
        -s|--skip_blast_step) SKIP_BLAST_STEP="$2" ; shift 2 ;;
        -z|--find_orf121) FIND_ORF121="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) echo "Invalid Option" ; exit 1 ;;
    esac
done
if [ -z "$FA" ]; then
    echo "Error: The -f or --fasta option is required."
    exit 1
fi

if [ -z "$TER" ]; then
    echo "Error: The -t or --ter_model option is required."
    exit 1
fi

if [ -z "$ORI" ]; then
    echo "Error: The -o or --ori_model option is required."
    exit 1
fi

echo "##### Provided PARAMETERS: ######"
echo "##### FASTA: ${FA} ######"
echo "##### TER: $TER ######"
echo "##### SKIP_BLAST_STEP: $SKIP_BLAST_STEP ######"
echo "##### FIND_ORF121: $FIND_ORF121 ######"
# Hardcoded reference directory used by your script
GENOMES_DIR="/net/nfs-bio/fast/def-burrusvi/20260323_isoIS91/ensembl_bacteria/Release_62/fasta"
MAP_FILE="/net/nfs-bio/fast/def-burrusvi/20260323_isoIS91/ensembl_bacteria/Release_62/genome_map.tsv"

echo "##### Generating Genome File Path Map... #####"
if [[ -s "$MAP_FILE" ]]; then
    echo "##### Found existing genome map at ${MAP_FILE}. Skipping generation step. #####"
else
    echo "##### Genome map not found or empty. Generating new file... #####"

    # Find the protein files once
    find "$GENOMES_DIR" -path "*/pep/*.pep.all.fa.gz" | while read -r pep_path; do
        # Extract the directory containing the pep file (e.g., .../bradyrhizobium_sp_.../pep)
        pep_dir=$(dirname "$pep_path")

        # Navigate directly to the sibling dna directory and resolve the wildcard
        # Ensembl adds a trailing underscore after the assembly version string in the DNA file name
        dna_path_array=("$pep_dir/../dna/"*.dna.toplevel.fa.gz)

        # Verify that the file exists and is not an unexpanded wildcard string
        if [[ -f "${dna_path_array[0]}" ]]; then
            echo -e "${pep_path}\t${dna_path_array[0]}" >> "$MAP_FILE"
        fi
    done
    echo "##### Genome file mapping complete: saved to ${MAP_FILE} #####"
fi

export GENOME_MAP_FILE="$MAP_FILE"


for f in 1000 2000 10000
do
  for c in 99 90 80 70 60
  do
    echo "running coverage=$c and flanking=$f on ${FA}"
    /net/nfs-bio/fast/def-burrusvi/programs/burrus_lab/find_isoIS.sh \
    -f ${FA} \
    -r $PWD -n ${SLURM_CPUS_PER_TASK} -c $c -l $f \
    -t ${TER} \
    -o ${ORI} \
    -s ${SKIP_BLAST_STEP} -z ${FIND_ORF121}
  done
done

