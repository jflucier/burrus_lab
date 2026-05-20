#!/bin/bash

IS_FA=""
ROOT_PATH=""
MODEL=""

# --- Parse Arguments ---
usage() {
    echo "Usage: $0 -f <is_fasta> -o <out_path> -n <model_name>"
    exit 1
}

while getopts "f:o:n:" opt; do
  case $opt in
    f) IS_FA="$OPTARG" ;;
    o) ROOT_PATH="$OPTARG" ;;
    n) MODEL="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Check Mandatory Args ---
if [[ -z "$IS_FA" || -z "$ROOT_PATH" || -z "$MODEL" ]]; then
    echo "Error: -f (fasta) and -r (root_path) and -n (model name) are mandatory."
    usage
fi

# --- Logic for dynamic naming ---
export BASE_PATH="${ROOT_PATH}/models"

echo "##### generating fasta #####"
mkdir -p ${BASE_PATH}
module purge
ml StdEnv/2020 emboss/6.6.0
seqret ${IS_FA} ${BASE_PATH}/${MODEL}.fa

echo "##### running mafft #####"
ml StdEnv/2020 mafft/7.471
mafft --auto ${BASE_PATH}/${MODEL}.fa > ${BASE_PATH}/${MODEL}_aln.fasta

# Predict consensus structure and output Stockholm format
echo "##### running RNAalifold #####"
ml viennarna/2.5.1
RNAalifold --aln-stk=${MODEL}_aln ${BASE_PATH}/${MODEL}_aln.fasta
mv ${MODEL}_aln.stk ${BASE_PATH}/
mv alirna.ps ${BASE_PATH}/${MODEL}_aln.ps
sed -n '1,/^\/\//p' ${BASE_PATH}/${MODEL}_aln.stk > ${BASE_PATH}/${MODEL}_aln.clean.stk

# 2. Add an explicit name tag to the header if it's missing
echo "##### cleaning stockholm file #####"
if ! grep -q "#=GF ID" ${BASE_PATH}/${MODEL}_aln.clean.stk; then
    sed -i "2i #=GF ID ${MODEL}" ${BASE_PATH}/${MODEL}_aln.clean.stk
fi

echo "##### running cmbuild and cmcalibrate #####"
module purge
ml  StdEnv/2023  gcc/12.3 infernal/1.1.5
cmbuild ${BASE_PATH}/${MODEL}.cm ${BASE_PATH}/${MODEL}_aln.clean.stk
cmcalibrate ${BASE_PATH}/${MODEL}.cm

echo "##### done #####"





