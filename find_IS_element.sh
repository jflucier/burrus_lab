#!/bin/bash

# --- Set Defaults ---
export NCORES=32
export COVERAGE=95
export PIDENT=95
export GENOMES="/fast/def-burrusvi/20260323_isoIS91/ensembl_bacteria/Release_62/fasta"
export SKIP_BLAST_STEP="false"

# --- Initialize Mandatory Variables as Empty ---
IS_FA=""
ROOT_PATH=""

# --- Parse Arguments ---
usage() {
    echo "Usage: $0 -f <is_fasta> -r <root_path> [-n ncores] [-c coverage] [-p pident] [-g genomes_path] [-s skip_blast_step]"
    exit 1
}

while getopts "f:r:n:c:g:t:o:l:s:p:" opt; do
  case $opt in
    f) IS_FA="$OPTARG" ;;
    r) ROOT_PATH="$OPTARG" ;;
    n) NCORES="$OPTARG" ;;
    c) COVERAGE="$OPTARG" ;;
    p) PIDENT="$OPTARG" ;;
    g) GENOMES="$OPTARG" ;;
    s) SKIP_BLAST_STEP="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Check Mandatory Args ---
if [[ -z "$IS_FA" || -z "$ROOT_PATH" ]]; then
    echo "Error: -f (tnpa_fasta) and -r (root_path) are mandatory."
    usage
fi

# --- Logic for dynamic naming ---
b=$(basename "${IS_FA}")
n=${b%.@(fa|txt)}
export BASE_PATH="${ROOT_PATH}/is_element_${n}_qcov${COVERAGE}_pident${PIDENT}"

# --- Export everything else ---
export NCORES COVERAGE GENOMES IS_FA
export BLASTOUT="${BASE_PATH}/blast"
export XTRACTOUT="${BASE_PATH}/extract"

mkdir -p "$BLASTOUT" "$XTRACTOUT"

echo "Running with $NCORES cores at $COVERAGE% coverage and pident $PIDENT%."
echo "FASTA: $IS_FA"
echo "Base Path: $BASE_PATH"

# Define the processing function
do_parallel_blast() {
    local fa=$1

    local b=${fa%.fa.gz}
    local n=$(basename $b)
    local tmp_fa="${b}_tmp.fa"
    local tmp_db="${b}_db"
    local out_file="${BLASTOUT}/blast.qcov${COVERAGE}/${n}.is_element_hits.txt"

    # 1. Prepare
    gunzip -c "$fa" > "$tmp_fa"
    makeblastdb -in "$tmp_fa" -dbtype nucl -out "$tmp_db" -logfile /dev/null

    # 2. BLAST (using 1 thread per blast call, parallel handles the rest)
    blastn \
      -query "${BLASTOUT}/in.fa" \
      -db "$tmp_db" \
      -num_threads 1 \
      -qcov_hsp_perc ${COVERAGE} \
      -perc_identity ${PIDENT} \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qcovs stitle" \
      -out "$out_file"

    # 3. Cleanup
    rm "$tmp_fa" "$tmp_db".n*
}
export -f do_parallel_blast

extract_and_feature() {
    # Parallel passes the whole line as $1
    local line="$1"

    # Use read to split the tab-separated line into variables
    # This must match your 24-column TSV structure exactly
    #qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qcovs stitle
    IFS=$'\t' read -r blast_fname qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qcovs stitle <<< "$line"
    local n=$(basename "$blast_fname" .dna.toplevel.is_element_hits.txt)

    # 1. Find Genome (Lowercase logic included for safety)
    local fa_gz=$(find "${GENOMES}" -path "*/*/dna/${n}*.dna.toplevel.fa.gz" | head -n 1)
    [[ -z "$fa_gz" ]] && { echo "Error: Genome $n not found" >&2; return 1; }

    # 2. Prepare Temp (Use BASHPID for thread safety)
    local tmp_fa="tmp_${n}_${BASHPID}.fa"
    gunzip -c "$fa_gz" > "$tmp_fa"
    samtools faidx "$tmp_fa"

    local max_len=$(awk -v s="${sseqid}" '$1==s {print $2}' "${tmp_fa}.fai")
    if [[ -z "$max_len" || ! "$max_len" =~ ^[0-9]+$ ]]; then
        rm -f "$tmp_fa"*
        echo "Error: chr length not defined " >&2;
        return 1
    fi

    local genomic_min genomic_max
    if [ "$sstart" -le "$send" ]; then
        genomic_min=$sstart
        genomic_max=$send
    else
        genomic_min=$send
        genomic_max=$sstart
    fi

    local exp_start=$genomic_min
    local exp_end=$genomic_max

    local rel_start rel_end strand_label
    if [[ "$sstrand" == "plus" ]]; then
        # On Plus strand: local start = (genomic start of gene - genomic start of window) + 1
        rel_start=$(( sstart - exp_start + 1 ))
        rel_end=$(( send - exp_start + 1 ))
        local strand_label="plus"
    else
        # On Minus strand (RC): local start = (genomic end of window - genomic end of gene) + 1
        # Because RC flips the sequence, the genomic END of the gene is now closer to the local START.
        rel_start=$(( exp_end - sstart + 1 ))
        rel_end=$(( exp_end - send + 1 ))
        local strand_label="minus"
    fi

    local header="${n}_${sseqid}_chrlen${max_len}_${sstart}-${send}_pident${pident}_qcov${qcovs}_strand${strand_label}_relstart${rel_start}_relend${rel_end}"
    local full_out="${XTRACTOUT}/is_element_seqs/${n}_chrlen${max_len}_seqs.fasta"

    # 4. Extract
    if [ "$strand_label" == "plus" ]; then
        samtools faidx "$tmp_fa" "${sseqid}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    else
        samtools faidx --reverse-complement "$tmp_fa" "${sseqid}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    fi

    # 5. Output Feature Row
    echo -e "${blast_fname}\t${header}\t${qseqid}\t${sseqid}\t${max_len}\t${sstart}\t${send}\t${sstrand}\t${pident}\t${qcovs}\t${rel_start}\t${rel_end}"

    # 6. Cleanup
    rm -f "$tmp_fa" "$tmp_fa".fai
}
export -f extract_and_feature

if [[ ${SKIP_BLAST_STEP} != "false" ]]; then
  echo "##### Skipping tnpA blast step #####"
else
  echo "##### finding isoform using BLAST #####"
  module purge
  ml StdEnv/2020 emboss/6.6.0
  seqret ${IS_FA} ${BLASTOUT}/in.fa

  module purge 2>/dev/null
  ml StdEnv/2020 gcc/9.3.0 blast+/2.14.0 2>/dev/null

  mkdir -p ${BLASTOUT}/blast.qcov${COVERAGE}
  find "${BLASTOUT}/blast.qcov${COVERAGE}/" -type f -name "*.txt" -delete
  total_files=$(find ${GENOMES}/ -path "*/dna/*.dna.toplevel.fa.gz" | wc -l)
  find ${GENOMES}/ -path "*/dna/*.dna.toplevel.fa.gz" | \
      pv -l -s "$total_files" | \
      parallel --jobs ${NCORES} do_parallel_blast {}

  #qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qcovs stitle
  HEADER="filename\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore\tqcovs\tstitle"
  echo -e "$HEADER" > "${BLASTOUT}/blast.qcov${COVERAGE}.tsv"
  echo -e "$HEADER" > "${BLASTOUT}/blast.qcov${COVERAGE}.filtered.tsv"

  find "${BLASTOUT}/blast.qcov${COVERAGE}" -name "*.txt" -not -empty -exec awk -v OFS="\t" '{ print FILENAME, $0 }' {} + >> "${BLASTOUT}/blast.qcov${COVERAGE}.tsv"
  tail -n +2 "${BLASTOUT}/blast.qcov${COVERAGE}.tsv" | \
      sort -t$'\t' -k1,1 -k14,14rn -k4,4rn | \
      awk -F'\t' '!seen[$1]++' >> "${BLASTOUT}/blast.qcov${COVERAGE}.filtered.tsv"
fi

#### Step 2: extract tnpA sequence with nt buffer ####
echo "##### extracting isoform sequences from blast hits #####"
module purge 2>/dev/null
ml StdEnv/2020 samtools/1.17 2>/dev/null
FILTERED_TSV=${BLASTOUT}/blast.qcov${COVERAGE}.filtered.tsv
mkdir -p ${XTRACTOUT}/is_element_seqs/
rm -f ${XTRACTOUT}/is_element_seqs/*

# "${blast_fname}\t${header}\t${qseqid}\t${sseqid}\t${max_len}\t${gstart}\t${gend}\t${gstrand}\t${pident}\t${qcov}\t${rel_start}\t${rel_end}"
echo -e "filename\theader\tis_name\tchr\tchr_len\tstart\tend\tstrand\tident\tcoverage\trel_start\trel_end" > "${XTRACTOUT}/is_element_seqs.features.tsv"
# Use Parallel to run the extraction
total_lines=$(tail -n +2 "$FILTERED_TSV" | wc -l)
tail -n +2 "$FILTERED_TSV" | \
pv -l -s "$total_lines" | \
parallel --jobs "$NCORES" extract_and_feature {} >> "${XTRACTOUT}/is_element_seqs.features.tsv"

# gen sequence tsv for db import
find "${XTRACTOUT}/is_element_seqs/" -type f -name "*_seqs.fasta" -exec cat {} + > "${XTRACTOUT}/is_element_seqs.fasta"

echo -e "tnpa_seqsig\tseq" > "${XTRACTOUT}/is_element_seqs.fastas.tsv"
awk '/^>/ { if (header) print header "\t" seq; header = substr($0,2); seq = ""; next } { seq = seq $0 } END { if (header) print header "\t" seq }' "${XTRACTOUT}/is_element_seqs.fasta" >> "${XTRACTOUT}/is_element_seqs.fastas.tsv"

echo "##### done #####"





