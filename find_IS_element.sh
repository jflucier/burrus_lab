#!/bin/bash

# --- Set Defaults ---
NCORES=32
COVERAGE=75
GENOMES="/fast/def-burrusvi/20260323_isoIS91/ensembl_bacteria/Release_62/fasta"
SKIP_BLAST_STEP="false"

# --- Initialize Mandatory Variables as Empty ---
IS_FA=""
ROOT_PATH=""

# --- Parse Arguments ---
usage() {
    echo "Usage: $0 -f <is_fasta> -r <root_path> [-n ncores] [-c coverage] [-g genomes_path] [-s skip_blast_step]"
    exit 1
}

while getopts "f:r:n:c:g:t:o:l:s:" opt; do
  case $opt in
    f) IS_FA="$OPTARG" ;;
    r) ROOT_PATH="$OPTARG" ;;
    n) NCORES="$OPTARG" ;;
    c) COVERAGE="$OPTARG" ;;
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
n=${b%.fa}
export BASE_PATH="${ROOT_PATH}/is_element_${n}_qcov${COVERAGE}"

# --- Export everything else ---
export NCORES COVERAGE GENOMES IS_FA
export BLASTOUT="${BASE_PATH}/blast"
export XTRACTOUT="${BASE_PATH}/extract"

mkdir -p "$BLASTOUT" "$XTRACTOUT"

echo "Running with $NCORES cores at $COVERAGE% coverage."
echo "FASTA: $IS_FA"
echo "Base Path: $BASE_PATH"

# Define the processing function
do_parallel_blast() {
    local fa=$1

    local b=${fa%.fa.gz}
    local n=$(basename $b)
    local tmp_fa="${b}_tmp.fa"
    local tmp_db="${b}_db"
    local out_file="${BLASTOUT}/${n}.is_element_hits.txt"

    # 1. Prepare
    gunzip -c "$fa" > "$tmp_fa"
    makeblastdb -in "$tmp_fa" -dbtype nucl -out "$tmp_db" -logfile /dev/null

    # 2. BLAST (using 1 thread per blast call, parallel handles the rest)
    blastn \
      -query "$IS_FA" \
      -db "$tmp_db" \
      -num_threads 1 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" | \
    awk -v OFS="\t" '{
      # 2. Extract stitle (Column 15)
      full_title = $14;
      for (i=15; i<=NF; i++) full_title = full_title " " $i;

      split(full_title, words, " ");
      loc_str = words[2];

      # 3. Parse Location (assembly:chr:start:end:strand)
      n = split(loc_str, loc, ":");
      strand = (loc[n] == "1") ? "+" : "-";
      g_end  = loc[n-1];
      g_start = loc[n-2];
      chr    = loc[n-3];
      ass    = loc[n-4];

      # 5. Output all columns
      # $0 includes original 15 columns, then we append our 8 extracted ones
      print $0, ass, chr, g_start, g_end, strand;
    }' > "$out_file"

    # 3. Cleanup
    rm "$tmp_fa" "$tmp_db".n*
}
export -f do_parallel_blast

extract_and_feature() {
    # Parallel passes the whole line as $1
    local line="$1"

    # Use read to split the tab-separated line into variables
    # This must match your 24-column TSV structure exactly
    IFS=$'\t' read -r blast_fname qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle ass chr gstart gend gstrand <<< "$line"
    local n=$(basename "$blast_fname" .is_element_hits.txt)

    # 1. Find Genome (Lowercase logic included for safety)
    local fa_gz=$(find "${GENOMES}" -path "*/dna/${n}*.dna.toplevel.fa.gz" | head -n 1)
    [[ -z "$fa_gz" ]] && { echo "Error: Genome $n not found" >&2; return 1; }

    # 2. Prepare Temp (Use BASHPID for thread safety)
    local tmp_fa="tmp_${n}_${BASHPID}.fa"
    gunzip -c "$fa_gz" > "$tmp_fa"
    samtools faidx "$tmp_fa"

    local max_len=$(awk -v s="${chr}" '$1==s {print $2}' "${tmp_fa}.fai")
    if [[ -z "$max_len" || ! "$max_len" =~ ^[0-9]+$ ]]; then
        rm -f "$tmp_fa"*
        echo "Error: chr length not defined " >&2;
        return 1
    fi

    local exp_start=$gstart
    [[ "$exp_start" -lt 1 ]] && exp_start=1

    local exp_end=$gend
    [[ "$exp_end" -gt "$max_len" ]] && exp_end="$max_len"

    local rel_start rel_end
    if [[ "$gstrand" == "+" ]]; then
        # On Plus strand: local start = (genomic start of gene - genomic start of window) + 1
        rel_start=$(( gstart - exp_start + 1 ))
        rel_end=$(( gend - exp_start + 1 ))
        local strand_label="plus"
    else
        # On Minus strand (RC): local start = (genomic end of window - genomic end of gene) + 1
        # Because RC flips the sequence, the genomic END of the gene is now closer to the local START.
        rel_start=$(( exp_end - gend + 1 ))
        rel_end=$(( exp_end - gstart + 1 ))
        local strand_label="minus"
    fi

    local header="${n}_${chr}_chrlen${max_len}_${gstart}-${gend}_flank${FLANKING_SEQ_LEN}k_pident${pident}_qcov${qcov}_strand${strand_label}_relstart${rel_start}_relend${rel_end}"
    local full_out="${XTRACTOUT}/is_element_seqs/${n}_chrlen${max_len}_seqs.fasta"

    # 4. Extract
    if [ "$strand" == "plus" ]; then
        samtools faidx "$tmp_fa" "${chr}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    else
        samtools faidx --reverse-complement "$tmp_fa" "${chr}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    fi

    # 5. Output Feature Row
    echo -e "${blast_fname}\t${header}\t${qid}\t${assembly}\t${chr}\t${max_len}\t${gstart}\t${gend}\t${gstrand}\t${pident}\t${qcov}\t${rel_start}\t${rel_end}"

    # 6. Cleanup
    rm -f "$tmp_fa" "$tmp_fa".fai
}
export -f extract_and_feature

if [[ ${SKIP_BLAST_STEP} != "false" ]]; then
  echo "##### Skipping tnpA blast step #####"
else
  echo "##### finding isoform using BLAST #####"
  module purge 2>/dev/null
  ml StdEnv/2020 gcc/9.3.0 blast+/2.14.0 2>/dev/null

  mkdir -p ${BLASTOUT}/blast.qcov${COVERAGE}
  find "${BLASTOUT}/blast.qcov${COVERAGE}/" -type f -name "*.txt" -delete
  total_files=$(find ${GENOMES}/ -path "*/pep/*.pep.all.fa.gz" | wc -l)
  find ${GENOMES}/ -path "*/pep/*.pep.all.fa.gz" | \
      pv -l -s "$total_files" | \
      parallel --jobs ${NCORES} do_parallel_blast {}

  #qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle ass, chr, g_start, g_end, strand
  HEADER="filename\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\tstitle\tassembly\tchr\tstart\tend\tstrand"
  echo -e "$HEADER" > "${BASE_PATH}/blast.qcov${COVERAGE}.tsv"
  echo -e "$HEADER" > "${BASE_PATH}/blast.qcov${COVERAGE}.filtered.tsv"

  find "${BLASTOUT}" -name "*.txt" -not -empty -exec awk -v OFS="\t" '{ print FILENAME, $0 }' {} + >> "${BASE_PATH}/blast.qcov${COVERAGE}.tsv"
  tail -n +2 "${BASE_PATH}/blast.qcov${COVERAGE}.tsv" | \
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

echo -e "filename\ttnpa_seqsig\ttnpA_hit\tassembly\tchr\tchr_len\tstart\tend\tstrand\tident\tcoverage\tgene_name\tgene_desc\trel_start\trel_end" > "${XTRACTOUT}/is_element_seqs.features.tsv"
# Use Parallel to run the extraction
total_lines=$(tail -n +2 "$FILTERED_TSV" | wc -l)
tail -n +2 "$FILTERED_TSV" | \
pv -l -s "$total_lines" | \
parallel --jobs "$NCORES" extract_and_feature {} >> "${XTRACTOUT}/is_element_seqs.features.tsv"

# gen sequence tsv for db import
echo -e "tnpa_seqsig\tseq" > ${XTRACTOUT}/is_element_seqs.fastas.tsv
awk '/^>/ { if (header) print header "\t" seq; header = substr($0,2); seq = ""; next } { seq = seq $0 } END { if (header) print header "\t" seq }' ${XTRACTOUT}/is_element_seqs/*_seqs.fasta >> ${XTRACTOUT}/is_element_seqs.fastas.tsv

echo "##### done #####"





