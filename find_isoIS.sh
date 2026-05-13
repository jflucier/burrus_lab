#!/bin/bash

# --- Set Defaults ---
NCORES=32
COVERAGE=75
GENOMES="/fast/def-burrusvi/20260323_isoIS91/ensembl_bacteria/Release_62/fasta"
SKIP_BLAST_STEP="false"

# --- Initialize Mandatory Variables as Empty ---
TNPA_FA=""
ROOT_PATH=""

# --- Parse Arguments ---
usage() {
    echo "Usage: $0 -f <tnpa_fasta> -r <root_path> [-n ncores] [-c coverage] [-g genomes_path] [-t teris_model] [-o oriis_model] [-s skip_blast_step]"
    exit 1
}

while getopts "f:r:n:c:g:t:o:l:s:" opt; do
  case $opt in
    f) TNPA_FA="$OPTARG" ;;
    r) ROOT_PATH="$OPTARG" ;;
    n) NCORES="$OPTARG" ;;
    c) COVERAGE="$OPTARG" ;;
    g) GENOMES="$OPTARG" ;;
    t) TERIS_MODEL="$OPTARG" ;;
    o) ORIIS_MODEL="$OPTARG" ;;
    l) FLANKING_SEQ_LEN="$OPTARG" ;;
    s) SKIP_BLAST_STEP="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Check Mandatory Args ---
if [[ -z "$TNPA_FA" || -z "$ROOT_PATH" ]]; then
    echo "Error: -f (tnpa_fasta) and -r (root_path) are mandatory."
    usage
fi

# --- Logic for dynamic naming ---
b=$(basename "${TNPA_FA}")
n=${b%.fa}
export BASE_PATH="${ROOT_PATH}/${n}_qcov${COVERAGE}"

# --- Set Model Defaults (if not provided via flags) ---
export TERIS_MODEL=${TERIS_MODEL:-"/fast/def-burrusvi/20260323_isoIS91/models/isoforms_IS91_terIS_nt_aln.cm"}
export ORIIS_MODEL=${ORIIS_MODEL:-"/fast/def-burrusvi/20260323_isoIS91/models/isoforms_IS91_oriIS_nt_aln.cm"}
export ORF121_BLASTDB=${ORF121_BLASTDB:-"/fast/def-burrusvi/20260323_isoIS91/IS_fasta/isoforms_IS91_orf121_aa.fa"}
export FLANKING_SEQ_LEN=${FLANKING_SEQ_LEN:-"2000"}


# --- Export everything else ---
export NCORES COVERAGE GENOMES TNPA_FA
export TNPAOUT="${BASE_PATH}/tnpA_search"
export XTRACTOUT="${BASE_PATH}/tnpA_seqs_${FLANKING_SEQ_LEN}"
export TERISOUT="${BASE_PATH}/terIS_search_${FLANKING_SEQ_LEN}"
export ORIISOUT="${BASE_PATH}/oriIS_search_${FLANKING_SEQ_LEN}"
export ORF121OUT="${BASE_PATH}/orf121_search_${FLANKING_SEQ_LEN}"

mkdir -p "$TNPAOUT" "$TERISOUT" "$ORIISOUT" "$ORF121OUT"

echo "Running with $NCORES cores at $COVERAGE% coverage."
echo "FASTA: $TNPA_FA"
echo "Base Path: $BASE_PATH"
echo "Flanking tnpA seq length: $FLANKING_SEQ_LEN"

# Define the processing function
#do_parallel_blast() {
#    local fa=$1
#
#    local b=${fa%.fa.gz}
#    local n=$(basename $b)
#    local tmp_fa="${b}_tmp.fa"
#    local tmp_db="${b}_db"
#    local out_file="${TNPAOUT}/blast.qcov${COVERAGE}/${n}_tnpa_hits.txt"
#
#    # 1. Prepare
#    gunzip -c "$fa" > "$tmp_fa"
#    makeblastdb -in "$tmp_fa" -dbtype nucl -out "$tmp_db" -logfile /dev/null
#
#    # 2. BLAST (using 1 thread per blast call, parallel handles the rest)
#    tblastn \
#      -query "$TNPA_FA" \
#      -db "$tmp_db" \
#      -num_threads 1 \
#      -db_gencode 11 \
#      -qcov_hsp_perc ${COVERAGE} \
#      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
#      -out "$out_file"
#
#    # 3. Cleanup
#    rm "$tmp_fa" "$tmp_db".n*
#}
#export -f do_parallel_blast

do_parallel_blast_pep() {
    local fa=$1
    local b=${fa%.all.fa.gz}
    local n=$(basename "$b")
    local tmp_fa="${b}_tmp.fa"
    local tmp_db="${b}_db"
    local out_file="${TNPAOUT}/blast.qcov${COVERAGE}/${n}_tnpa_hits.txt"

    # 1. Prepare (Proteins don't need translation tables)
    gunzip -c "$fa" > "$tmp_fa"
    makeblastdb -in "$tmp_fa" -dbtype prot -out "$tmp_db" -parse_seqids -logfile /dev/null

    # 2. BLASTP (Query: AA, DB: AA)
    # We include 'slen' to calculate subject coverage later
    blastp \
    -query "$TNPA_FA" \
    -db "$tmp_db" \
    -num_threads 1 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen gaps stitle" | \
    awk -v file_name="$fa" -v min_cov="$COVERAGE" -v OFS="\t" '{
        true_length = $4 - $14;
        cov = (true_length / $13) * 100;
        if (cov >= min_cov) {
           # 2. Extract stitle (Column 15)
            full_title = $15;
            for (i=16; i<=NF; i++) full_title = full_title " " $i;

            split(full_title, words, " ");
            loc_str = words[2];

            # 3. Parse Location (assembly:chr:start:end:strand)
            n = split(loc_str, loc, ":");
            strand = (loc[n] == "1") ? "+" : "-";
            g_end  = loc[n-1];
            g_start = loc[n-2];
            chr    = loc[n-3];
            ass    = loc[n-4];

            # 4. Extract Gene Name and Description correctly
            # Note the use of g[1] and d[1] for the captured groups
            gene_name = "N/A";
            if (match(full_title, /gene:([^ ]+)/, g)) { gene_name = g[1]; }

            gene_desc = "N/A";
            if (match(full_title, /description:(.+)$/, d)) { gene_desc = d[1]; }

            # 5. Output all columns
            # $0 includes original 15 columns, then we append our 8 extracted ones
            print file_name, $0, sprintf("%.4f", cov), ass, chr, g_start, g_end, strand, gene_name, gene_desc;

        }
    }' > "$out_file"

    # 3. Cleanup
    rm "$tmp_fa" "$tmp_db".p*
}
export -f do_parallel_blast_pep

extract_and_feature() {
    # Parallel passes the whole line as $1
    local line="$1"

    # Use read to split the tab-separated line into variables
    # This must match your 24-column TSV structure exactly
    IFS=$'\t' read -r fa_gz blast_fname qid sid pident len mis gap qstart qend sstart send eval bit slen gaps stitle qcov assembly chr gstart gend gstrand gene_name gene_desc <<< "$line"
    local n=$(basename "$blast_fname" .pep_tnpa_hits.txt)

    # 1. Find Genome (Lowercase logic included for safety)
#    local fa_gz=$(find "${GENOMES}" -path "*/dna/${n}*.dna.toplevel.fa.gz" | head -n 1)
    [[ -z "$fa_gz" ]] && { echo "Error: Genome $n not found" >&2; return 1; }

    # 2. Prepare Temp (Use BASHPID for thread safety)
    local tmp_fa="tmp_${n}_${BASHPID}.fa"
    gunzip -c "$fa_gz" > "$tmp_fa"
    samtools faidx "$tmp_fa"

    # 3. Coordinate Logic
    local strand="plus"
    [[ "$gstrand" == "-" ]] && strand="minus"

    # Get max length from index
    local max_len=$(awk -v s="${chr}" '$1==s {print $2}' "${tmp_fa}.fai")

    # Prevent math errors if chr not found in this specific genome
    if [[ -z "$max_len" || ! "$max_len" =~ ^[0-9]+$ ]]; then
         echo "Error: Chr $chr not found in $n" >&2
         rm -f "$tmp_fa"*
         return 1
    fi

    local exp_start=$(( gstart - FLANKING_SEQ_LEN ))
    local rel_start=$FLANKING_SEQ_LEN
    if [ "$exp_start" -lt 1 ]; then
        exp_start=1
        rel_start=${gstart}
    fi

    local exp_end=$(( gend + FLANKING_SEQ_LEN ))
    local rel_end=$(( rel_start + gend - gstart ))
    if [ "$exp_end" -gt "$max_len" ]; then
        exp_end=$max_len
    fi

    local header="${n}_${chr}_chrlen${max_len}_${gstart}-${gend}_flank${FLANKING_SEQ_LEN}k_pident${pident}_qcov${qcov}_strand${strand}_relstart${rel_start}_relend${rel_end}"
    local full_out="${XTRACTOUT}/out/${n}_chrlen${max_len}_seqs.fasta"

    # 4. Extract
    if [ "$strand" == "plus" ]; then
        samtools faidx "$tmp_fa" "${chr}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    else
        samtools faidx --reverse-complement "$tmp_fa" "${chr}:${exp_start}-${exp_end}" | sed "1s/.*/\>${header}/" >> "$full_out"
    fi

    # 5. Output Feature Row
    echo -e "${blast_fname}\t${header}\t${qid}\t${assembly}\t${chr}\t${max_len}\t${gstart}\t${gend}\t${gstrand}\t${pident}\t${qcov}\t${gene_name}\t${gene_desc}\t${rel_start}\t${rel_end}"

    # 6. Cleanup
    rm -f "$tmp_fa" "$tmp_fa".fai
}
export -f extract_and_feature

do_teris_search() {
    local fa=$1
    local n=$(basename "$fa" .fasta)
    local out_file="${TERISOUT}/cmsearch/${n}_structural_hits"

    # 1. Extract max_len directly from the FASTA filename string using Bash regex
    if [[ "$n" =~ _chrlen([0-9]+) ]]; then
        local max_len="${BASH_REMATCH[1]}"
    else
        echo "Error: Could not parse max_len from filename: $fa" >&2
        return 1
    fi

    # 2. Convert base pairs to Megabases (Mb) and divide by 2 for --toponly
    local half_contig_mb=$(echo "scale=6; ($max_len / 1000000) / 2" | bc)

    # cmsearch is the tool for searching one CM against many sequences
    #--max --toponly -Z <half_of_contig_Mb>
    cmsearch \
      --cpu 1 \
      --max \
      --toponly \
      -Z "${half_contig_mb}" \
      -o ${out_file}.out \
      --tblout ${out_file}.tbl --fmt 3 \
      ${TERIS_MODEL} "$fa" > /dev/null

    # 2. Check if hits were found (check .tbl file)
    # If no hits, delete all associated files and return
    if ! grep -qv '^#' "${out_file}.tbl"; then
        rm -f "${out_file}.tbl" "${out_file}.out"
        return
    fi

    perl -e '
    open(my $FH, "<'${out_file}'.tbl ");
    my @lines = <$FH>;
    chomp(@lines);
    foreach my $l (@lines){
      if($l !~ /^\#/){
        my @fields = split(/\s+/,$l);
        my $str = join("\t",@fields);
        print "$str\n";
      }
    }
    ' > ${out_file}.tsv
}
export -f do_teris_search

do_oriis_search() {
    local fa=$1
    local n=$(basename "$fa" .fasta)
    local out_file="${ORIISOUT}/cmsearch/${n}_structural_hits"

        # 1. Extract max_len directly from the FASTA filename string using Bash regex
    if [[ "$n" =~ _chrlen([0-9]+) ]]; then
        local max_len="${BASH_REMATCH[1]}"
    else
        echo "Error: Could not parse max_len from filename: $fa" >&2
        return 1
    fi

    # 2. Convert base pairs to Megabases (Mb) and divide by 2 for --toponly
    local half_contig_mb=$(echo "scale=6; ($max_len / 1000000) / 2" | bc)

    # cmsearch is the tool for searching one CM against many sequences
    cmsearch \
      --cpu 1 \
      --max \
      --toponly \
      -Z "${half_contig_mb}" \
      -o ${out_file}.out \
      --tblout ${out_file}.tbl --fmt 3 \
      ${ORIIS_MODEL} "$fa" > /dev/null

    # 2. Check if hits were found (check .tbl file)
    # If no hits, delete all associated files and return
    if ! grep -qv '^#' "${out_file}.tbl"; then
        rm -f "${out_file}.tbl" "${out_file}.out"
        return
    fi

    perl -e '
    open(my $FH, "<'${out_file}'.tbl ");
    my @lines = <$FH>;
    chomp(@lines);
    foreach my $l (@lines){
      if($l !~ /^\#/){
        my @fields = split(/\s+/,$l);
        my $str = join("\t",@fields);
        print "$str\n";
      }
    }
    ' > ${out_file}.tsv
}
export -f do_oriis_search

do_translation() {
    local fa=$1
    local n=$(basename "$fa" .fasta)
    local out_file="${ORF121OUT}/translations/${n}"
    local blastout_file="${ORF121OUT}/translations_best/${n}"

    getorf -table 11 -minsize 50 -find 0 "$fa" ${out_file}.fa 2>/dev/null
    perl -e '
        $/ = ">"; # Set record separator to >
        while (<>) {
          chomp;
          next if !$_;
          my ($header, @seq_lines) = split(/\n/);
          next if $header =~ /\(REVERSE SENSE\)/;

          my ($rel_start) = $header =~ /_relstart(\d+)_relend/;
          die "Error: rel_start not found in filename $n_val" unless defined $rel_start;

          if ($header =~ /\[\s*(\d+)\s*-\s*(\d+)\s*\]/) {
            my $orig_start = $1;
            my $orig_end   = $2;

            if ($orig_end < $rel_start) {
              my $seq = join("", @seq_lines);

              if ($seq =~ /M/) {
                my $offset = index($seq, "M"); # 0-based index of first M
                # 3. Update the sequence
                $seq = substr($seq, $offset);
                if (length($seq) >= 30) {
                  my $new_start = $orig_start + ($offset * 3);
                  $header =~ s/\[\s*\d+\s*-/\[$new_start -/;

                  print ">$header\n$seq\n";
                }
              }
            }
          }
        }
    ' ${out_file}.fa > ${out_file}.clean.fa

    perl -e '
        $/ = ">";
        while (<>) {
            chomp;
            next unless $_;
            my ($full_header, @lines) = split /\n/;
            my $seq = join "", @lines;

            if ($full_header =~ /^(\S+)\s+\[(\d+)\s+-\s+(\d+)\]/) {
                my $name = $1;
                my $orig_start = $2;
                my $orig_end = $3;

                # Remove the trailing underscore and digits
                my $orig_name = $name;
                $name =~ s/_\d+$//;

                my $len = length($seq);
                print "$orig_name\t$name\t$orig_start\t$orig_end\t$len\t$seq\n";
            }
        }
    ' "${out_file}.clean.fa" | sort -t$'\t' -k1,1 -k4,4rn > "${out_file}.tsv"

    ## need to check if "${out_file}.clean.fa" is empty before running blast
    if [[ -s "${out_file}.clean.fa" ]]; then
      blastp -query "${out_file}.clean.fa" \
      -db ${ORF121_BLASTDB} \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq" \
      -max_target_seqs 10 | \
      awk -F'\t' 'BEGIN { OFS="\t" }
      {
        # 1. Prepare names
        orig_name = $1;
        trimmed_id = $1;
        sub(/_[0-9]+$/, "", trimmed_id);

        # 2. Calculate coverage: (Alignment Length / Subject Length) * 100
        current_cov = ($4 / $13) * 100;

        # 3. Use trimmed_id as the key to find the best hit across all fragments
        if (current_cov > max_cov[trimmed_id]) {
            max_cov[trimmed_id] = current_cov;

            # 4. Reconstruct line: $1=Orig, $2=Trimmed, $3...rest of line
            # We shift original columns 2..NF to 3..NF+1
            line = orig_name "\t" trimmed_id;
            for (i=2; i<=NF; i++) {
                line = line "\t" $i;
            }

            best_row[trimmed_id] = line;
        }
      }
      END {
        for (q in best_row) {
            print best_row[q]
        }
      }' > "${blastout_file}.blast_hits.tsv"
    fi
}
export -f do_translation

##### Step 1: find tnpA with blast #####

### old method blasting prot against genome. Changed to remove truncated hits
#mkdir -p ${TNPAOUT}/blast.qcov${COVERAGE}
#total_files=$(find ${GENOMES}/ -path "*/dna/*.dna.toplevel.fa.gz" | wc -l)
#find ${GENOMES}/ -path "*/dna/*.dna.toplevel.fa.gz" | \
#pv -l -s "$total_files" | \
#parallel --jobs ${NCORES} do_parallel_blast {}
#
#echo -e "filename\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs" >${TNPAOUT}/blast.qcov${COVERAGE}.tsv
#find ${TNPAOUT}/blast.qcov${COVERAGE} -name "*.txt" -not -empty -exec awk '{OFS="\t"; print FILENAME, $0}' {} + >> ${TNPAOUT}/blast.qcov${COVERAGE}.tsv
#
## filter duplicate hits
#echo -e "filename\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs" >${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv
#tail -n +2 ${TNPAOUT}/blast.qcov${COVERAGE}.tsv | \
#sort -t$'\t' -k1,1 -k4,4rn -k14,14rn | \
#awk -F'\t' '!seen[$1]++' >> ${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv

if [[ ${SKIP_BLAST_STEP} != "false" ]]; then
  echo "##### Skipping tnpA blast step #####"
else
  echo "##### finding isoform using BLAST #####"
  module purge 2>/dev/null
  ml StdEnv/2020 gcc/9.3.0 blast+/2.14.0 2>/dev/null

  mkdir -p ${TNPAOUT}/blast.qcov${COVERAGE}
  find "${TNPAOUT}/blast.qcov${COVERAGE}/" -type f -name "*.txt" -delete
  total_files=$(find ${GENOMES}/ -path "*/pep/*.pep.all.fa.gz" | wc -l)
  find ${GENOMES}/ -path "*/pep/*.pep.all.fa.gz" | \
      pv -l -s "$total_files" | \
      parallel --jobs ${NCORES} do_parallel_blast_pep {}

  HEADER="filename\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen\tgaps\tstitle\tsubject_cov\tassembly\tchr\tstart\tend\tstrand\tgene_name\tgene_desc"
  echo -e "$HEADER" > "${TNPAOUT}/blast.qcov${COVERAGE}.tsv"
  echo -e "$HEADER" > "${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv"

  find "${TNPAOUT}/blast.qcov${COVERAGE}" -name "*.txt" -not -empty -exec awk -v OFS="\t" '{ print FILENAME, $0 }' {} + >> "${TNPAOUT}/blast.qcov${COVERAGE}.tsv"
  tail -n +2 "${TNPAOUT}/blast.qcov${COVERAGE}.tsv" | \
      sort -t$'\t' -k2,2 -k17,17rn -k4,4rn | \
      awk -F'\t' '!seen[$1]++' >> "${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv"
fi

#### Step 2: extract tnpA sequence with nt buffer ####
echo "##### extracting isoform sequences from blast hits #####"
module purge 2>/dev/null
ml StdEnv/2020 samtools/1.17 2>/dev/null
FILTERED_TSV=${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv
mkdir -p ${XTRACTOUT}/out/
rm -f ${XTRACTOUT}/out/*

echo -e "filename\ttnpa_seqsig\ttnpA_hit\tassembly\tchr\tchr_len\tstart\tend\tstrand\tident\tcoverage\tgene_name\tgene_desc\trel_start\trel_end" > "${XTRACTOUT}/tnpA_seqs.features.tsv"
# Use Parallel to run the extraction
total_lines=$(tail -n +2 "$FILTERED_TSV" | wc -l)
tail -n +2 "$FILTERED_TSV" | \
pv -l -s "$total_lines" | \
parallel --jobs "$NCORES" extract_and_feature {} >> "${XTRACTOUT}/tnpA_seqs.features.tsv"

# gen sequence tsv for db import
echo -e "tnpa_seqsig\tseq" > ${XTRACTOUT}/tnpA_seqs.fastas.tsv
awk '/^>/ { if (header) print header "\t" seq; header = substr($0,2); seq = ""; next } { seq = seq $0 } END { if (header) print header "\t" seq }' ${XTRACTOUT}/out/*_seqs.fasta >> ${XTRACTOUT}/tnpA_seqs.fastas.tsv

### RUN ONCE: lets find terIS et oriIS in extracted seqs ####
#module purge
#ml StdEnv/2020 mafft/7.471
#mkdir terIS_search
#mafft --auto IS_fasta/isoforms_IS91_terIS_nt.fa > terIS_search/isoforms_IS91_terIS_nt_aln.fasta
#mkdir oriIS_search
#mafft --auto IS_fasta/isoforms_IS91_oriIS_nt.fa > oriIS_search/isoforms_IS91_oriIS_nt_aln.fasta
#
## Predict consensus structure and output Stockholm format
#ml viennarna/2.5.1
#RNAalifold --aln-stk=isoforms_IS91_terIS_nt_aln terIS_search/isoforms_IS91_terIS_nt_aln.fasta
#mv RNAalifold_results.stk terIS_search/isoforms_IS91_terIS_nt_aln.stk
#sed -n '1,/^\/\//p' terIS_search/isoforms_IS91_terIS_nt_aln.stk > terIS_search/isoforms_IS91_terIS_nt_aln.clean.stk
#RNAalifold --aln-stk=isoforms_IS91_oriIS_nt_aln oriIS_search/isoforms_IS91_oriIS_nt_aln.fasta
#mv isoforms_IS91_oriIS_nt_aln.stk oriIS_search/isoforms_IS91_oriIS_nt_aln.stk
#sed -n '1,/^\/\//p' oriIS_search/isoforms_IS91_oriIS_nt_aln.stk > oriIS_search/isoforms_IS91_oriIS_nt_aln.clean.stk
#
## 2. Add an explicit name tag to the header if it's missing
#if ! grep -q "#=GF ID" terIS_search/isoforms_IS91_terIS_nt_aln.clean.stk; then
#    sed -i '2i #=GF ID isoforms_IS91_terIS_nt' terIS_search/isoforms_IS91_terIS_nt_aln.clean.stk
#fi
#if ! grep -q "#=GF ID" oriIS_search/isoforms_IS91_oriIS_nt_aln.clean.stk; then
#    sed -i '2i #=GF ID isoforms_IS91_oriIS_nt' oriIS_search/isoforms_IS91_oriIS_nt_aln.clean.stk
#fi
#
#module purge
#ml  StdEnv/2023  gcc/12.3 infernal/1.1.5
## 1. Build the Covariance Model
#cmbuild terIS_search/isoforms_IS91_terIS_nt_aln.cm terIS_search/isoforms_IS91_terIS_nt_aln.clean.stk
#cmbuild oriIS_search/isoforms_IS91_oriIS_nt_aln.cm oriIS_search/isoforms_IS91_oriIS_nt_aln.clean.stk
## 2. Calibrate
#cmcalibrate terIS_search/isoforms_IS91_terIS_nt_aln.cm
#cmcalibrate oriIS_search/isoforms_IS91_oriIS_nt_aln.cm

##### Step 3: finding terIS #####
echo "##### finding terIS #####"
echo "Using Model: ${TERIS_MODEL}"
module purge 2>/dev/null
ml  StdEnv/2023  gcc/12.3 infernal/1.1.5 2>/dev/null
mkdir -p ${TERISOUT}/cmsearch
rm -f ${TERISOUT}/cmsearch/*

total_files=$(find ${XTRACTOUT}/out/ -name "*.fasta" | wc -l)
find ${XTRACTOUT}/out/ -name "*.fasta" | \
pv -l -s "$total_files" | \
parallel --jobs ${NCORES} --progress do_teris_search {}

echo -e "target_name\taccession\tquery_name\tquery_accession\tmodel_type\tmodel_from_coord\tmodel_to_coord\ttarget_from_coord\ttarget_to_coord\tstrand\ttrunc\tpass\tgc\tbias\tscore\tevalue\tinc\tmdl_len\tseq_len\tdescription" > ${TERISOUT}/all.teris.tsv
cat ${TERISOUT}/cmsearch/*.tsv >> ${TERISOUT}/all.teris.tsv

##### Step 4: find oriIS #####
echo "##### finding oriIS #####"
echo "Using Model: ${ORIIS_MODEL}"
mkdir -p ${ORIISOUT}/cmsearch
rm -f ${ORIISOUT}/cmsearch/*

total_files=$(find ${XTRACTOUT}/out/ -name "*.fasta" | wc -l)
find ${XTRACTOUT}/out/ -name "*.fasta" | \
pv -l -s "$total_files" | \
parallel --jobs ${NCORES} --progress do_oriis_search {}

echo -e "target_name\taccession\tquery_name\tquery_accession\tmodel_type\tmodel_from_coord\tmodel_to_coord\ttarget_from_coord\ttarget_to_coord\tstrand\ttrunc\tpass\tgc\tbias\tscore\tevalue\tinc\tmdl_len\tseq_len\tdescription" > ${ORIISOUT}/all.oriis.tsv
cat ${ORIISOUT}/cmsearch/*.tsv >> ${ORIISOUT}/all.oriis.tsv

##### Step 5: find orf121 #####
echo "##### finding orf121 #####"
module purge 2>/dev/null
ml StdEnv/2020 gcc/9.3.0 blast+/2.14.0 emboss/6.6.0 2>/dev/null

# RUN ONCE: make blast db of orf121 for best hit filtering
#echo "##### finding orf121 using BLAST #####"
#module purge 2>/dev/null
#ml StdEnv/2020 gcc/9.3.0 blast+/2.14.0 2>/dev/null
#makeblastdb -in "IS_fasta/isoforms_IS91_orf121_aa.fa" -dbtype prot -logfile /dev/null

mkdir -p ${ORF121OUT}/translations ${ORF121OUT}/translations_best
rm -f ${ORF121OUT}/translations/* ${ORF121OUT}/translations_best/*

total_files=$(find ${XTRACTOUT}/out/ -name "*.fasta" | wc -l)
find ${XTRACTOUT}/out/ -name "*.fasta" | \
pv -l -s "$total_files" | \
parallel --jobs ${NCORES} --progress do_translation {}

echo -e "orf121_name\ttarget_name\tstart\tend\tlength\tsequence" > ${ORF121OUT}/all.orf121.tsv
cat ${ORF121OUT}/translations/*.tsv >> ${ORF121OUT}/all.orf121.tsv

echo -e "orf121_name\ttarget_name\torf121_hit\tpident\ttarget_length\tmismatch\tgapopen\ttarget_start\ttarget_end\torf121_hit_start\torf121_hit_end\tevalue\tbitscore\torf121_hit_length\ttarget_align\torf121_hit_align" > ${ORF121OUT}/best_hits.orf121.tsv
cat ${ORF121OUT}/translations_best/*.tsv >> ${ORF121OUT}/best_hits.orf121.tsv

##### Step 6: generate signature report #####
# import blast results
sqlitedb="${BASE_PATH}/${n}_qcov${COVERAGE}_${FLANKING_SEQ_LEN}.db"
rm -f ${sqlitedb}
#sqlite3 ${sqlitedb} <<EOF
#.mode tabs
#.import "${TNPAOUT}/blast.qcov${COVERAGE}.filtered.tsv" tnpA_blast
#.quit
#EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import "${XTRACTOUT}/tnpA_seqs.features.tsv" tnpA_features
.quit
EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import ${XTRACTOUT}/tnpA_seqs.fastas.tsv tnpA_seqs
.quit
EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import "${TERISOUT}/all.teris.tsv" terIS
.quit
EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import "${ORIISOUT}/all.oriis.tsv" oriIS
.quit
EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import "${ORF121OUT}/all.orf121.tsv" orf121
.quit
EOF

sqlite3 ${sqlitedb} <<EOF
.mode tabs
.import "${ORF121OUT}/best_hits.orf121.tsv" orf121_best
.quit
EOF

sqlite3 ${sqlitedb} '.headers on' '.separator "\t"' "
SELECT
  f.tnpa_seqsig,
  f.tnpA_hit,
  f.assembly tnpA_hit_assembly,
  f.chr tnpA_hit_chr,
  f.start tnpA_hit_start,
  f.end tnpA_hit_end,
  f.strand tnpA_hit_strand,
  f.ident tnpA_hit_identity,
  f.coverage tnpA_hit_coverage,
  f.gene_name tnpA_hit_genename,
  f.gene_desc tnpA_hit_genedesc,
  t.target_from_coord terIS_start,
  t.target_to_coord terIS_end,
  orf.start orf121_pred_start,
  orf.end orf121_pred_end,
  orf.length orf121_pred_length,
  f.rel_start tnpA_rel_start,
  f.rel_end tnpA_rel_end,
  o.target_from_coord oriIS_start,
  o.target_to_coord oriIS_end,
  orf.sequence orf121_pred_sequence,
  b.pident orf121_hit_pident,
  b.mismatch orf121_hit_mismatch,
  b.gapopen orf121_hit_gapopen,
  b.evalue orf121_hit_evalue,
  b.bitscore orf121_hit_bitscore,
  substr(s.seq, t.target_from_coord, (t.target_to_coord - t.target_from_coord + 1)) as terIS_seq,
  substr(
      s.seq,
      orf.start,
      length(orf.sequence) * 3
  ) AS orf121_nt_seq,
  substr(s.seq,
      CASE
          WHEN f.strand = '+' THEN
              (f.start - CASE WHEN f.start - ${FLANKING_SEQ_LEN} < 1 THEN 1 ELSE f.start - ${FLANKING_SEQ_LEN} END + 1)
          WHEN f.strand = '-' THEN
              (CASE WHEN (f.end - f.start + ${FLANKING_SEQ_LEN} + ${FLANKING_SEQ_LEN} + 1) = LENGTH(s.seq) THEN ${FLANKING_SEQ_LEN} ELSE (LENGTH(s.seq) - f.end) END + 1)
      END,
      (f.end - f.start + 1)
  ) AS tnpA_nt_seq,
  substr(s.seq, o.target_from_coord, (o.target_to_coord - o.target_from_coord + 1)) as oriIS_seq,
  substr(s.seq, t.target_from_coord, (o.target_to_coord - t.target_from_coord + 1)) AS full_seq
FROM tnpA_features f
  JOIN terIS t ON t.target_name = f.tnpa_seqsig
      AND t.strand = '+'
      AND CAST(t.target_from_coord AS INT) <= CAST(f.rel_start AS INT) + 10
  JOIN oriIS o ON o.target_name = f.tnpa_seqsig
      AND o.strand = '+'
      AND CAST(o.target_from_coord AS INT) >= CAST(f.rel_end AS INT) - 10
  JOIN orf121_best b on b.target_name=f.tnpa_seqsig
  JOIN orf121 orf on orf.orf121_name=b.orf121_name
  JOIN tnpA_seqs s on s.tnpa_seqsig=f.tnpa_seqsig
WHERE
  -- (t.target_from_coord IS NULL OR (CAST(f.rel_start AS INT) >= CAST(t.target_from_coord AS INT) - 10 AND CAST(orf.end AS INT) > CAST(t.target_from_coord AS INT)))
  CAST(o.target_from_coord AS INT) >= CAST(f.rel_end AS INT) - 10
  AND CAST(orf.end AS INT) <= CAST(f.rel_start AS INT) + 10
  AND (t.target_from_coord IS NULL OR (CAST(f.rel_start AS INT) >= CAST(t.target_from_coord AS INT) - 10 AND CAST(orf.end AS INT) > CAST(t.target_from_coord AS INT)))
" > ${BASE_PATH}/signatures_report_qcov${COVERAGE}_${FLANKING_SEQ_LEN}.tsv

echo "##### zipping results in ${ROOT_PATH}/${n}_qcov${COVERAGE}.zip #####"
cd ${ROOT_PATH}
rm -f ${n}_qcov${COVERAGE}_${FLANKING_SEQ_LEN}.zip
zip -r ${n}_qcov${COVERAGE}_${FLANKING_SEQ_LEN}.zip \
${n}_qcov${COVERAGE}/*_${FLANKING_SEQ_LEN}.db \
${n}_qcov${COVERAGE}/*_${FLANKING_SEQ_LEN}.tsv \
${n}_qcov${COVERAGE}/tnpA_seqs_${FLANKING_SEQ_LEN}/*.tsv \
${n}_qcov${COVERAGE}/orf121_search_${FLANKING_SEQ_LEN}/*.tsv \
${n}_qcov${COVERAGE}/oriIS_search_${FLANKING_SEQ_LEN}/*.tsv \
${n}_qcov${COVERAGE}/terIS_search_${FLANKING_SEQ_LEN}/*.tsv

cd -
echo "##### done #####"





