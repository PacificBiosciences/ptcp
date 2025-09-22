#!/bin/bash

set -e  # Exit on error

function print_usage {
    echo "Usage: $0 -i INPUT.bam -s SAMPLE_NAME -x SEX -r REFERENCE.fasta -t TRGT_BED -p PARAPHASE_CONFIG -g GENOME_VERSION -q QC_BED -o OUTPUT_DIR"
    echo "  -i: Input HiFi BAM file"
    echo "  -s: Sample name"
    echo "  -x: Sample sex (M or F)"
    echo "  -r: Reference genome FASTA"
    echo "  -t: TRGT BED file of tandem repeat regions"
    echo "  -p: Paraphase config YAML file"
    echo "  -g: Genome version (e.g., 38)"
    echo "  -q: QC BED file"
    echo "  -o: Output directory"
    echo "  -h: Show this help message"
    exit 1
}

function check_file_exists {
    local file_path="$1"
    local file_description="$2"
    if [ ! -f "$file_path" ]; then
        echo "Error: $file_description does not exist: $file_path"
        exit 1
    fi
}

THREADS=$(nproc)

while getopts "i:s:x:r:t:p:g:q:o:h" opt; do
    case ${opt} in
        i )
            INPUT_BAM=$OPTARG
            ;;
        s )
            SAMPLE_NAME=$OPTARG
            ;;
        x )
            SEX=$OPTARG
            ;;
        r )
            REF_FASTA=$OPTARG
            ;;
        t )
            TRGT_BED=$OPTARG
            ;;
        p )
            PARAPHASE_CONFIG=$OPTARG
            ;;
        g )
            GENOME_VERSION=$OPTARG
            ;;
        q )
            QC_BED=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        h )
            print_usage
            ;;
        \? )
            print_usage
            ;;
    esac
done

if [ -z "$INPUT_BAM" ] || [ -z "$SAMPLE_NAME" ] || [ -z "$SEX" ] || [ -z "$REF_FASTA" ] || \
   [ -z "$TRGT_BED" ] || [ -z "$PARAPHASE_CONFIG" ] || [ -z "$GENOME_VERSION" ] || \
   [ -z "$QC_BED" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required parameters"
    print_usage
fi

if [ "$SEX" != "M" ] && [ "$SEX" != "F" ]; then
    echo "Error: Sex must be either 'M' or 'F'"
    exit 1
fi

check_file_exists "$INPUT_BAM" "Input BAM file"
check_file_exists "$REF_FASTA" "Reference FASTA file"
check_file_exists "$TRGT_BED" "TRGT BED file"
check_file_exists "$PARAPHASE_CONFIG" "Paraphase config file"
check_file_exists "$QC_BED" "QC BED file"

# Setup directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/ptcp_qc"

LOG_FILE="$OUTPUT_DIR/$SAMPLE_NAME.ptcp.log"
echo "Starting PTCP direct processing at $(date)" | tee -a "$LOG_FILE"
echo "Sample: $SAMPLE_NAME, Sex: $SEX" | tee -a "$LOG_FILE"

# Prepare input
echo "Checking input BAM file..." | tee -a "$LOG_FILE"
samtools quickcheck -u -vvvvv "$INPUT_BAM" 2>&1 | tee -a "$LOG_FILE"

if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Creating reference FASTA index..." | tee -a "$LOG_FILE"
    samtools faidx "$REF_FASTA"
fi

# Step 1: Align reads to reference using pbmm2
echo "Aligning reads to reference with pbmm2..." | tee -a "$LOG_FILE"
MAPPED_BAM="$OUTPUT_DIR/${SAMPLE_NAME}.mapped.bam"
pbmm2 align \
    --log-level INFO \
    --alarms "$OUTPUT_DIR/${SAMPLE_NAME}.alarms.json" \
    --num-threads $THREADS \
    --sort \
    --preset HiFi \
    --report-json "$OUTPUT_DIR/${SAMPLE_NAME}.mapping_stats_report.json" \
    --sample "$SAMPLE_NAME" \
    -N 1 --unmapped \
    "$REF_FASTA" \
    "$INPUT_BAM" \
    "$MAPPED_BAM" 2>&1 | tee -a "$LOG_FILE"

pbindex --num-threads $THREADS "$MAPPED_BAM" 2>&1 | tee -a "$LOG_FILE"

# Step 2: TRGT processing
echo "Processing with TRGT..." | tee -a "$LOG_FILE"

# 2a: Extract reads overlapping repeats
echo "Extracting reads overlapping repeat regions..." | tee -a "$LOG_FILE"
EXPANDED_BED="$OUTPUT_DIR/${SAMPLE_NAME}.expanded.bed"
REPEATS_BAM="$OUTPUT_DIR/${SAMPLE_NAME}.repeats.bam"

awk -v sz=3000 \
    '{if ($2 < sz) { print $1"\t0\t"$3 + sz} else {print $1"\t"$2 - sz"\t"$3 + sz} }' \
    "$TRGT_BED" > "$EXPANDED_BED"

samtools view \
    --write-index \
    --use-index \
    --target-file "$EXPANDED_BED" \
    --output "${REPEATS_BAM}##idx##${REPEATS_BAM}.bai" \
    "$MAPPED_BAM" 2>&1 | tee -a "$LOG_FILE"

rm -rf "$EXPANDED_BED"

# 2b: TRGT genotyping
echo "Genotyping tandem repeats with TRGT..." | tee -a "$LOG_FILE"
TRGT_OUT_PREFIX="$OUTPUT_DIR/${SAMPLE_NAME}.trgt"
KARYOTYPE=$([ "$SEX" == "M" ] && echo "XY" || echo "XX")

trgt --verbose genotype \
    --threads $THREADS \
    --preset targeted \
    --genome "$REF_FASTA" \
    --karyotype "$KARYOTYPE" \
    --reads "$MAPPED_BAM" \
    --repeats "$TRGT_BED" \
    --output-prefix "$TRGT_OUT_PREFIX" 2>&1 | tee -a "$LOG_FILE"

gunzip -c "${TRGT_OUT_PREFIX}.vcf.gz" > "${TRGT_OUT_PREFIX}.vcf"
rm -rf "${TRGT_OUT_PREFIX}.vcf.gz"
samtools sort -o "${TRGT_OUT_PREFIX}.sorted.spanning.bam" "${TRGT_OUT_PREFIX}.spanning.bam"
rm -rf "${TRGT_OUT_PREFIX}.spanning.bam"
samtools index "${TRGT_OUT_PREFIX}.sorted.spanning.bam"

# 2c: Create plots for TRGT results
echo "Creating visualization plots for TRGT results..." | tee -a "$LOG_FILE"

# Function to create TRGT plots
create_trgt_plots() {
    local show=$1
    local plot_type=$2
    local output_prefix="$OUTPUT_DIR/${SAMPLE_NAME}.${show}_${plot_type}"

    grep -v "^#" "${TRGT_OUT_PREFIX}.vcf" | awk -F"[=:\t;]" '{print $9}' | while read -r i; do
        echo "Creating $show $plot_type plot for repeat $i" | tee -a "$LOG_FILE"
        trgt plot \
            --repeat-id "$i" \
            --genome "$REF_FASTA" \
            --repeats "$TRGT_BED" \
            --spanning-reads "${TRGT_OUT_PREFIX}.sorted.spanning.bam" \
            --vcf "${TRGT_OUT_PREFIX}.vcf" \
            --show "$show" \
            --plot-type "$plot_type" \
            --image "${output_prefix}.$(printf '%s' "$i").trgt_plot.svg" 2>&1 | tee -a "$LOG_FILE"
    done

    # Use zip command instead of python zipfile module to avoid user ID issues
    python3 -m zipfile -c "${output_prefix}.trgt_plots.zip" "${output_prefix}"*.trgt_plot.svg
    rm -f "${output_prefix}"*.trgt_plot.svg
}

# Create motifs and methylation plots in both allele and waterfall views
create_trgt_plots "motifs" "allele"
create_trgt_plots "meth" "allele"
create_trgt_plots "motifs" "waterfall"
create_trgt_plots "meth" "waterfall"

# Step 3: Paraphase processing
echo "Processing with Paraphase..." | tee -a "$LOG_FILE"

# 3a: Run paraphase
PARAPHASE_OUT_DIR="$OUTPUT_DIR/${SAMPLE_NAME}_paraphase"
echo "Running Paraphase to phase haplotypes..." | tee -a "$LOG_FILE"

paraphase \
    --bam "$MAPPED_BAM" \
    --reference "$REF_FASTA" \
    --out "$PARAPHASE_OUT_DIR" \
    --genome "$GENOME_VERSION" \
    --threads $THREADS \
    --config "$PARAPHASE_CONFIG" \
    --write-nocalls-in-vcf \
    --targeted 2>&1 | tee -a "$LOG_FILE"

# 3b: F8 inversion calling
echo "Calling F8 inversions..." | tee -a "$LOG_FILE"
f8_inversion.py \
    --bam "${PARAPHASE_OUT_DIR}/${SAMPLE_NAME}.paraphase.bam" \
    --prefix "$SAMPLE_NAME" \
    --json \
    --sex "$SEX" \
    --out "$OUTPUT_DIR/" 2>&1 | tee -a "$LOG_FILE"

# Step 4: Create manifest file for ptcp-qc
echo "Creating manifest for ptcp-qc analysis..." | tee -a "$LOG_FILE"

cat <<EOF > "$OUTPUT_DIR/ptcp_qc/input_manifest.json"
{
  "sample_names": ["$SAMPLE_NAME"],
  "mapped_bams": ["$MAPPED_BAM"],
  "trgt_vcfs": ["${TRGT_OUT_PREFIX}.vcf"],
  "trgt_spanning_bams": ["${TRGT_OUT_PREFIX}.sorted.spanning.bam"],
  "paraphase_bams": ["${PARAPHASE_OUT_DIR}/${SAMPLE_NAME}.paraphase.bam"],
  "paraphase_jsons": ["${PARAPHASE_OUT_DIR}/${SAMPLE_NAME}.paraphase.json"],
  "f8_jsons": ["$OUTPUT_DIR/${SAMPLE_NAME}.f8inversion.json"]
}
EOF

# Step 5: Run ptcp-qc
echo "Running ptcp-qc analysis..." | tee -a "$LOG_FILE"

ptcp-qc -v analyze \
    --manifest "$OUTPUT_DIR/ptcp_qc/input_manifest.json" \
    --targets-bed "$QC_BED" \
    --output-prefix "$OUTPUT_DIR/ptcp_qc/qc" \
    --threads $THREADS \
    --json 2>&1 | tee -a "$LOG_FILE"

mv "$OUTPUT_DIR/ptcp_qc/qc.${SAMPLE_NAME}.json" "$OUTPUT_DIR/${SAMPLE_NAME}.qc.json"
rm -rf "$OUTPUT_DIR/ptcp_qc"

echo "PTCP direct processing completed at $(date)" | tee -a "$LOG_FILE"
echo "All results are in $OUTPUT_DIR"
echo "Done!"
