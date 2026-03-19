#!/usr/bin/env bash
# ================================================================
#        RNA-seq QC Pipeline: Qualimap RNA-seq + BAMQC
#        Loops over all BAM files and aggregates results via MultiQC
# ================================================================

set -euo pipefail

# ================================
# USAGE
# ================================
# ./run_qualimap.sh <input_dir> <output_dir>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# ================================
# CREATE OUTPUT DIRECTORY
# ================================
mkdir -p "$OUTPUT_DIR"

# ================================
# PATH TO ANNOTATION (adjust to your system)
# ================================
GTF_FILE="/home/focyte/hg38Index/GCF_000001405.40/gencode.v49.annotation.gtf"

# ================================
# LOOP OVER BAM FILES
# ================================
for bam_file in "$INPUT_DIR"/*.bam; do
    # Skip if no BAM files exist
    [ -e "$bam_file" ] || { echo "No BAM files found in $INPUT_DIR"; exit 1; }

    sample_name=$(basename "$bam_file" .bam)
    sample_outdir="$OUTPUT_DIR/$sample_name"
    bamqc_outdir="$sample_outdir/bamqc"

    mkdir -p "$sample_outdir" "$bamqc_outdir"

    echo "----------------------------------------"
    echo "Processing sample: $sample_name"
    echo "Input BAM: $bam_file"
    echo "Output directory: $sample_outdir"
    echo "----------------------------------------"

    # ================================
    # RUN QUALIMAP RNA-seq QC
    # ================================
    qualimap rnaseq \
        -bam "$bam_file" \
        --java-mem-size=12000M \
        -gtf "$GTF_FILE" \
        -outdir "$sample_outdir" \
        -s \
        -outformat PDF \
        --paired

    # ================================
    # RUN QUALIMAP BAMQC
    # ================================
    qualimap bamqc \
        -bam "$bam_file" \
        --java-mem-size=12000M \
        -outdir "$bamqc_outdir" \
        -outformat PDF

    # ================================
    # CHECK SUCCESS
    # ================================
    if [ $? -eq 0 ]; then
        echo "✅ Finished Qualimap for $sample_name"
    else
        echo "❌ Error running Qualimap for $sample_name" >&2
    fi
done

# ================================
# AGGREGATE RESULTS WITH MULTIQC
# ================================
echo "========================================"
echo "Running MultiQC to aggregate Qualimap results..."
multiqc "$OUTPUT_DIR" -n multiqc_qualimap_report.html

echo "✅ MultiQC report generated: $OUTPUT_DIR/multiqc_qualimap_report.html"
echo "Pipeline completed successfully!"
