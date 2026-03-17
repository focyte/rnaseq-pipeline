#!/usr/bin/env bash
set -euo pipefail

## ================================

## CONFIGURATION

## ================================

THREADS=12

GENOME_INDEX="../hg38Index/GCF_000001405.40/GCA_000001405.15_GRCh38_no_alt_analysis_set" ## Replace with your index set
SPLICE_SITES="../hg38Index/GCF_000001405.40/human_splice_sites.txt" ## Replace with your splice site file
GTF_FILE="../hg38Index/GCF_000001405.40/gencode.v49.annotation.gtf" ## Replace with your gtf file

RAW_DIR="00_raw_fastq" # Add your fastq files here
PRETRIM_QC="01_fastqc_pretrim"
TRIM_DIR="02_trimmed_fastq"
POSTTRIM_QC="03_fastqc_posttrim"
SAM_DIR="06_sam"
BAM_DIR="07_bam"
SORTED_BAM_DIR="08_bam_sorted"
COUNTS_DIR="10_featurecounts"
LOG_DIR="logs"
MULTIQC_DIR="multiqc_reports"

mkdir -p $PRETRIM_QC $TRIM_DIR $POSTTRIM_QC \
         $SAM_DIR $BAM_DIR $SORTED_BAM_DIR $COUNTS_DIR \
         $LOG_DIR $MULTIQC_DIR


## CONTROL: Set starting step

# 1 = Pre-trim FastQC
# 2 = Trimming
# 3 = Mapping
# 4 = SAM to BAM
# 5 = Sort & Index BAM
# 6 = featureCounts

START_STEP=1  # Change this to skip earlier steps

## ================================

## STEP 1: FastQC BEFORE trimming

## ================================

if [ "$START_STEP" -le 1 ]; then
    echo "=== Step 1: FastQC BEFORE trimming ==="
    fastqc -t $THREADS -o $PRETRIM_QC $RAW_DIR/*.fastq.gz
    multiqc -o $MULTIQC_DIR $PRETRIM_QC
fi

## ================================

## STEP 2: Fastp trimming

## ================================

if [ "$START_STEP" -le 2 ]; then
    echo "=== Step 2: Trimming with fastp ==="
    for fq1 in $RAW_DIR/*_R1_001.fastq.gz; do # Change this based on the format of your file
        fq2=${fq1/_R1_001.fastq.gz/_R2_001.fastq.gz}

        base1=$(basename "$fq1" .fastq.gz)
        base2=$(basename "$fq2" .fastq.gz)
        sample_prefix=${base1%_R1_001}

        fastp -i "$fq1" -I "$fq2" \
              -o "$TRIM_DIR/${base1}_val_1.fq.gz" \
              -O "$TRIM_DIR/${base2}_val_2.fq.gz" \
              -w $THREADS \
              --detect_adapter_for_pe \
              --length_required 30 \
              --html "$LOG_DIR/${sample_prefix}_fastp.html" \
              --json "$LOG_DIR/${sample_prefix}_fastp.json"
    done
    fastqc -t $THREADS -o $POSTTRIM_QC $TRIM_DIR/*.fq.gz
    multiqc -o $MULTIQC_DIR $POSTTRIM_QC
fi

## ================================

## STEP 3: Mapping with HISAT2

## ================================

if [ "$START_STEP" -le 3 ]; then
    echo "=== Step 3: Mapping with HISAT2 ==="
    for fq1 in $TRIM_DIR/*_val_1.fq.gz; do
        fq2=${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}
        sample=$(basename "$fq1" | sed 's/_R1_001_val_1.fq.gz//')
        
        hisat2 -x $GENOME_INDEX \
            --known-splicesite-infile $SPLICE_SITES \
            -p $THREADS \
            -1 "$fq1" \
            -2 "$fq2" \
            -S "$SAM_DIR/${sample}.sam" \
            2> "$LOG_DIR/${sample}_hisat2.log"
    done
fi

## ================================

## STEP 4: SAM to BAM

## ================================

if [ "$START_STEP" -le 4 ]; then
    echo "=== Step 4: SAM to BAM ==="
    for sam in $SAM_DIR/*.sam; do
        bam=$BAM_DIR/$(basename "$sam" .sam).bam
        samtools view -bS "$sam" > "$bam"
    done
fi

## ================================

## STEP 5: Sort & Index BAM

## ================================

if [ "$START_STEP" -le 5 ]; then
    echo "=== Step 5: Sort & Index BAM ==="
    for bam in $BAM_DIR/*.bam; do
        sorted=$SORTED_BAM_DIR/$(basename "$bam" .bam)_sorted.bam
        samtools sort "$bam" -o "$sorted"
        samtools index "$sorted"
    done
fi

## ================================

## STEP 6: featureCounts

## ================================

if [ "$START_STEP" -le 6 ]; then
    echo "=== Step 6: featureCounts ==="
    for sorted_bam_file in $SORTED_BAM_DIR/*_sorted.bam; do
        sample=$(basename "$sorted_bam_file" _sorted.bam)
        featureCounts -T $THREADS -p --countReadPairs \
            -t exon -a "$GTF_FILE" \
            -o "$COUNTS_DIR/${sample}_featurecounts.txt" \
            -g gene_id "$sorted_bam_file"
    done
fi

echo "=== Pipeline completed successfully! ==="
