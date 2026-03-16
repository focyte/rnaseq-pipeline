#!/usr/bin/env bash
set -euo pipefail

## ================================

## CONFIGURATION

## ================================

THREADS=12
GENOME_INDEX="../hg38Index/GCA_000001405.15_GRCh38_full_analysis_set"
SPLICE_SITES="../hg38Index/human_splice_sites.txt"
GTF_FILE="../hg38Index/GCA_000001405.15_GRCh38_primary_only.gtf"
RAW_DIR="00_raw_fastq"
PRETRIM_QC="01_fastqc_pretrim"
TRIM_DIR="02_trimmed_fastq"
POSTTRIM_QC="03_fastqc_posttrim"
UMI_DIR="04_umi_extracted_fastq"
POSTUMI_QC="05_fastqc_postumi"
SAM_DIR="06_sam"
BAM_DIR="07_bam"
SORTED_BAM_DIR="08_bam_sorted"
DEDUP_BAM_DIR="09_bam_dedup"
COUNTS_DIR="10_featurecounts"
LOG_DIR="logs"
MULTIQC_DIR="multiqc_reports"

mkdir -p $PRETRIM_QC $TRIM_DIR $POSTTRIM_QC $UMI_DIR $POSTUMI_QC \
         $SAM_DIR $BAM_DIR $SORTED_BAM_DIR $DEDUP_BAM_DIR $COUNTS_DIR \
         $LOG_DIR $MULTIQC_DIR

## ================================

## CONTROL: Set starting step

## ================================

# 1 = Pre-trim FastQC
# 2 = Trimming
# 3 = UMI extraction
# 4 = Mapping
# 5 = SAM to BAM
# 6 = Sort & Index BAM
# 7 = Deduplication
# 8 = featureCounts

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
    for fq1 in $RAW_DIR/*_R1_001.fastq.gz; do #change this based on the format of your file
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

## STEP 3: UMI extraction

## ================================

if [ "$START_STEP" -le 3 ]; then
    echo "=== Step 3: Extract UMIs ==="
    for fq1 in $TRIM_DIR/*_R1_001_val_1.fq.gz; do
        fq2=${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}
        sample=$(basename "$fq1" _R1_001_val_1.fq.gz)

        umi_tools extract \
            -I "$fq2" \
            --read2-in="$fq1" \
            --bc-pattern=NNNNNNNNCCCCCC \
            --stdout="$UMI_DIR/${sample}_R2_umi.fastq.gz" \
            --read2-out="$UMI_DIR/${sample}_R1_umi.fastq.gz"
    done
    fastqc -t $THREADS -o $POSTUMI_QC $UMI_DIR/*.fastq.gz
    multiqc -o $MULTIQC_DIR $POSTUMI_QC
fi

## ================================

## STEP 4: Mapping with HISAT2

## ================================

if [ "$START_STEP" -le 4 ]; then
    echo "=== Step 4: Mapping with HISAT2 ==="
    for fq1 in $UMI_DIR/*_R1_umi.fastq.gz; do
        fq2=${fq1/_R1_/_R2_}
        sample=$(basename "$fq1" | sed 's/_R1_umi.fastq.gz//')
        
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

## STEP 5: SAM to BAM

## ================================

if [ "$START_STEP" -le 5 ]; then
    echo "=== Step 5: SAM to BAM ==="
    for sam in $SAM_DIR/*.sam; do
        bam=$BAM_DIR/$(basename "$sam" .sam).bam
        samtools view -bS "$sam" > "$bam"
    done
fi

## ================================

## STEP 6: Sort & Index BAM

## ================================

if [ "$START_STEP" -le 6 ]; then
    echo "=== Step 6: Sort & Index BAM ==="
    for bam in $BAM_DIR/*.bam; do
        sorted=$SORTED_BAM_DIR/$(basename "$bam" .bam)_sorted.bam
        samtools sort "$bam" -o "$sorted"
        samtools index "$sorted"
    done
fi

## ================================

## STEP 7: Deduplication (per chromosome, including chrM)

## ================================

if [ "$START_STEP" -le 7 ]; then
    set -euo pipefail
    shopt -s nullglob

    echo "=== Step 7: Deduplication (per chromosome, including chrM) ==="
    for bam in $SORTED_BAM_DIR/*_sorted.bam; do
        sample=$(basename "$bam" _sorted.bam)
        tmpdir="$DEDUP_BAM_DIR/tmp_${sample}"
        mkdir -p "$tmpdir"

        # Dedup per contig (keep all references except '*')
        for chr in $(samtools idxstats "$bam" | cut -f1 | grep -v '\*'); do
            echo "  Processing $sample : $chr"
            chr_bam="$tmpdir/${sample}_${chr}.bam"
            chr_dedup="$tmpdir/${sample}_${chr}_dedup.bam"

            # Extract contig -> BAM + index
            samtools view -b -h "$bam" "$chr" -o "$chr_bam"
            samtools index "$chr_bam"

            # Deduplicate
            umi_tools dedup \
                --paired \
                --chimeric-pairs=discard \
                --unpaired-reads=discard \
                -I "$chr_bam" \
                -S "$chr_dedup"

            # Clean per-contig input
            rm -f "$chr_bam" "$chr_bam.bai"
        done

        # Merge per-chromosome dedup BAMs
        outbam="$DEDUP_BAM_DIR/${sample}_sorted_dedup.bam"
        echo "  Merging per-chromosome BAMs for $sample -> $outbam"
        samtools merge -f "$outbam" "$tmpdir/${sample}_"*_dedup.bam
        samtools index "$outbam"

        # Clean temp folder safely
        rm -rf "$tmpdir"
    done
fi

## ================================

## STEP 8: featureCounts

## ================================

if [ "$START_STEP" -le 8 ]; then
    echo "=== Step 8: featureCounts ==="
    for sorted_dedup_bam_file in $DEDUP_BAM_DIR/*_sorted_dedup.bam; do
        sample=$(basename "$sorted_dedup_bam_file" _sorted_dedup.bam)
        featureCounts -T $THREADS -p --countReadPairs \
            -t exon -a "$GTF_FILE" \
            -o "$COUNTS_DIR/${sample}_featurecounts.txt" \
            -g gene_id "$sorted_dedup_bam_file"
    done
fi

echo "=== Pipeline completed successfully! ==="
