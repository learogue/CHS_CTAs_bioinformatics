#!/bin/bash
# -----------------------------------------------------------------------------------------------------------------
# bash script : RNAseq cell lines pipeline
# Author : Léa ROGUE
# Date : 27-06-2025
# Description : This script performs a complete and reproducible RNA-seq processing pipeline for CHS cell lines. It 
# automates steps to generate gene expression matrices from raw FASTQ files.
#   1. Initial quality control of raw FASTQ files using FastQC and MultiQC
#   2. Adapter trimming and quality filtering with Cutadapt
#   3. Post-trimming quality assessment
#   4. Alignment of reads to the GRCh38 human reference genome using STAR
#   5. Gene-level quantification using featureCounts with paired-end mode
# -----------------------------------------------------------------------------------------------------------------

# Paths
RAW=data/raw
TRIM=data/trim
STAR_INDEX=data/star
GTF=data/Homo_sapiens.GRCh38.114.gtf
QC_INIT=results/qc/qc_init
QC_POST=results/qc/qc_post_trim
COUNTS=results/counts

mkdir -p $TRIM $QC_INIT $QC_POST $COUNTS

# Sample list
SAMPLES=("JJ012_S1" "CH2879_S2")

# Initial FastQC
echo "Running initial FastQC..."
fastqc $RAW/*.fastq.gz -t 6 -o $QC_INIT
multiqc $QC_INIT -o $QC_INIT/multiqc_report

# Cutadapt trimming
echo "Trimming adapters..."

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW}/${SAMPLE}_R1_001.fastq.gz"
    R2="${RAW}/${SAMPLE}_R2_001.fastq.gz"
    
    OUT1="${TRIM}/trim_${SAMPLE}_R1.fastq.gz"
    OUT2="${TRIM}/trim_${SAMPLE}_R2.fastq.gz"
    LOG="${TRIM}/cutadapt_${SAMPLE}.log"

    cutadapt \
        -a $ADAPTER_R1 -A $ADAPTER_R2 \
        -o $OUT1 -p $OUT2 \
        -q 20 -m 33 -j 6 \
        $R1 $R2 > $LOG 2>&1
done

# FastQC after trimming
echo "Running FastQC post-trim..."
fastqc $TRIM/*.gz -t 6 -o $QC_POST
multiqc $QC_POST -o $QC_POST/multiqc_report

# STAR alignment
echo "Running STAR..."

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${TRIM}/trim_${SAMPLE}_R1.fastq.gz"
    R2="${TRIM}/trim_${SAMPLE}_R2.fastq.gz"

    STAR \
        --runThreadN 6 \
        --genomeDir $STAR_INDEX \
        --readFilesIn $R1 $R2 \
        --readFilesCommand zcat \
        --outFileNamePrefix results/star/${SAMPLE}_ \
        --outSAMtype BAM SortedByCoordinate
done

# featureCounts
echo "Running featureCounts..."

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/star/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    OUT="${COUNTS}/${SAMPLE}_counts.tab"

    featureCounts -T 6 -p --countReadPairs -s 0 \
        -a $GTF -o $OUT $BAM
done

echo "Analysis completed successfully."
