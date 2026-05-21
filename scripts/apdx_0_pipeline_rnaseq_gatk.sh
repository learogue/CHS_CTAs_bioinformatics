#!/bin/bash
# ----------------------------------------------------------------------------------------------------------------------
# Script : RNA-seq Variant Calling Pipeline (GATK)
# Author : Léa ROGUE
# Date : 18/09/2025
# Usage : ./apdx_0_pipeline_rnaseq_gatk.sh.sh SAMPLE_ID R1.fastq.gz R2.fastq.gz
#     Input : SAMPLE_ID, paired-end FASTQ files (R1/R2)
#     Output : Aligned BAM files, processed BAM for variant calling, raw VCF, and filtered VCF stored in results
# Description : This script runs a complete RNA-seq variant-calling pipeline following GATK best practices.  
# It performs :
#     1. STAR 2-pass alignment (using GRCh38 reference)
#     2. Read group addition and duplicate marking (Picard + GATK)
#     3. SplitNCigarReads processing for RNA-seq intron-aware variant calling
#     4. Variant calling using GATK HaplotypeCaller
#     5. Hard-filtering of variants for quality control
#     6. Variant annotation
# ----------------------------------------------------------------------------------------------------------------------

SAMPLE=$1
R1=$2
R2=$3
REF=../data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=../data/Homo_sapiens.GRCh38.114.gtf
STAR_INDEX=../data/star

mkdir -p ../results/${SAMPLE}

# STAR 2-pass alignment
echo "STAR 1-pass alignment"
STAR --runThreadN 16 \
     --genomeDir ${STAR_INDEX} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --outFileNamePrefix ../results/${SAMPLE}/${SAMPLE}.1pass. \
     --outSAMtype BAM Unsorted

echo "STAR 2-pass alignment"
STAR --runThreadN 16 \
     --genomeDir ${STAR_INDEX} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --sjdbFileChrStartEnd ../results/${SAMPLE}/${SAMPLE}.1pass.SJ.out.tab \
     --outFileNamePrefix ../results/${SAMPLE}/${SAMPLE}.2pass. \
     --outSAMtype BAM SortedByCoordinate

BAM_2PASS=../results/${SAMPLE}/${SAMPLE}.2pass.Aligned.sortedByCoord.out.bam

# Add read groups & mark duplicates
echo "Adding read groups and marking duplicates"
java -jar /home/ubuntu/tools/picard.jar AddOrReplaceReadGroups \
      I=${BAM_2PASS} \
      O=../results/${SAMPLE}/${SAMPLE}.RG.bam \
      RGID=group1 RGLB=library1 RGPL=Illumina RGPU=unit1 RGSM=${SAMPLE}

gatk MarkDuplicatesSpark \
      -I ../results/${SAMPLE}/${SAMPLE}.RG.bam \
      -O ../results/${SAMPLE}/${SAMPLE}.dedup.bam

samtools index ../results/${SAMPLE}/${SAMPLE}.dedup.bam

# SplitNCigarReads
echo "SplitNCigarReads"
gatk SplitNCigarReads \
      -R ${REF} \
      -I ../results/${SAMPLE}/${SAMPLE}.dedup.bam \
      -O ../results/${SAMPLE}/${SAMPLE}.split.bam

BAM_FOR_VARIANT=../results/${SAMPLE}/${SAMPLE}.split.bam

# HaplotypeCaller
echo "Variant calling with HaplotypeCaller"
gatk HaplotypeCaller \
      -R ${REF} \
      -I ${BAM_FOR_VARIANT} \
      -O ../results/${SAMPLE}/${SAMPLE}.raw.vcf \
      --dont-use-soft-clipped-bases \
      --standard-min-confidence-threshold-for-calling 20

# Hard-filtering
echo "Hard-filtering"
gatk VariantFiltration \
      -R ${REF} \
      -V ../results/${SAMPLE}/${SAMPLE}.raw.vcf \
      --filter-expression "QD < 2.0" --filter-name "LowQD" \
      --filter-expression "FS > 30.0" --filter-name "StrandBias" \
      --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
      --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
      -O ../results/${SAMPLE}/${SAMPLE}.filtered.vcf

# Variant annotation
echo "VEP annotation"
vep -i ../results/${SAMPLE}/${SAMPLE}.filtered.vcf \
    -o ../results/${SAMPLE}/${SAMPLE}.vep.vcf \
    --cache \
    --offline \
    --assembly GRCh38 \
    --vcf \
    --pick \
    --hgvs \
    --protein \
    --uniprot

echo "Finished"
