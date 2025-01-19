#!/bin/bash

################################################################################
# Tumor-Normal Mutation Analysis Pipeline using GATK Mutect2
#
# This script follows the GATK best practices workflow for somatic mutation
# analysis. It processes paired-end whole-genome sequencing (WGS) data for
# tumor and normal samples. Key steps include reference preparation,
# preprocessing, alignment, and mutation calling.
#
# Author: [Vinay Rajput]
# email : srm.vinay0005@gmail.com
# Repository: [https://github.com/vinayrajput0005/]
#
# Requirements:
# - GATK (v4 or later)
# - Samtools
# - picard
# - wget
# - BWA
# - Reference genome (hg38 used in this example)
# - MultiQC
#
# Note: Paths in this script are examples. Adjust to your directory structure.
################################################################################


echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################

# download reference files
wget -P /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/hg38.fa


# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/hg38.fa O=/media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/hg38.dict


# download known sites files for BQSR from GATK resource bundle
wget -P /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


################################################### Mutect2 files (TO BE DOWNLOADED ONLY ONCE) ##########################################################

echo "Download Mutect2 supporting files..."

# gnomAD
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz /media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files/
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi /media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files/

# PoN
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz /media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files/
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi /media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files/

# to create your own panel of normals: https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA

# intervals list
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list /media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files/

###################################################### VARIANT CALLING STEPS ####################################################################


# Set directories path
ref=/media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/hg38.fa
known_sites=/media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
project_dir=/media/ncim/16TB_Data/Students/Vinay/pupil
aligned_reads=$project_dir/aligned
reads=$project_dir/reads
results=$project_dir/results
mutect2_supporting_files=/media/ncim/16TB_Data/Students/Vinay/pupil/mutect2_supporting_files


# -------------------
# STEP 1: QC - Run fastqc
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz -o ${reads}/
fastqc ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz -o ${reads}/

fastqc ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz -o ${reads}/
fastqc ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz -o ${reads}/

#for i in *.gz; do fastqc $i -o fastqc/; done

# No trimming required, quality looks okay.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference
bwa index ${ref}

# BWA alignment
bwa mem -t 30 -R "@RG\tID:PA220KH-lib09-P19-Tumor_S2\tPL:ILLUMINA\tSM:PA220KH-lib09-P19-Tumor_S2" ${ref} ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2.paired.sam

bwa mem -t 30 -R "@RG\tID:PA221MH-lib09-P19-Norm_S1\tPL:ILLUMINA\tSM:PA221MH-lib09-P19-Norm_S1" ${ref} ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz > ${aligned_reads}/PA221MH-lib09-P19-Norm_S1.paired.sam


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2.paired.sam -O ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_reads.bam

gatk MarkDuplicatesSpark -I ${aligned_reads}/PA221MH-lib09-P19-Norm_S1.paired.sam -O ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_reads.bam


# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_recal_data.table

gatk BaseRecalibrator -I ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_recal_data.table



# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_recal_data.table -O ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_bqsr_reads.bam

gatk ApplyBQSR -I ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_recal_data.table -O ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_bqsr_reads.bam



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_bqsr_reads.bam  O=${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_alignment_metrics.txt

gatk/gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_insert_size_histogram.pdf


gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_bqsr_reads.bam  O=${aligned_reads}/PA221MH-lib09-P19-Norm_S1_alignment_metrics.txt

gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/PA221MH-lib09-P19-Norm_S1_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/PA221MH-lib09-P19-Norm_S1_insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Somatic Variants - Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk mutect2 caller"


gatk Mutect2 -R ${ref} \
     -I ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_bqsr_reads.bam \
     -I ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_bqsr_reads.bam \
     -tumor PA220KH-lib09-P19-Tumor_S2 \
     -normal PA221MH-lib09-P19-Norm_S1 \
     --germline-resource ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
     --panel-of-normals ${mutect2_supporting_files}/1000g_pon.hg38.vcf.gz \
     -O ${results}/PA_somatic_variants_mutect2.vcf.gz \
     --f1r2-tar-gz ${results}/PA_f1r2.tar.gz \
     --native-pair-hmm-threads 30


# ----------------------------------------------
# STEP 7: Estimate cross-sample contamination
# ----------------------------------------------

# GetPileupSummaries
# Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results are used with CalculateContamination.

# Intervals: wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list

echo "STEP 7: Estimate cross-sample contamination"

# tumor
gatk GetPileupSummaries \
    --java-options '-Xmx50G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/PA220KH-lib09-P19-Tumor_S2_sorted_dedup_bqsr_reads.bam \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/PA220KH-lib09-P19-Tumor_S2_getpileupsummaries.table

# normal
gatk GetPileupSummaries \
    --java-options '-Xmx50G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/PA221MH-lib09-P19-Norm_S1_sorted_dedup_bqsr_reads.bam  \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/PA221MH-lib09-P19-Norm_S1_getpileupsummaries.table



# Calculate contamination
gatk CalculateContamination \
    -I ${results}/PA220KH-lib09-P19-Tumor_S2_getpileupsummaries.table \
    -matched ${results}/PA221MH-lib09-P19-Norm_S1_getpileupsummaries.table \
    -O ${results}/PA_pair_calculatecontamination.table


# ----------------------------------------------
# STEP 8: Estimate read orientation artifacts
# ----------------------------------------------

echo "STEP 8: Estimate read orientation artifacts"

# read orientation
gatk LearnReadOrientationModel \
    -I ${results}/PA_f1r2.tar.gz \
    -O ${results}/read-orientation-model.tar.gz


# ----------------------------------------------
# STEP 9: Filter Variants Called By Mutect2
# ----------------------------------------------

echo "STEP 9: Filter Variants"
gatk FilterMutectCalls \
        -V ${results}/PA_somatic_variants_mutect2.vcf.gz \
        -R ${ref} \
        --contamination-table ${results}/PA_pair_calculatecontamination.table \
        --ob-priors ${results}/read-orientation-model.tar.gz \
        -O ${results}/PA_somatic_variants_filtered_mutect2.vcf


# ----------------------------------------------
# STEP 10: Annotate Variants - Funcotator
# ----------------------------------------------

echo "STEP 10: Annotate Variants"

# Annotate using Funcotator
gatk Funcotator \
    --variant ${results}/PA_somatic_variants_filtered_mutect2.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path /media/ncim/16TB_Data/Students/Vinay/pupil/reference_files/fruncotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
    --output ${results}/PA_somatic_variants_functotated.vcf \
    --output-file-format MAF
