# GATK Mutect2 Pipeline for Somatic Mutation Calling in Tumor-Normal Samples

This repository contains a comprehensive pipeline designed for the analysis of somatic mutations in tumor-normal pair sequencing data using GATK (Genome Analysis Toolkit). The pipeline utilizes Mutect2 for somatic variant calling and GATK tools for subsequent filtering and annotation of detected mutations. It is suitable for whole genome sequencing (WGS) data analysis and can be easily adapted for other sequencing types.

## Components:

 - Preprocessing: Quality control, trimming, and alignment of raw sequencing reads.
 - Variant Calling: Uses Mutect2 for tumor-normal pair somatic variant detection.
 - Variant Filtering: Implements best practice filters for somatic variants to ensure accurate results.
 - Annotation: Adds biological annotations, including gene names, mutation types, and pathogenicity predictions.
 - Visualization: Generates basic visual outputs, such as variant plots or summary reports, for easy interpretation of results.

## Prerequisites

The following software tools and libraries are required to run the pipeline:

- GATK (v4 or later)
- Samtools
- Picard
- BWA
- FastQC
- MultiQC
- Wget
- Reference genome (hg38 used in this example)

## Pipeline Overview

The pipeline consists of the following key steps:

1. **Prepare reference files**:
    - Download reference genome and associated files.
    - Create required index and dictionary files for the reference.

2. **Download supporting files for Mutect2**:
    - Download the GnomAD database, Panel of Normals (PoN), and target interval list.

3. **Quality Control (QC)**:
    - Run FastQC on raw FASTQ files to check the quality of sequencing data.

4. **Mapping**:
    - Map reads to the reference genome using BWA-MEM.

5. **Mark Duplicates and Sort**:
    - Use GATK's MarkDuplicatesSpark to mark duplicate reads and sort them.

6. **Base Quality Recalibration**:
    - Perform base quality score recalibration (BQSR) using known variant sites and recalibration files.

7. **Collect Alignment Metrics**:
    - Generate alignment and insert size metrics using GATK's CollectAlignmentSummaryMetrics and CollectInsertSizeMetrics.

8. **Variant Calling**:
    - Call somatic mutations using GATK's Mutect2, providing tumor and normal samples, known germline resources, and Panel of Normals.

9. **Cross-Sample Contamination Estimation**:
    - Estimate contamination in the tumor sample using GetPileupSummaries and CalculateContamination.

## Usage

### Clone the Repository

To use the pipeline, clone the repository to your local machine:

git clone https://github.com/yourusername/Tumor-Normal-Mutation-Analysis.git
cd Tumor-Normal-Mutation-Analysis

### Setting Up Reference Files

Ensure that all reference files are downloaded and placed in the appropriate directories. Adjust the paths in the script to match your local directory structure.
## Example reference download step:

  wget -P /path/to/reference/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

This command downloads the hg38 reference genome from UCSC and places it in the specified directory (/path/to/reference/). You can replace the path with the location where you'd like to store the reference files on your system.

### Running the Pipeline

1. Modify the script to set paths to your reference files, reads, and result directories.
2. Make the script executable:

     chmod a+x tumor_normal_mutation_analysis.sh

3. Run the pipeline:

    ./tumor_normal_mutation_analysis.sh

### Input Files

> FASTQ files for tumor and normal samples (*_R1_001.fastq.gz, *_R2_001.fastq.gz).
> Reference genome (hg38) and associated files (e.g., hg38.fa, hg38.dict).
> Known sites files (e.g., dbsnp138.vcf) for base quality recalibration.

### Output Files

The pipeline generates various output files:
  > Alignment and Duplicate Metrics: Metrics for alignment quality and insert size.
  > Somatic Variant Calls: VCF file containing somatic mutations detected by Mutect2 (PA_somatic_variants_mutect2.vcf.gz).
  > Contamination Estimates: Tables for cross-sample contamination estimates (getpileupsummaries.table).

## Optional: Custom Panel of Normals (PoN)

You can create your own Panel of Normals using GATK's CreateSomaticPanelOfNormals. Refer to the official documentation for instructions.

### Acknowledgments

  > GATK: Genome Analysis Toolkit used for variant calling and processing.
  > Broad Institute: Development of reference files and resources for somatic mutation analysis.

## Feel free to adjust the content, especially regarding the installation steps and usage, depending on your specific project structure.
