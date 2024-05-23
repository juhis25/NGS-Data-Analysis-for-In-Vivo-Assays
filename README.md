# NGS-Data-Analysis-for-In-Vivo-Assays
Pipeline for analyzing NGS data from in-vivo assays and evaluating gene editing tool efficacy

# NGS Data Analysis for In Vivo Assays

This repository contains a bioinformatics pipeline for analyzing Next-Generation Sequencing (NGS) data from in-vivo assays. The pipeline evaluates gene editing tool efficacy, linking genomic alterations with nanoparticle delivery methods.

## Overview

The pipeline includes the following steps:
1. Quality Control
2. Trimming
3. Alignment
4. Post-Alignment Processing
5. Variant Calling
6. Variant Filtering
7. Variant Annotation
8. Reporting and Visualization

## Requirements

- **Software:**
  - FastQC
  - Trimmomatic
  - BWA
  - SAMtools
  - Picard
  - GATK
  - ANNOVAR

- **R Packages:**
  - Rsamtools
  - GenomicAlignments
  - VariantAnnotation
  - ggplot2
  - GenomicFeatures
  - ensembldb

## Setup

### Software Installation

1. **FastQC**:
   - Download from [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
   - Follow the installation instructions provided on the website.

2. **Trimmomatic**:
   - Download from [Trimmomatic GitHub](http://www.usadellab.org/cms/?page=trimmomatic).
   - Follow the installation instructions provided on the website.

3. **BWA**:
   - Install using package manager (e.g., `apt-get install bwa` for Debian-based systems).

4. **SAMtools**:
   - Install using package manager (e.g., `apt-get install samtools` for Debian-based systems).

5. **Picard**:
   - Download from [Picard GitHub](https://broadinstitute.github.io/picard/).
   - Follow the installation instructions provided on the website.

6. **GATK**:
   - Download from [GATK website](https://gatk.broadinstitute.org/).
   - Follow the installation instructions provided on the website.

7. **ANNOVAR**:
   - Download from [ANNOVAR website](https://annovar.openbioinformatics.org/en/latest/).
   - Follow the installation instructions provided on the website.

### R Package Installation

Run the following script in R to install the required packages:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("Rsamtools", "GenomicAlignments", "VariantAnnotation", "GenomicFeatures", "ensembldb"))

# Install CRAN package
install.packages("ggplot2")


