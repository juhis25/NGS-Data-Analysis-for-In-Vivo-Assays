#!/bin/bash

# Create directories for storing data
mkdir -p data/fastq
mkdir -p data/reference

# Change directory to data/fastq
cd data/fastq

# Download example data from SRA
fastq-dump --split-files SRR12345678

# Rename files for consistency (if necessary)
mv SRR12345678_1.fastq sample_R1.fastq
mv SRR12345678_2.fastq sample_R2.fastq

# Download the reference genome
cd ../reference
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh37.75.dna.primary_assembly.fa hg19.fasta
