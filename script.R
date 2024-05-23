# Load necessary libraries
library(Rsamtools)
library(GenomicAlignments)
library(VariantAnnotation)
library(ggplot2)
library(GenomicFeatures)
library(ensembldb)

# Set file paths
fastq_files <- list.files(path = "data/fastq", pattern = "*.fastq.gz", full.names = TRUE)
ref_genome <- "data/reference/hg19.fasta"
output_dir <- "output/"

# Create output directories if not exist
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(file.path(output_dir, "fastqc"))) dir.create(file.path(output_dir, "fastqc"))
if (!dir.exists(file.path(output_dir, "plots"))) dir.create(file.path(output_dir, "plots"))

# Step 1: Quality Control
system("fastqc -o output/fastqc/ data/fastq/*.fastq.gz")

# Step 2: Trimming
system("trimmomatic PE -phred33 data/fastq/sample_R1.fastq.gz data/fastq/sample_R2.fastq.gz output/trimmed_R1_paired.fastq.gz output/trimmed_R1_unpaired.fastq.gz output/trimmed_R2_paired.fastq.gz output/trimmed_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36")

# Step 3: Alignment
system("bwa mem data/reference/hg19.fasta output/trimmed_R1_paired.fastq.gz output/trimmed_R2_paired.fastq.gz > output/aligned_reads.sam")

# Convert SAM to BAM, sort and index
system("samtools view -Sb output/aligned_reads.sam | samtools sort -o output/aligned_reads_sorted.bam")
system("samtools index output/aligned_reads_sorted.bam")

# Step 4: Mark Duplicates
system("picard MarkDuplicates I=output/aligned_reads_sorted.bam O=output/aligned_reads_sorted_dedup.bam M=output/markdup_metrics.txt")
system("samtools index output/aligned_reads_sorted_dedup.bam")

# Step 5: Variant Calling
system("gatk HaplotypeCaller -R data/reference/hg19.fasta -I output/aligned_reads_sorted_dedup.bam -O output/raw_variants.vcf")

# Step 6: Variant Filtering
system("gatk VariantFiltration -R data/reference/hg19.fasta -V output/raw_variants.vcf -O output/filtered_variants.vcf --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filter-name \"my_snp_filter\"")

# Step 7: Variant Annotation
system("table_annovar.pl output/filtered_variants.vcf humandb/ -buildver hg19 -out output/annotated_variants -remove -protocol refGene,cytoBand,exac03,avsnp150 -operation g,r,f,f -nastring .")

# Step 8: Visualization

# Load VCF file
vcf <- readVcf("output/annotated_variants.hg19_multianno.vcf", "hg19")

# Variant Density Plot
plot_variant_density <- function(vcf) {
  df <- as.data.frame(rowRanges(vcf))
  p <- ggplot(df, aes(x=start, fill=seqnames)) +
    geom_histogram(binwidth=1e6) +
    theme_minimal() +
    labs(title="Variant Density", x="Genomic Position", y="Variant Count")
  ggsave(filename = file.path(output_dir, "plots", "variant_density.png"), plot = p)
  return(p)
}

# Save Variant Density Plot
plot_variant_density(vcf)

# Bar Plot of Variant Types
plot_variant_types <- function(vcf) {
  info <- info(vcf)
  types <- table(info$TYPE)
  df <- as.data.frame(types)
  colnames(df) <- c("VariantType", "Count")
  p <- ggplot(df, aes(x=VariantType, y=Count, fill=VariantType)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(title="Variant Types", x="Type", y="Count")
  ggsave(filename = file.path(output_dir, "plots", "variant_types.png"), plot = p)
  return(p)
}

# Save Variant Types Plot
plot_variant_types(vcf)

# Chromosome-wise Variant Distribution
plot_chromosome_distribution <- function(vcf) {
  df <- as.data.frame(rowRanges(vcf))
  df$chromosome <- as.factor(seqnames(df))
  p <- ggplot(df, aes(x=chromosome, fill=chromosome)) +
    geom_bar() +
    theme_minimal() +
    labs(title="Chromosome-wise Variant Distribution", x="Chromosome", y="Variant Count")
  ggsave(filename = file.path(output_dir, "plots", "chromosome_distribution.png"), plot = p)
  return(p)
}

# Save Chromosome Distribution Plot
plot_chromosome_distribution(vcf)

# Variant Quality Scores Distribution
plot_quality_scores <- function(vcf) {
  df <- as.data.frame(rowRanges(vcf))
  df$QUAL <- qual(vcf)
  p <- ggplot(df, aes(x=QUAL)) +
    geom_histogram(binwidth=10, fill="blue", color="black", alpha=0.7) +
    theme_minimal() +
    labs(title="Variant Quality Scores", x="Quality Score", y="Count")
  ggsave(filename = file.path(output_dir, "plots", "quality_scores.png"), plot = p)
  return(p)
}

# Save Quality Scores Plot
plot_quality_scores(vcf)

# Read Depth Distribution
plot_depth_distribution <- function(vcf) {
  df <- as.data.frame(rowRanges(vcf))
  df$DP <- geno(vcf)$DP
  p <- ggplot(df, aes(x=DP)) +
    geom_histogram(binwidth=10, fill="green", color="black", alpha=0.7) +
    theme_minimal() +
    labs(title="Read Depth Distribution", x="Read Depth", y="Count")
  ggsave(filename = file.path(output_dir, "plots", "depth_distribution.png"), plot = p)
  return(p)
}

# Save Depth Distribution Plot
plot_depth_distribution(vcf)
