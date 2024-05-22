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
vcf <- readVcf("output/annotated_variants.hg19_multianno.vcf", "hg19")
plot_variant_density <- function(vcf) {
  df <- as.data.frame(rowRanges(vcf))
  ggplot(df, aes(x=start, fill=seqnames)) +
    geom_histogram(binwidth=1e6) +
    theme_minimal() +
    labs(title="Variant Density", x="Genomic Position", y="Variant Count")
}
plot_variant_density(vcf)

# Save plot
ggsave("output/variant_density.png")
