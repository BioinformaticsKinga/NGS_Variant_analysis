# Step 1: Quality check of raw data using FASTQC
# First, we perform a quality check on the raw FASTQ files using FASTQC. This tool generates reports showing various quality metrics,
# such as the overall read quality, GC content, adapter contamination, and other potential issues with sequencing.
fastqc b1.fastq.gz b2.fastq.gz -o fastqc_reports/

# Step 2: Trimming low-quality bases using Trimmomatic
# After quality checking, we trim low-quality bases and adapter sequences from the raw reads. Trimmomatic removes unwanted sequences,
# like adapters, and trims off poor-quality bases from both ends of the reads. This step helps improve the quality of the data before alignment.
trimmomatic PE -phred33 \
  b1_R1.fastq.gz b1_R2.fastq.gz \
  b1_R1_trimmed.fastq.gz b1_R1_unpaired.fastq.gz \
  b1_R2_trimmed.fastq.gz b1_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Aligning reads to the reference genome using BWA
# In this step, we align the cleaned and trimmed reads to a reference genome. This helps us identify the exact genomic regions
# where the reads originate. BWA (Burrows-Wheeler Aligner) is a tool used to efficiently map short DNA reads to a reference genome.
bwa mem -t 4 reference_genome.fasta b1_R1_trimmed.fastq.gz b1_R2_trimmed.fastq.gz > b1_aligned.sam

# Step 4: Converting and sorting the SAM file into BAM format
# Once the reads are aligned, we convert the SAM file (which is a text file) into BAM format (a binary format), which is more efficient for processing.
# We also sort the BAM file by the coordinates of the reads, which is a necessary step before variant calling.
samtools view -bS b1_aligned.sam | samtools sort -o b1_sorted.bam
samtools index b1_sorted.bam

# Step 5: Detecting structural variants and indels with DELLY
# DELLY is a tool used to detect structural variants (SVs) and indels (insertions and deletions) from the aligned BAM files.
# Here, we focus on detecting deletions (DEL) and inversions (INV). If we wanted to find other types of structural variants,
# we would modify the variant type in the DELLY command.
delly call -t DEL -g reference_genome.fasta -o b1_deletion.vcf b1_sorted.bam
delly call -t INV -g reference_genome.fasta -o b1_inversion.vcf b1_sorted.bam

# Step 6: Annotating variants using ANNOVAR
# After detecting the variants, we use ANNOVAR to annotate them. Annotation adds functional information to the variants, such as
# which genes are affected and whether the variants have clinical significance or any known effects.
# ANNOVAR can also provide information on how the variant might affect gene function.
./table_annovar.pl b1_deletion.vcf humandb/ -buildver hg19 -out b1_annotated -remove -protocol refGene,cytoBand,exac03 -operation g,r,f -nastring .

# Step 7: Visualizing and interpreting results in R
# Finally, we visualize the annotated variants using R. We use ggplot2 to create a bar plot that shows how the variants are distributed
# across different functional categories (such as "synonymous", "nonsynonymous", etc.). This helps us better understand the impact
# of the variants on the genes and their functions.
library(ggplot2)
variants <- read.table("b1_annotated.hg19_multianno.txt", header=TRUE, sep="\t")
ggplot(variants, aes(x=Func.refGene)) +
  geom_bar() +
  theme_minimal() +
  xlab("Gene Function") +
  ylab("Number of Variants")

