# NGS Variant analysis
This pipeline processes raw sequencing data to detect and annotate variants, helping interpret their biological relevance.

Workflow summary:

Quality Control (FASTQC):
Check raw data quality (FASTQ files) for issues like contamination and low-quality bases.

Read Trimming (Trimmomatic):
Trim low-quality bases and adapter sequences to improve data quality.

Alignment (BWA):
Align cleaned reads to a reference genome to map their positions.

SAM to BAM Conversion (Samtools):
Convert SAM to BAM, sort by genomic coordinates, and index for efficient processing.

Variant Detection (DELLY):
Detect structural variants and indels from the aligned BAM file.

Variant Annotation (ANNOVAR):
Annotate variants to add functional information (affected genes, clinical significance).

Visualization (R):
Visualize variant distribution and impact using ggplot2.
