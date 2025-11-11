# Snakemake_pipeline
## Ribo-seq Snakemake pipeline

**Pipeline overview:**
- Download sequencing data (if needed)
    - You need to manually create a file called "samples.links" containing the links to download the Ribo-seq data.
    - The links should contain the following four arguments: <LIB_ID>_<LANE_ID>_<READ_ID>_<RUN_ID>.fastq.gz. For example, "PF007_L1_R1_001.fastq.gz"
- Combine fastq files according to <LIB_ID>. For example, "PF007_L1_R1_001.fastq.gz" and "PF007_L2_R1_001.fastq.gz" will be combined.
- Perform FastQC
- Perform Trim Galore
- Perform "UMI-tools whitelist"
- Perform "UMI-tools extract", size filtering, quality filtering, and 2nt (from left) trimming
- Perform FastQC
- Sequential mapping using STAR
    1) mouse rRNA
    2) human rRNA
    3) mouse tRNA
    4) mouse cDNA
- Split mouse tRNA unmapped reads per barcode
- Prepare RSEM reference (rsem-prepare-reference)
- Perform mouse genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)
- Perform "UMI-tools dedup" on transcriptome-projected .bam
- Keep only forward reads of transcriptome-projected .bam
- Perform RSEM (rsem-calculate-expression)
- Create BigWig
- Extract QC statistics
- Plot QC


## RNA-seq Snakemake pipeline
