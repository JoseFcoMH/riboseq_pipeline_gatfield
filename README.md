## Ribo-seq Snakemake pipeline
**Pipeline overview for mouse data:**
- Download sequencing data (if needed)
    - You need to manually create a file called "samples.links" containing the links to download the Ribo-seq data.
    - The links should contain the following four arguments: LibID_LaneID_ReadID_RunID.fastq.gz. For example, "PF007_L1_R1_001.fastq.gz"
- Combine .fastq files according to LibID. For example, "PF007_L1_R1_001.fastq.gz" and "PF007_L2_R1_001.fastq.gz" will be combined as "PF007.fastq.gz".
- Perform FastQC
- Perform Trim Galore (remove adapters and do size filtering)
- Perform "UMI-tools whitelist" (identify barcodes)
- Perform "UMI-tools extract" (filter reads without a barcode), size filtering (according to monosome/disome size), quality filtering, and 2nt (from left) trimming
- Perform FastQC
- Sequential mapping using STAR
    1) mouse rRNA
    2) human rRNA
    3) mouse tRNA
    4) mouse cDNA
- Split mouse tRNA unmapped reads (.fastq) per barcode
- Prepare RSEM reference (rsem-prepare-reference) for mouse genome
- Perform mouse genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)
- Perform "UMI-tools dedup" on transcriptome-projected .bam
- Keep only forward reads of transcriptome-projected .bam
- Perform RSEM (rsem-calculate-expression)
- Create BigWig
- Extract QC statistics
- Plot QC figures

**Pipeline overview for human data:**
Same step as for mouse data, except the mapping.
- Sequential mapping using STAR
    1) human rRNA
    2) human tRNA
    3) humman cDNA
- Split human tRNA unmapped reads (.fastq) per barcode
- Prepare RSEM reference (rsem-prepare-reference) for human genome
- Perform human genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)


## RNA-seq Snakemake pipeline
