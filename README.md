In the lab of David Gatfield, we developped pipelines to process and map Ribo-seq and RNA-seq data. The Ribo-seq pipeline was based on the previous [`pipeline`](https://github.com/gatfieldlab/pipeline) [(Arpat et al. 2020)](https://pubmed.ncbi.nlm.nih.gov/32703885/). 

The pipeline includes a Snakefile, a config.yaml, software commands, and R and Python scripts. 

##  How to run the pipeline:
You can use the pipeline by cloning the repository:

`git clone https://github.com/gatfieldlab/Snakemake_pipeline.git`

and copying all files present in [`RiboSeq/mouse`](RiboSeq/mouse) in your project directory:

`cp Snakemake_pipeline/RiboSeq/mouse/* /path/to/your/project`

**Check the prequisites and edit the config.yaml**

You can create and activate your [`Snakemake`](https://anaconda.org/bioconda/snakemake) [`conda`](https://conda.io/docs/) environment:

`conda env create -f env.yaml`

`conda activate myenv`

`export PATH=$PATH:<path/to/>Snakemake_pipeline/script` # to make scripts executable from everywhere

Then, you can prepare the reference files:

`bash prepare_refs.sh`

And finally, you can run the pipeline:

`snakemake --cores n`

If you need any help, please contact virginie.ricci@unil.ch.

##  Prerequisites:
- [samtools](https://www.htslib.org/)
- [STAR](https://github.com/alexdobin/STAR)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [UMI-tools](https://github.com/CGATOxford/UMI-tools)
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
- [seqtk](https://github.com/lh3/seqtk)
- [Entrez-Direct: e-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- Modify the software versions in `config.yaml` according to your conda environment
- Edit `config.yaml` according to your dataset and reference files


# Ribo-seq Snakemake pipeline
## Pipeline overview for mouse data:
- Download sequencing data (if needed)
    - You need to manually create a file called "samples.links" containing the links to download the Ribo-seq data.
    - The links should contain the following **four** arguments: LibID_LaneID_ReadID_RunID.fastq.gz. For example, "PF007_L1_R1_001.fastq.gz"
    - The LibID should either start with **'PF'** for monosome-seq data or **'PD'** for disome-seq data.
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
- Split mouse tRNA unmapped reads (.fastq) by barcode
- Prepare RSEM reference (rsem-prepare-reference) for mouse genome
- Perform mouse genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)
- Perform "UMI-tools dedup" on transcriptome-projected .bam
- Keep only forward reads of transcriptome-projected .bam
- Perform RSEM (rsem-calculate-expression)
- Create BigWig
- Extract QC stats
- Plot QC figures

## Pipeline overview for human data:

-> Same step as for mouse data, except the mapping.
- Sequential mapping using STAR
    1) human rRNA
    2) human tRNA
    3) humman cDNA
- Split human tRNA unmapped reads (.fastq) by barcode
- Prepare RSEM reference (rsem-prepare-reference) for human genome
- Perform human genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)


# RNA-seq Snakemake pipeline
## Pipeline overview for mouse data:
- Download sequencing data (if needed)
    - You need to manually create a file called "samples.links" containing the links to download the Ribo-seq data.
    - The links should contain the following four arguments: LibID_LaneID_ReadID_RunID.fastq.gz. For example, "PF007_L1_R1_001.fastq.gz"
- Combine .fastq files according to LibID. For example, "PF007_L1_R1_001.fastq.gz" and "PF007_L2_R1_001.fastq.gz" will be combined as "PF007.fastq.gz".
- Perform FastQC
- Perform Trim Galore (remove adapters and do size filtering)
- Prepare RSEM reference (rsem-prepare-reference) for mouse genome
- Perform mouse genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)
- Perform RSEM (rsem-calculate-expression)
- Create BigWig
- Extract QC stats
- Plot QC figures
