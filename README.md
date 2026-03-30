In the lab of David Gatfield, we developped pipelines to process and map Ribo-seq and RNA-seq data. The Ribo-seq pipeline was based on the previous [`pipeline`](https://github.com/gatfieldlab/pipeline) [(Arpat et al. 2020)](https://pubmed.ncbi.nlm.nih.gov/32703885/) and is made for multiplexed sequencing data. 

The pipeline includes a Snakefile, a config.yaml, software commands, and R and Python scripts. 

##  How to run the pipeline:
Clone the repository:

`git clone https://github.com/gatfieldlab/Snakemake_pipeline.git`

Copy all files present in [`RiboSeq/mouse`](RiboSeq/mouse) in your project directory:

`cp Snakemake_pipeline/RiboSeq/mouse/* /path/to/your/project`

`cd /path/to/your/project`

Edit the config.yaml of the pipeline:

`vi config.yaml`

Create your [`Snakemake`](https://anaconda.org/bioconda/snakemake) [`conda`](https://conda.io/docs/) environment:

`cp env.yaml /path/to/your/project`

`conda env create -f env.yaml --prefix snakepipe`

`conda rename -p /path/to/snakepipe snakepipe` # if needed

Make scripts executable from anywhere:

`chmod 755 -R <path/to/>Snakemake_pipeline/script` # make scripts executable

`export PATH=$PATH:<path/to/>Snakemake_pipeline/script` # make scripts available

Prepare the reference files:

`vi script/prepare_refs_[mouse/human].sh`

`bash script/prepare_refs_[mouse/human].sh`

Install BBTools (https://github.com/bbushnell/BBTools) to get demuxbyname.sh:

`wget https://sourceforge.net/projects/bbmap/files/latest/download -O BBTools.tar.gz`

`tar -xzf BBTools.tar.gz`

`cd bbmap`

`./bbversion.sh`

Run the pipeline:

`snakemake --cores n`

To run the pipeline on a cluster (e.g. Curnagl cluster of UNIL), you need to open a tmux session, activate your conda environment and run your snakemake command frontend. **Do not forget to create and edit a workflow profile.**

`vi WorkFlowProfile/config.yaml`

`tmux new-session -s snakepipe_session` # `tmux attach-session -t snakepipe_session` to attach and `tmux detach` to detach

`conda activate snakepipe`

`snakemake --workflow-profile WorkFlowProfile/config.yaml --configfile config.yaml`

If you need any help, please contact virginie.ricci@unil.ch.

##  Reminder:
- All the software needed are listed in `env.yaml`
- Edit `config.yaml`
- Edit `script/prepare_refs_[mouse/human].sh`
- Edit `WorkFlowProfile/config.yaml` if you are using a cluster



# Ribo-seq Snakemake pipeline
## Pipeline overview for mouse data:
- Download sequencing data **manually** and store it in raw_data/
    - The files must be **.fastq.gz** and contain the **Library ID in the first element when split by '_'** : LibID_LaneID_ReadID_RunID.fastq.gz. For example, "PF007_L1_R1_001.fastq.gz"
    - The LibID should start with **'PF'** for monosome-seq data and **'PD'** for disome-seq data.
- Combine .fastq files according to **LibID**. For example, "PF007_L1_R1_001.fastq.gz" and "PF007_L2_R1_001.fastq.gz" will be combined as "PF007.fastq.gz".
- Perform FastQC
- Perform Trim Galore (remove adapters and apply size filtering)
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
    3) human cDNA
- Split human tRNA unmapped reads (.fastq) by barcode
- Prepare RSEM reference (rsem-prepare-reference) for human genome
- Perform human genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)


# RNA-seq Snakemake pipeline
## Pipeline overview for mouse data:
- Download sequencing data **manually** and store it in raw_data/
    - The files must be **.fastq.gz** and contain the **Library ID in the first element when split by '_'** : LibID_LaneID_ReadID_RunID.fastq.gz. For example, "RNA007_L1_R1_001.fastq.gz"
- Combine .fastq files according to **LibID**. For example, "RNA007_L1_R1_001.fastq.gz" and "RNA007_L2_R1_001.fastq.gz" will be combined as "RNA007.fastq.gz".
- Perform FastQC
- Perform Trim Galore (remove adapters and apply size filtering)
- Prepare RSEM reference (rsem-prepare-reference) for mouse genome
- Perform mouse genome mapping using STAR with RSEM parameters (mapping on genome and projection on transcriptome)
- Perform RSEM (rsem-calculate-expression)
- Create BigWig
- Extract QC stats
- Plot QC figures
