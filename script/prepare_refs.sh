#!/bin/bash

### LICENCE
# Copyright (c) 2024-2025 Virginie Ricci

# This file is part of Snakemake_pipeline (https://github.com/gatfieldlab/Snakemake_pipeline/).

# Snakemake_pipeline is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.

# Snakemake_pipeline is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Snakemake_pipeline. 
# If not, see <https://www.gnu.org/licenses/>.
### LICENCE 

# module load samtools/1.9
# module load STAR/2.7.11b
# module load python3/3.7.3
# module load R/4.4.0



### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
export PATH=$PATH:/home/vricci/miniconda3/bin/
# efetch, esearch, ...

export PATH=$PATH:<path/to/>Snakemake_pipeline/script
# convert_ensembl_ids.py and Get_rRNA_tRNA.R


GRC=GRCm39
organism=Mmusculus
folder_ncbi=mus_musculus
file_ncbi=Mus_musculus
cur_release=111
path_ref=/data/databases/mouse/


#tRNA_link=https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.fa # previous link, same file
tRNA_link=https://gtrnadb.org/genomes/eukaryota/Mmusc39/mm39-tRNAs.fa 

preribo_45S='NR_046233.2'
### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 


# URL positions
BASE_URL=ftp://ftp.ensembl.org/pub/release-

CDNA_DIR=/fasta/${folder_ncbi}/cdna/
CDNA_FILE=${file_ncbi}.${GRC}.cdna.all.fa.gz
CDNA_OUTPUT=${organism}.${GRC}.${cur_release}.cdna.ensembl.fa

DNA_DIR=/fasta/${folder_ncbi}/dna/
DNA_FILE=${file_ncbi}.${GRC}.dna.primary_assembly.fa.gz
DNA_OUTPUT=${organism}.${GRC}.${cur_release}.dna.ensembl.fa

GTF_DIR=/gtf/${folder_ncbi}/
GTF_FILE=${file_ncbi}.${GRC}.${cur_release}.chr.gtf.gz
GTF_OUTPUT=${organism}.${GRC}.${cur_release}.gtf

GFF_DIR=/gff3/${folder_ncbi}/
GFF_FILE=${file_ncbi}.${GRC}.${cur_release}.chr.gff3.gz
GFF_OUTPUT=${organism}.${GRC}.${cur_release}.gff3

STAR_DIR=${path_ref}star/${organism}.${GRC}.${cur_release}

cur_log_file=${PWD}/fetch_ensembl_rel-${cur_release}.log

# checks
samtools_ver=$(samtools --version | head -1 | cut -d' ' -f2)
star_ver=$(STAR --version)

# logs
echo $cur_log_file
echo Fetching ENSEMBL release ${cur_release} > ${cur_log_file}
echo samtools: ${samtools_ver} >> ${cur_log_file} 
echo STAR: ${star_ver} >> ${cur_log_file}


# cDNA
cur_cdna_url=${BASE_URL}${cur_release}${CDNA_DIR}${CDNA_FILE}
cur_cdna_output=${path_ref}fasta/${CDNA_OUTPUT/${cur_release}/${cur_release}}
echo Downloading cDNA sequences from ${cur_cdna_url}
mkdir -p ${path_ref}fasta
if [ ! -f ${cur_cdna_output} ]; then
  wget --progress=dot:giga --append-output=${cur_log_file} -O - ${cur_cdna_url} \
    | gzip -dc | convert_ensembl_ids.py > ${cur_cdna_output}
    samtools faidx ${cur_cdna_output}
fi


# DNA
cur_dna_url=${BASE_URL}${cur_release}${DNA_DIR}${DNA_FILE}
cur_dna_output=${path_ref}fasta/${DNA_OUTPUT/${cur_release}/${cur_release}}
echo Downloading DNA sequences from ${cur_dna_url}
if [ ! -f ${cur_dna_output} ]; then
  wget --progress=dot:giga --append-output=${cur_log_file} -O - ${cur_dna_url} \
    | gzip -dc > ${cur_dna_output}
    samtools faidx ${cur_dna_output}
    cat ${cur_dna_output}.fai | awk '{print $1,$2}' > ${cur_dna_output}.fai.tsv
fi


# GTF
cur_gtf_url=${BASE_URL}${cur_release}${GTF_DIR}${GTF_FILE/${cur_release}/${cur_release}}
cur_gtf_output=${path_ref}gtf/${GTF_OUTPUT/${cur_release}/${cur_release}}
export cur_gtf_sorted=${cur_gtf_output/\.gtf/\.sorted\.gtf}
echo Downloading GTF from ${cur_gtf_url}
mkdir -p ${path_ref}gtf
wget --progress=dot:giga --append-output=${cur_log_file} -O - ${cur_gtf_url} \
  | gzip -dc > ${cur_gtf_output}
echo Sorting GTF
head -n 100 ${cur_gtf_output} | grep -e '^#' > ${cur_gtf_sorted}
header_len=$((1 + $(wc -l < ${cur_gtf_sorted})))
echo   .. splitting into chr files
tail -n +${header_len} ${cur_gtf_output} \
  | awk '{print $0 > "temp_sort_split_"$1}'
temp_files=$(ls temp_sort_split_* | sort)
ordered_temp_files=$( echo ${temp_files//_/ } \
  | join -11 -24 ${path_ref}chrom_order.txt - \
  | sort -k2,2n | awk '{ print "temp_sort_split_"$1}' )
echo   .. sorting split files
echo ${ordered_temp_files} \
  | xargs -P 6 -I '{}' bash -c 'sort -k4,4n -o {} {}'
echo   .. joining and cleaning split files
echo ${ordered_temp_files} \
  | xargs -P 1 -I '{}' bash -c 'rm {}'



# GFF
cur_gff_url=${BASE_URL}${cur_release}${GFF_DIR}${GFF_FILE/${cur_release}/${cur_release}}
cur_gff_output=${path_ref}gtf/${GFF_OUTPUT/${cur_release}/${cur_release}}
export cur_gff_sorted=${cur_gff_output/\.gff3/\.sorted\.gff3}
echo Downloading GFF from ${cur_gff_url}
mkdir -p ${path_ref}gtf
if [ ! -f ${cur_gff_output} ]; then
  wget --progress=dot:giga --append-output=${cur_log_file} -O - ${cur_gff_url} \
    | gzip -dc > ${cur_gff_output}
fi




### rRNA and tRNA
echo Get rRNA and tRNA
Rscript Get_rRNA_tRNA.R $organism $GRC $cur_release ${path_ref}fasta/

# rRNA + mt-rRNA
cat ${path_ref}fasta/${organism}.${GRC}.${cur_release}.rrna.ensembl.fa ${path_ref}fasta/${organism}.${GRC}.${cur_release}.mt-rrna.ensembl.fa > ${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-rrna.ensembl.fa


# tRNA + mt-tRNA
echo "Are you sure  $tRNA_link  is the right file? This should be manually checked!"
cur_trna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.trna.ensembl.fa
if [ ! -f ${cur_trna_output} ]; then
  wget --progress=dot:giga --append-output=${cur_log_file} --no-check-certificate -O - $tRNA_link > ${cur_trna_output}
fi

cat ${path_ref}fasta/${organism}.${GRC}.${cur_release}.trna.ensembl.fa ${path_ref}fasta/${organism}.${GRC}.${cur_release}.mt-trna.ensembl.fa > ${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-trna.ensembl.fa


# 45S pre-ribosomal rRNA included in mouse rRNA
echo Are you sure the mouse 45S pre-ribosomal rRNA is $preribo_45S ? This should be manually checked!
# https://www.ncbi.nlm.nih.gov/nuccore/NR_046233.2
PreRibo_45S_fa=${path_ref}fasta/${organism}.${GRC}.${cur_release}.PreRibo_45S.fa
touch $PreRibo_45S_fa
efetch -db nuccore -id $preribo_45S -format fasta > $PreRibo_45S_fa

cat ${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-rrna.ensembl.fa $PreRibo_45S_fa > ${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-rrna.PreRibo45s.ensembl.fa



# STAR reference building
echo Building ${organism} STAR all-rRNA PreRibo database for release-${cur_release}
cur_star_dir=${STAR_DIR/${cur_release}/${cur_release}}/all_rRNA/
mkdir -p ${cur_star_dir}
cur_rrna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-rrna.ensembl.fa
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $cur_star_dir \
     --genomeFastaFiles $cur_rrna_output \
     --genomeSAindexNbases 6


echo Building ${organism} STAR all-rRNA+45S PreRibo database for release-${cur_release}
cur_star_dir=${STAR_DIR/${cur_release}/${cur_release}}/all_rRNA_PreRibo/
mkdir -p ${cur_star_dir}
cur_rrna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-rrna.PreRibo45s.ensembl.fa
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $cur_star_dir \
     --genomeFastaFiles $cur_rrna_output \
     --genomeSAindexNbases 6


echo Building ${organism} STAR all-tRNA database for release-${cur_release}
cur_star_dir=${STAR_DIR/${cur_release}/${cur_release}}/all_tRNA/
mkdir -p ${cur_star_dir}
cur_trna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.all-trna.ensembl.fa
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $cur_star_dir \
     --genomeFastaFiles $cur_trna_output \
     --genomeSAindexNbases 6


echo Building ${organism} STAR cDNA database for release-${cur_release}
cur_star_dir=${STAR_DIR/${cur_release}/${cur_release}}/cDNA/
mkdir -p ${cur_star_dir}
cur_cdna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.cdna.ensembl.fa
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $cur_star_dir \
     --genomeFastaFiles $cur_cdna_output \
     --limitGenomeGenerateRAM 81027662090 \
     --genomeSAindexNbases 12


echo Building ${organism} STAR genome database for release-${cur_release}
cur_star_dir=${STAR_DIR/${cur_release}/${cur_release}}
mkdir -p ${cur_star_dir}
cur_dna_output=${path_ref}fasta/${organism}.${GRC}.${cur_release}.dna.ensembl.fa
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $cur_star_dir \
     --genomeFastaFiles $cur_dna_output --sjdbGTFfile $cur_gtf_output \
     --sjdbOverhang 100
 cat ${STAR_DIR}/Log.out >> ${cur_log_file} && rm ${STAR_DIR}/Log.out



# Get seqkit stats of fasta files
seqkit stats ${cur_rrna_output} > ${path_ref}fasta/${organism}.${GRC}.${cur_release}.seqkit_stats.txt
seqkit stats ${cur_trna_output} >> ${path_ref}fasta/${organism}.${GRC}.${cur_release}.seqkit_stats.txt
seqkit stats ${cur_cdna_output} >> ${path_ref}fasta/${organism}.${GRC}.${cur_release}.seqkit_stats.txt
seqkit stats ${cur_dna_output} >> ${path_ref}fasta/${organism}.${GRC}.${cur_release}.seqkit_stats.txt
cut -f1,2 ${cur_dna_output}.fai > ${path_ref}fasta/${organism}.${GRC}.${cur_release}.chromsizes.txt 

