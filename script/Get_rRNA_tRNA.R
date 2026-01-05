### LICENCE
# Copyright 2024-2025 Virginie Ricci

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


library(data.table)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# organism $GRC $cur_release $OUT
organism = args[1]
GRC = args[2]
cur_release = args[3]
path_fasta = args[4]



####################################################################
####################################################################
####################################################################

dataset_rRNA_tRNA = paste0(tolower(organism), '_gene_ensembl')
mart <- useEnsembl('ensembl', dataset = dataset_rRNA_tRNA, version = cur_release)
# version 111 is the currently released version from Jan 2024

####################################################################
### rRNA
rRNA = biomaRt::getBM(values="rRNA", 
               filters="biotype", 
               attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
               mart = mart)
head(rRNA); dim(rRNA)

rRNA_seq = getSequence(id=rRNA$ensembl_gene_id, type='ensembl_gene_id', seqType='gene_exon_intron', mart=mart)
head(rRNA_seq); dim(rRNA_seq)

rRNA_final = left_join(rRNA_seq, rRNA)
head(rRNA_final); dim(rRNA_final)

exportFASTA(rRNA_seq, paste0(path_fasta, organism, '.', GRC, '.', cur_release, '.rrna.ensembl.fa'))

####################################################################
### mt-rRNA
mt_rRNA = biomaRt::getBM(values="Mt_rRNA", 
                      filters="biotype", 
                      attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
                      mart = mart)
head(mt_rRNA); dim(mt_rRNA)

mt_rRNA_seq = getSequence(id=mt_rRNA$ensembl_gene_id, type='ensembl_gene_id', seqType='gene_exon_intron', mart=mart)
head(mt_rRNA_seq); dim(mt_rRNA_seq)

mt_rRNA_final = left_join(mt_rRNA_seq, mt_rRNA)
head(mt_rRNA_final); dim(mt_rRNA_final)

exportFASTA(mt_rRNA_seq, paste0(path_fasta, organism, '.', GRC, '.', cur_release, '.mt-rrna.ensembl.fa'))

####################################################################
### mt-tRNA
mt_tRNA = biomaRt::getBM(values="Mt_tRNA", 
                         filters="biotype", 
                         attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
                         mart = mart)
head(mt_tRNA); dim(mt_tRNA)

mt_tRNA_seq = getSequence(id=mt_tRNA$ensembl_gene_id, type='ensembl_gene_id', seqType='gene_exon_intron', mart=mart)
head(mt_tRNA_seq); dim(mt_tRNA_seq)

mt_tRNA_final = left_join(mt_tRNA_seq, mt_tRNA)
head(mt_tRNA_final); dim(mt_tRNA_final)

exportFASTA(mt_tRNA_seq, paste0(path_fasta, organism, '.', GRC, '.', cur_release, '.mt-trna.ensembl.fa'))






