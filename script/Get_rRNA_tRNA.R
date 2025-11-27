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






