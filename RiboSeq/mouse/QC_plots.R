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

library(ggplot2)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(ggpubr)

args<-commandArgs(TRUE)

work_dir=args[1]
# work_dir='/data/vricci/projects/Julien/'
raw_data=paste0(work_dir, 'raw_data/')
stats_data=paste0(work_dir, 'stats/')
path_plot=paste0(stats_data, 'plots/')
system(paste0('mkdir -p ', path_plot))

samples_df = read.csv(list.files(path=raw_data, pattern='_data.txt', full.names=T), header=FALSE, sep='\t')
names(samples_df) = 'Sample'
head(samples_df); dim(samples_df)
samples_df_uniq = unique(do.call(rbind, strsplit(samples_df$Sample, '_'))[,1])
samples_df_uniq = samples_df_uniq[c(grep('PF', samples_df_uniq), grep('PD', samples_df_uniq))]
samples_df_uniq

stats_samples = data.frame(matrix(nrow=0, ncol=4))
# for (s in samples_df$Sample){
#     stats_s = read.csv(paste0(stats_data, s, '_R1_stats.txt'), header=FALSE, sep='\t')
for (s in samples_df_uniq){
    stats_s = read.csv(paste0(stats_data, s, '_stats.txt'), header=FALSE, sep='\t')
    names(stats_s) = c('Sample', 'Step', 'Type', 'Counts')
    head(stats_s); dim(stats_s)
    stats_samples = rbind(stats_samples, stats_s)
    
}
stats_samples$Sample = str_replace(stats_samples$Sample, '_R1', '')
Sample = unique(sort(stats_samples$Sample))
stats_samples$Sample = factor(stats_samples$Sample, levels=c(Sample[grep('PF', Sample)], Sample[grep('PD', Sample)]))
Sample = unique(sort(stats_samples$Sample))
stats_samples$Sample = factor(stats_samples$Sample, levels=c(Sample[grep('PF', Sample)], Sample[grep('PD', Sample)]))
head(stats_samples); dim(stats_samples)
stats_samples$labels = paste0(stats_samples$Step, '_', stats_samples$Type)

stats_samples = stats_samples[!is.na(stats_samples$Counts),]
stats_samples$Counts = as.numeric(stats_samples$Counts)
stats_samples_BC = stats_samples[stats_samples$Type == 'splitseq_to_map' & stats_samples$Counts < 1000,]$Step
stats_samples = stats_samples[! stats_samples$Step %in% stats_samples_BC,]
head(stats_samples); dim(stats_samples)


### raw_fastq
raw_df = stats_samples[stats_samples$Step == 'raw_fastq',]
head(raw_df); dim(raw_df)


### cutadapt
cutadapt_df = stats_samples[stats_samples$Step == 'cutadapt',]
cutadapt_df = cutadapt_df[cutadapt_df$Type== 'input',]
head(cutadapt_df); dim(cutadapt_df)


### umi_tools
umi_tools_df = stats_samples[stats_samples$Step == 'umi_tools',]
umi_tools_df = umi_tools_df[umi_tools_df$Type %in% c('input', 'output'),]
head(umi_tools_df); dim(umi_tools_df)


### consume
consume_df = stats_samples[stats_samples$Step == 'consume',]
head(consume_df); dim(consume_df)


### fastq_quality_filter
fastq_quality_filter_df = stats_samples[stats_samples$Step == 'fastq_quality_filter',]
fastq_quality_filter_df = fastq_quality_filter_df[fastq_quality_filter_df$Type != 'low_qual',]
head(fastq_quality_filter_df); dim(fastq_quality_filter_df)

for_verif = fastq_quality_filter_df[fastq_quality_filter_df$Type == 'seq_to_map',]
for_verif_tbl = for_verif$Counts
names(for_verif_tbl) = for_verif$Sample
for_verif_tbl = sort(for_verif_tbl)
for_verif_tbl


### STAR mouse_rRNA, human_rRNA, mouse_tRNA
STAR_df = stats_samples[stats_samples$Step %in% c('STAR_mouse_rRNA', 'STAR_human_rRNA', 'STAR_mouse_tRNA'),]
STAR_df$Step = str_replace(STAR_df$Step, 'STAR_', '')
head(STAR_df); dim(STAR_df)

STAR_df_simple = STAR_df[STAR_df$Type %in% c('uniq_map', 'multi_map', 'unmapped'),]
head(STAR_df_simple); dim(STAR_df_simple)

STAR_df_detailed = STAR_df[STAR_df$Type %in% c('uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric'),]
head(STAR_df_detailed); dim(STAR_df_detailed)


### STAR mouse cDNA
STAR_cDNA_df = stats_samples[stats_samples$Step == 'STAR_mouse_cDNA',]
STAR_cDNA_df$Step = str_replace(STAR_cDNA_df$Step, 'STAR_', '')
head(STAR_cDNA_df); dim(STAR_cDNA_df)

STAR_cDNA_df_simple = STAR_cDNA_df[STAR_cDNA_df$Type %in% c('uniq_map', 'multi_map', 'unmapped'),]
head(STAR_cDNA_df_simple); dim(STAR_cDNA_df_simple)

STAR_cDNA_df_detailed = STAR_cDNA_df[STAR_cDNA_df$Type %in% c('uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric'),]
head(STAR_cDNA_df_detailed); dim(STAR_cDNA_df_detailed)


### STAR unmapped mouse tRNA >> what will map to genome
STAR_unmapped_mouse_tRNA_df = stats_samples[stats_samples$Type == 'splitseq_to_map',]
head(STAR_unmapped_mouse_tRNA_df); dim(STAR_unmapped_mouse_tRNA_df)


### STAR genome
STAR_genome_df = stats_samples[grep('^[ATCG]', stats_samples$Step),]
STAR_genome_df = STAR_genome_df[STAR_genome_df$Type != 'splitseq_to_map', ]
head(STAR_genome_df); dim(STAR_genome_df)

STAR_genome_df_simple = STAR_genome_df[STAR_genome_df$Type %in% c('uniq_map', 'multi_map', 'unmapped'),]
head(STAR_genome_df_simple); dim(STAR_genome_df_simple)

STAR_genome_df_detailed = STAR_genome_df[STAR_genome_df$Type %in% c('uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric'),]
head(STAR_genome_df_detailed); dim(STAR_genome_df_detailed)


### 

my_theme = theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
    strip.text.y.right = element_text(angle = 0,  size=5), 
    strip.background = element_rect(colour="black", fill="white"), legend.title=element_blank())


################################################
# raw_fastq input == cutadapt input
# cutadapt output == umi_tools input
# umi_tools output - consume filtered_out == fastq_quality_filter size_ok

colsBar = brewer.pal(n = 8, name = "Dark2")[c(1:5)]
names(colsBar) = c('input', 'cutadapt', 'umi_tools', 'consume', 'fastq_quality_filter')

df_Preprocessing = do.call("rbind", list(cutadapt_df, umi_tools_df, fastq_quality_filter_df))
head(df_Preprocessing); dim(df_Preprocessing)

df_Preprocessing[df_Preprocessing$labels == 'cutadapt_input', ]$labels = 'input'
df_Preprocessing[df_Preprocessing$labels == 'umi_tools_input', ]$labels = 'cutadapt'
df_Preprocessing[df_Preprocessing$labels == 'umi_tools_output', ]$labels = 'umi_tools'
df_Preprocessing[df_Preprocessing$labels == 'fastq_quality_filter_size_ok', ]$labels = 'consume'
df_Preprocessing[df_Preprocessing$labels == 'fastq_quality_filter_seq_to_map', ]$labels = 'fastq_quality_filter'
head(df_Preprocessing); dim(df_Preprocessing)

df_Preprocessing$Type = factor(df_Preprocessing$labels, levels=names(colsBar))

max_y = max(df_Preprocessing[df_Preprocessing$labels == 'fastq_quality_filter', ]$Counts)


print('01_Preprocessing.pdf')
pdf(file=paste0(path_plot, '01_Preprocessing.pdf'))
ggplot(df_Preprocessing, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='dodge2') +
    geom_text(aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    color="black", size=3, angle=90, position = position_dodge(width = .9)) +
    scale_fill_manual(values=colsBar, 'Step') +
    xlab('Library') +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    ggtitle('Number of sequences') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



################################################
# mapping steps

colsMap = c('grey', brewer.pal(n = 8, name = "Dark2")[c(1:4)], 'red')
names(colsMap) = c('unmapped', 'mouse_rRNA', 'human_rRNA', 'mouse_tRNA', 'mouse_cDNA', 'mouse_genome')
colsMap

colsMapType = c('seagreen3', 'palevioletred1', 'cornflowerblue', 'tan3', 'mediumvioletred', 'darkblue', 'seagreen', 'grey')
names(colsMapType) = c('uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric', 'unmapped')
colsMapType

df_STAR_w_cDNA = rbind(STAR_df_simple, STAR_cDNA_df_simple)

df_STAR_w_cDNA$Type = factor(df_STAR_w_cDNA$Type, level= rev(names(colsMapType)))
df_STAR_w_cDNA$Step = factor(df_STAR_w_cDNA$Step, level = names(colsMap))

head(df_STAR_w_cDNA); dim(df_STAR_w_cDNA)


print('02_SequentialMapping.pdf')
pdf(file=paste0(path_plot, '02_SequentialMapping.pdf'))
ggplot() + 
    geom_bar(data=df_STAR_w_cDNA, aes(x=Step, y=Counts, fill=Type), stat='identity') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)') +
    facet_grid(~ factor(Sample)) +
    ylim(NA, max_y)

ggplot() +
    geom_bar(data=df_STAR_w_cDNA, aes(x=Step, y=Counts, fill=Type), stat='identity', position='fill') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)') +
    facet_grid(~ factor(Sample))

###

df_STAR_w_cDNA_overall = rbind(df_STAR_w_cDNA[df_STAR_w_cDNA$Type != 'unmapped',], df_STAR_w_cDNA[df_STAR_w_cDNA$Step == 'mouse_cDNA' & df_STAR_w_cDNA$Type == 'unmapped',])
head(df_STAR_w_cDNA_overall); dim(df_STAR_w_cDNA_overall)

df_STAR_w_cDNA_overall$Step = as.character(df_STAR_w_cDNA_overall$Step)
df_STAR_w_cDNA_overall$Type = as.character(df_STAR_w_cDNA_overall$Type)
df_STAR_w_cDNA_overall[df_STAR_w_cDNA_overall$Step == 'mouse_cDNA' & df_STAR_w_cDNA_overall$Type == 'unmapped', ]$Step = 'unmapped'
head(df_STAR_w_cDNA_overall); dim(df_STAR_w_cDNA_overall)

df_STAR_w_cDNA_overall[df_STAR_w_cDNA_overall$Type == 'uniq_map', ]$Type = 'mapped'
df_STAR_w_cDNA_overall[df_STAR_w_cDNA_overall$Type == 'multi_map', ]$Type = 'mapped'

df_STAR_w_cDNA_overall_aggr = aggregate(df_STAR_w_cDNA_overall$Counts, by=list(df_STAR_w_cDNA_overall$Sample, df_STAR_w_cDNA_overall$Step, df_STAR_w_cDNA_overall$Type), FUN=sum)
names(df_STAR_w_cDNA_overall_aggr) = c('Sample', 'Step', 'Type', 'Counts')
head(df_STAR_w_cDNA_overall_aggr); dim(df_STAR_w_cDNA_overall_aggr)

for_verif_tbl
to_verif = aggregate(df_STAR_w_cDNA_overall_aggr$Counts, by=list(df_STAR_w_cDNA_overall_aggr$Sample), FUN=sum)
names(to_verif) = c('Sample', 'Counts')
to_verif

to_verif_tbl = to_verif$Counts
names(to_verif_tbl) = to_verif$Sample
to_verif_tbl = sort(to_verif_tbl)

if (all(to_verif_tbl == for_verif_tbl) != TRUE){
    print('Problem in the code...')
}

df_STAR_w_cDNA_overall$Type = factor(df_STAR_w_cDNA_overall$Type, levels=c('unmapped', 'mapped'))
df_STAR_w_cDNA_overall$Step = factor(df_STAR_w_cDNA_overall$Step, levels=names(colsMap))

ggplot(df_STAR_w_cDNA_overall, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMap) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)', subtitle = 'Uniquely and multimapped reads') +
    facet_grid(~ factor(Sample), space='free', scale='free') +
    ylim(NA, max_y)

ggplot(df_STAR_w_cDNA_overall, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMap) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)', subtitle = 'Uniquely and multimapped reads') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



###
# 'uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric'
df_STAR_detailed_w_cDNA = rbind(STAR_df_detailed, STAR_cDNA_df_detailed)

df_STAR_detailed_w_cDNA$Type = factor(df_STAR_detailed_w_cDNA$Type, level= rev(names(colsMapType)))
df_STAR_detailed_w_cDNA$Step = factor(df_STAR_detailed_w_cDNA$Step, level = names(colsMap))
head(df_STAR_detailed_w_cDNA); dim(df_STAR_detailed_w_cDNA)


print('02_SequentialMapping_bis.pdf')
pdf(file=paste0(path_plot, '02_SequentialMapping_bis.pdf'))
ggplot() + 
    geom_bar(data=df_STAR_detailed_w_cDNA, aes(x=Step, y=Counts, fill=Type), stat='identity') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)') +
    facet_grid(~ factor(Sample)) +
    ylim(NA, max_y)

ggplot() +
    geom_bar(data=df_STAR_detailed_w_cDNA, aes(x=Step, y=Counts, fill=Type), stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping (w cDNA mapping step)') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



################################################
# STAR_mouse_tRNA unmapped == input for splitseq_to_map by barcode
head(STAR_unmapped_mouse_tRNA_df); dim(STAR_unmapped_mouse_tRNA_df)

max_y_axis_barcode = max(aggregate(STAR_unmapped_mouse_tRNA_df$Counts, by=list(STAR_unmapped_mouse_tRNA_df$Sample), FUN=sum)$x)
max_y_axis_barcode

colsBC = brewer.pal(n = length(unique(STAR_unmapped_mouse_tRNA_df$Step)), name = "Set3")[c(1:length(unique(STAR_unmapped_mouse_tRNA_df$Step)))]
names(colsBC) = unique(STAR_unmapped_mouse_tRNA_df$Step)
colsBC


print('03_SplitBarcode.pdf')
pdf(file=paste0(path_plot, '03_SplitBarcode.pdf'))
ggplot(STAR_unmapped_mouse_tRNA_df, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsBC) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    ggtitle('Split barcode (tRNA_unmapped)') +
    facet_grid(~ factor(Sample), space='free', scale='free') + ylim(c(0, max_y_axis_barcode))

ggplot(STAR_unmapped_mouse_tRNA_df, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsBC) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    ggtitle('Split barcode (tRNA_unmapped)') +
    facet_grid(~ factor(Sample), space='free', scale='free')


ggplot(STAR_unmapped_mouse_tRNA_df, aes(x=Step, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsBC) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Split barcode (tRNA_unmapped)') +
    facet_grid(~ factor(Sample))
dev.off()
###



################################################
# STAR_mouse genome == mapped per barcode
head(STAR_genome_df_simple); dim(STAR_genome_df_simple)
head(STAR_genome_df_detailed); dim(STAR_genome_df_detailed)

STAR_genome_df_simple$Type = factor(STAR_genome_df_simple$Type, levels=rev(names(colsMapType)))
STAR_genome_df_detailed$Type = factor(STAR_genome_df_detailed$Type, levels=rev(names(colsMapType)))

print('04_GenomeMapping.pdf')
pdf(file=paste0(path_plot, '04_GenomeMapping.pdf'))
ggplot(STAR_genome_df_simple, aes(x=Step, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_simple$Type))]) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample))

ggplot(STAR_genome_df_simple, aes(x=Step, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_simple$Type))]) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample))



STAR_genome_df_simple_aggr = aggregate(STAR_genome_df_simple$Counts, by=list(STAR_genome_df_simple$Sample, STAR_genome_df_simple$Type), FUN=sum)
names(STAR_genome_df_simple_aggr) = c('Sample', 'Type', 'Counts')
head(STAR_genome_df_simple_aggr); dim(STAR_genome_df_simple_aggr)

ggplot(STAR_genome_df_simple_aggr, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_simple_aggr$Type))]) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free') +
    ylim(NA, max_y)

ggplot(STAR_genome_df_simple_aggr, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_simple_aggr$Type))]) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



print('04_GenomeMapping_bis.pdf')
pdf(file=paste0(path_plot, '04_GenomeMapping_bis.pdf'))
ggplot(STAR_genome_df_detailed, aes(x=Step, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_detailed$Type))]) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample))

ggplot(STAR_genome_df_detailed, aes(x=Step, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_detailed$Type))]) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample))



STAR_genome_df_detailed_aggr = aggregate(STAR_genome_df_detailed$Counts, by=list(STAR_genome_df_detailed$Sample, STAR_genome_df_detailed$Type), FUN=sum)
names(STAR_genome_df_detailed_aggr) = c('Sample', 'Type', 'Counts')
head(STAR_genome_df_detailed_aggr); dim(STAR_genome_df_detailed_aggr)

ggplot(STAR_genome_df_detailed_aggr, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_detailed_aggr$Type))]) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free') +
    ylim(NA, max_y)

ggplot(STAR_genome_df_detailed_aggr, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType[as.character(unique(STAR_genome_df_detailed_aggr$Type))]) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



################################################
### Genome mapping per barcode
head(STAR_genome_df_simple); dim(STAR_genome_df_simple)

STAR_genome_df_simple_mapped = STAR_genome_df_simple[STAR_genome_df_simple$Type %in% c('uniq_map', 'multi_map'),]
head(STAR_genome_df_simple_mapped)

STAR_genome_df_simple_mapped_aggr = aggregate(STAR_genome_df_simple_mapped$Counts, by=list(STAR_genome_df_simple_mapped$Sample, STAR_genome_df_simple_mapped$Step), FUN=sum)
names(STAR_genome_df_simple_mapped_aggr) = c('Sample', 'Step', 'Counts')
head(STAR_genome_df_simple_mapped_aggr); dim(STAR_genome_df_simple_mapped_aggr)


print('04_GenomeMappingSplitBarcode.pdf')
pdf(file=paste0(path_plot, '04_GenomeMappingSplitBarcode.pdf'))
ggplot(STAR_genome_df_simple_mapped_aggr, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsBC) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free') + ylim(c(0, max_y_axis_barcode))

ggplot(STAR_genome_df_simple_mapped_aggr, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsBC) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    ggtitle('Genome mapping') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



################################################
# duplication rate

dedup_samples = data.frame(matrix(nrow=0, ncol=4))
# for (s in samples_df$Sample){
#     dedup_s = read.csv(paste0(stats_data, s, '_R1_dedup.dat'), header=FALSE, sep='\t')
for (s in samples_df_uniq){
    dedup_s = read.csv(paste0(stats_data, s, '_dedup.dat'), header=FALSE, sep='\t')
    names(dedup_s) = c('Sample', 'Step', 'Type', 'Counts')
    head(dedup_s); dim(dedup_s)
    dedup_samples = rbind(dedup_samples, dedup_s)
    
}
Dedup_df = dedup_samples[dedup_samples$Type == 'duplication_rate', ]
head(Dedup_df); dim(Dedup_df)
Dedup_df$Sample = factor(Dedup_df$Sample, levels=c(unique(Dedup_df$Sample[grep('PF', Dedup_df$Sample)]), unique(Dedup_df$Sample[grep('PD', Dedup_df$Sample)])))

print('05_DuplicationRate.pdf')
pdf(file=paste0(path_plot, '05_DuplicationRate.pdf'))
ggplot(Dedup_df, aes(x=Step, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    geom_text(aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    color="black", size=3, angle=90) +
    scale_fill_manual(values=colsBC) +
    my_theme + theme(legend.position="none",
    axis.title.x=element_blank()) +
    ylab('Rate') +
    ggtitle('Duplication rate (genome mapping)') +
    facet_grid(~ factor(Sample))
dev.off()
###


  
################################################
### Sequential mapping + genome mapping

STAR_genome_df_simple_aggr = aggregate(STAR_genome_df_simple$Counts, by=list(STAR_genome_df_simple$Sample, STAR_genome_df_simple$Type), FUN=sum)
names(STAR_genome_df_simple_aggr) = c('Sample', 'Type', 'Counts')
STAR_genome_df_simple_aggr$Step = 'mouse_genome'
STAR_genome_df_simple_aggr$labels = ''
head(STAR_genome_df_simple_aggr); dim(STAR_genome_df_simple_aggr)

df_STAR_w_genome = rbind(STAR_df_simple, STAR_genome_df_simple_aggr)
head(df_STAR_w_genome); dim(df_STAR_w_genome)

df_STAR_w_genome$Type = factor(df_STAR_w_genome$Type, level= rev(names(colsMapType)))
df_STAR_w_genome$Step = factor(df_STAR_w_genome$Step, level = names(colsMap))

head(df_STAR_w_genome); dim(df_STAR_w_genome)
rownames(df_STAR_w_genome)=NULL


print('06_AllMapping.pdf')
pdf(file=paste0(path_plot, '06_AllMapping.pdf'))
ggplot() + 
    geom_bar(data=df_STAR_w_genome, aes(x=Step, y=Counts, fill=Type), stat='identity') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping') +
    facet_grid(~ factor(Sample)) +
    ylim(NA, max_y)

ggplot() +
    geom_bar(data=df_STAR_w_genome, aes(x=Step, y=Counts, fill=Type), stat='identity', position='fill') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping') +
    facet_grid(~ factor(Sample))

###

df_STAR_w_genome_overall = rbind(df_STAR_w_genome[df_STAR_w_genome$Type != 'unmapped',], df_STAR_w_genome[df_STAR_w_genome$Step == 'mouse_genome' & df_STAR_w_genome$Type == 'unmapped',])
head(df_STAR_w_genome_overall); dim(df_STAR_w_genome_overall)

df_STAR_w_genome_overall$Step = as.character(df_STAR_w_genome_overall$Step)
df_STAR_w_genome_overall$Type = as.character(df_STAR_w_genome_overall$Type)
df_STAR_w_genome_overall[df_STAR_w_genome_overall$Step == 'mouse_genome' & df_STAR_w_genome_overall$Type == 'unmapped', ]$Step = 'unmapped'
head(df_STAR_w_genome_overall); dim(df_STAR_w_genome_overall)

df_STAR_w_genome_overall[df_STAR_w_genome_overall$Type == 'uniq_map', ]$Type = 'mapped'
df_STAR_w_genome_overall[df_STAR_w_genome_overall$Type == 'multi_map', ]$Type = 'mapped'

df_STAR_w_genome_overall_aggr = aggregate(df_STAR_w_genome_overall$Counts, by=list(df_STAR_w_genome_overall$Sample, df_STAR_w_genome_overall$Step, df_STAR_w_genome_overall$Type), FUN=sum)
names(df_STAR_w_genome_overall_aggr) = c('Sample', 'Step', 'Type', 'Counts')
head(df_STAR_w_genome_overall_aggr); dim(df_STAR_w_genome_overall_aggr)

for_verif_tbl
to_verif = aggregate(df_STAR_w_genome_overall_aggr$Counts, by=list(df_STAR_w_genome_overall_aggr$Sample), FUN=sum)
names(to_verif) = c('Sample', 'Counts')
to_verif

to_verif_tbl = to_verif$Counts
names(to_verif_tbl) = to_verif$Sample
to_verif_tbl = sort(to_verif_tbl)

if (all(to_verif_tbl == for_verif_tbl) != TRUE){
    print('Problem in the code...')
}

df_STAR_w_genome_overall$Type = factor(df_STAR_w_genome_overall$Type, levels=c('unmapped', 'mapped'))
df_STAR_w_genome_overall$Step = factor(df_STAR_w_genome_overall$Step, levels=names(colsMap))

ggplot(df_STAR_w_genome_overall, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMap) +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping', subtitle = 'Uniquely and multimapped reads') +
    facet_grid(~ factor(Sample), space='free', scale='free') +
    ylim(NA, max_y)

ggplot(df_STAR_w_genome_overall, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMap) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping', subtitle = 'Uniquely and multimapped reads') +
    facet_grid(~ factor(Sample), space='free', scale='free')
dev.off()
###



###
# 'uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric'
STAR_genome_df_detailed$Step = 'mouse_genome'
df_STAR_detailed_w_genome = rbind(STAR_df_detailed, STAR_genome_df_detailed)


df_STAR_detailed_w_genome$Type = factor(df_STAR_detailed_w_genome$Type, level= rev(names(colsMapType)))
df_STAR_detailed_w_genome$Step = factor(df_STAR_detailed_w_genome$Step, level = names(colsMap))
head(df_STAR_detailed_w_genome); dim(df_STAR_detailed_w_genome)


print('06_AllMapping_bis.pdf')
pdf(file=paste0(path_plot, '06_AllMapping_bis.pdf'))
ggplot() + 
    geom_bar(data=df_STAR_detailed_w_genome, aes(x=Step, y=Counts, fill=Type), stat='identity') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping') +
    facet_grid(~ factor(Sample)) +
    ylim(NA, max_y)

ggplot() +
    geom_bar(data=df_STAR_detailed_w_genome, aes(x=Step, y=Counts, fill=Type), stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom", 
    axis.title.x=element_blank(), legend.key.width= unit(2, "mm"), legend.key.size = unit(1, 'mm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Sequential mapping') +
    facet_grid(~ factor(Sample))
dev.off()
###



################################################
### Transcript types of genome mapping 

maptype_df=data.frame()
sample_BC = unique(paste0(STAR_genome_df$Sample, '_', STAR_genome_df$Step))
sample_BC = sample_BC[c(grep('PF', sample_BC), grep('PD', sample_BC))]

for (sb in sample_BC){
    info = file.info(paste0(stats_data, sb, '_mapType.dat'))
    if (info$size != 0L){
        maptype_s = read.csv(paste0(stats_data, sb, '_mapType.dat'), header=FALSE, sep='\t')
    names(maptype_s) = c('Type', 'Description', 'Counts')
    # maptype_s$Sample = paste(str_split(sb, '_', 6)[[1]][c(1:5)], collapse='_') # library
    # maptype_s$Step = str_split(sb, '_', 6)[[1]][6] # barcode
    maptype_s$Sample = str_split(sb, '_', 2)[[1]][1] # library
    maptype_s$Step = str_split(sb, '_', 2)[[1]][2] # barcode
    maptype_s = maptype_s[c('Sample', 'Step', 'Type', 'Description', 'Counts')]
    head(maptype_s); dim(maptype_s)

    maptype_df = rbind(maptype_df, maptype_s)
    }
}
maptype_df <- maptype_df[maptype_df$Counts > 100,]
maptype_df$Sample = factor(maptype_df$Sample, levels=c(unique(maptype_df$Sample[grep('PF', maptype_df$Sample)]), unique(maptype_df$Sample[grep('PD', maptype_df$Sample)])))
head(maptype_df); dim(maptype_df)

colsTxType = c(brewer.pal(n = 12, name = "Set3")[c(1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12)], 
    'blue', 'green', 'red', 'orange', 'purple', 'pink', 'brown', 'yellow')
names(colsTxType) = sort(unique(maptype_df$Description))

print('07_GenomeMappingTypes.pdf')
pdf(file=paste0(path_plot, '07_GenomeMappingTypes.pdf'), width=10)
ggplot(maptype_df, aes(x=Sample, y=Counts, fill=Description)) +
    geom_bar(stat='identity') +
    #coord_flip() +
    my_theme + theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    scale_fill_manual('', values=colsTxType) +
    ggtitle('Genome mapping - transcript types') +
    facet_grid(~ factor(Sample), space='free', scale='free')

ggplot(maptype_df, aes(x=Step, y=Counts, fill=Description)) +
    geom_bar(stat='identity') +
    #coord_flip() +
    my_theme + theme(axis.title.x=element_blank()) +
    scale_fill_manual('', values=colsTxType) +
    ggtitle('Genome mapping - transcript types') +
    facet_grid(~ factor(Sample))
dev.off()




################################################
### read length distribution - trimmed

trimmed_readLen_df=data.frame()
for (s in samples_df_uniq){
    trimmed_readLen_s = read.csv(paste0(stats_data, s, '_trimmed_readLenDist.dat'), header=FALSE, sep='\t')
    names(trimmed_readLen_s) = c('ReadLen', 'Counts')
    trimmed_readLen_s$Sample = s
    trimmed_readLen_s = trimmed_readLen_s[c('Sample', 'ReadLen', 'Counts')]
    head(trimmed_readLen_s); dim(trimmed_readLen_s)

    trimmed_readLen_df = rbind(trimmed_readLen_df, trimmed_readLen_s)

}
head(trimmed_readLen_df); dim(trimmed_readLen_df)
trimmed_readLen_df$Ribosome=''
trimmed_readLen_df[grep('^PF', trimmed_readLen_df$Sample),]$Ribosome = 'Monosome'
if (length(grep('^PD', trimmed_readLen_df$Sample)) > 0){
    trimmed_readLen_df[grep('^PD', trimmed_readLen_df$Sample),]$Ribosome = 'Disome'
}


print('08_ReadLen_trimmed.pdf')
pdf(file=paste0(path_plot, '08_ReadLen_trimmed.pdf'))
for (ribo in unique(trimmed_readLen_df$Ribosome)){
    ribo_df = trimmed_readLen_df[trimmed_readLen_df$Ribosome == ribo,]

    scale_x = scale_x_continuous(breaks=seq(min(ribo_df$ReadLen), max(ribo_df$ReadLen), 3))

    ribo_plot = ggplot(ribo_df, aes(x=ReadLen, y=Counts)) + 
    geom_bar(stat='identity', fill='royalblue') +
    scale_x +
    facet_wrap(. ~ factor(Sample),
                ncol=ceiling(length(unique(trimmed_readLen_df$Sample))/2)) + 
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - trimmed')

    plot(ribo_plot)
    }
dev.off()




################################################
### read length distribution - filtered

filtered_readLen_df=data.frame()
for (s in samples_df_uniq){
    filtered_readLen_s = read.csv(paste0(stats_data, s, '_filtered_readLenDist.dat'), header=FALSE, sep='\t')
    names(filtered_readLen_s) = c('ReadLen', 'Counts')
    filtered_readLen_s$Sample = s
    filtered_readLen_s = filtered_readLen_s[c('Sample', 'ReadLen', 'Counts')]
    head(filtered_readLen_s); dim(filtered_readLen_s)

    filtered_readLen_df = rbind(filtered_readLen_df, filtered_readLen_s)

}
head(filtered_readLen_df); dim(filtered_readLen_df)
filtered_readLen_df$Ribosome=''
filtered_readLen_df[grep('^PF', filtered_readLen_df$Sample),]$Ribosome = 'Monosome'
if (length(grep('^PD', filtered_readLen_df$Sample)) > 0){
    filtered_readLen_df[grep('^PD', filtered_readLen_df$Sample),]$Ribosome = 'Disome'
}


print('09_ReadLen_filtered.pdf')
pdf(file=paste0(path_plot, '09_ReadLen_filtered.pdf'))
for (ribo in rev(unique(filtered_readLen_df$Ribosome))){
    ribo_df = filtered_readLen_df[filtered_readLen_df$Ribosome == ribo,]

    scale_x = scale_x_continuous(breaks=seq(min(ribo_df$ReadLen), max(ribo_df$ReadLen), 3))

    ribo_plot = ggplot(ribo_df, aes(x=ReadLen, y=Counts)) + 
    geom_bar(stat='identity', fill='royalblue') +
    scale_x +
    facet_wrap(. ~ factor(Sample),
                ncol=ceiling(length(unique(filtered_readLen_df$Sample))/2)) + 
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - filtered')

    plot(ribo_plot)
    }
dev.off()




################################################
### read length distribution - mappedTx PER BARCODE

mappedTx_readLen_df=data.frame()
sample_BC = unique(paste0(STAR_genome_df$Sample, '_', STAR_genome_df$Step))
sample_BC = sample_BC[c(grep('PF', sample_BC), grep('PD', sample_BC))]

for (sb in sample_BC){
    print(sb)
    mappedTx_s = read.csv(paste0(stats_data, sb, '_mappedTx_readLenDist.dat'), header=FALSE, sep='\t')
    names(mappedTx_s) = c('ReadLen', 'Counts')
    # mappedTx_s$Sample = paste(str_split(sb, '_', 6)[[1]][c(1:5)], collapse='_') # library
    # mappedTx_s$Step = str_split(sb, '_', 6)[[1]][6] # barcode
    mappedTx_s$Sample = str_split(sb, '_', 2)[[1]][1] # library
    mappedTx_s$Step = str_split(sb, '_', 2)[[1]][2] # barcode
    mappedTx_s = mappedTx_s[c('Sample', 'Step', 'ReadLen', 'Counts')]
    head(mappedTx_s); dim(mappedTx_s)

    mappedTx_readLen_df = rbind(mappedTx_readLen_df, mappedTx_s)

}
head(mappedTx_readLen_df); dim(mappedTx_readLen_df)
mappedTx_readLen_df$Ribosome=''
mappedTx_readLen_df[grep('^PF', mappedTx_readLen_df$Sample),]$Ribosome = 'Monosome'
if (length(grep('^PD', mappedTx_readLen_df$Sample)) > 0){
    mappedTx_readLen_df[grep('^PD', mappedTx_readLen_df$Sample),]$Ribosome = 'Disome'
}

mappedTx_readLen_df = mappedTx_readLen_df[mappedTx_readLen_df$ReadLen != 0,]

print('10_ReadLen_onTranscriptome.pdf')
pdf(file=paste0(path_plot, '10_ReadLen_onTranscriptome.pdf'), height=9, width=12)
for (ribo in unique(mappedTx_readLen_df$Ribosome)){
    ribo_df = mappedTx_readLen_df[mappedTx_readLen_df$Ribosome == ribo,]

    scale_x = scale_x_continuous(breaks=seq(min(ribo_df$ReadLen), max(ribo_df$ReadLen), 3))

    ribo_plot = ggplot(ribo_df, aes(x=ReadLen, y=Counts)) + 
    geom_bar(stat='identity', fill='royalblue') +
    scale_x +
    facet_wrap(. ~ factor(Sample)+factor(Step),
                labeller = label_wrap_gen(multi_line=FALSE),
                ncol=ceiling(length(unique(paste0(mappedTx_readLen_df$Sample, '_', mappedTx_readLen_df$Step)))/4)) + 
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - projected on Transcriptome')

    plot(ribo_plot)
    }
dev.off()



################################################
### RSEM total counts
RSEM_df=data.frame()

for (s in samples_df_uniq){
    print(s)
    RSEM_s = read.csv(paste0(stats_data, s, '_RSEM.dat'), header=FALSE, sep='\t')
    names(RSEM_s) = c('Sample', 'Step', 'Type', 'Counts')
    head(RSEM_s); dim(RSEM_s)

    RSEM_df = rbind(RSEM_df, RSEM_s)

}
head(RSEM_df); dim(RSEM_df)
RSEM_df$Ribosome=''
RSEM_df[grep('^PF', RSEM_df$Sample),]$Ribosome = 'Monosome'
if (length(grep('^PD', RSEM_df$Sample)) > 0){
    RSEM_df[grep('^PD', RSEM_df$Sample),]$Ribosome = 'Disome'
}
RSEM_df = RSEM_df[!is.na(RSEM_df$Counts),]
head(RSEM_df); dim(RSEM_df)


print('11_RSEM_TotalExpectedCounts.pdf')
pdf(file=paste0(path_plot, '11_RSEM_TotalExpectedCounts.pdf'), height=9, width=12)
RSEM_plot = ggplot(RSEM_df, aes(x=Step, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    # geom_text(aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    # color="black", size=3, angle=90) +
    geom_text(data = RSEM_df[RSEM_df$Counts > mean(RSEM_df$Counts), ], aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    color="black", size=3, angle=90) +
    geom_text(data = RSEM_df[RSEM_df$Counts <= mean(RSEM_df$Counts), ], aes(y=Counts, label=Counts), vjust=0.5, hjust=-0.1, 
    color="black", size=3, angle=90) +
    scale_fill_manual(values=colsBC) +
    my_theme + theme(legend.position="none",
    axis.title.x=element_blank()) +
    ggtitle('RSEM total expected counts') +
    facet_grid(~ factor(Sample))

plot(RSEM_plot)  
dev.off()


################################################
