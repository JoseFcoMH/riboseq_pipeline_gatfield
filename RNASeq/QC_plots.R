library(ggplot2)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(RColorBrewer)
library(ggpubr)

args<-commandArgs(TRUE)

work_dir=args[1]
# work_dir = '/data/vricci/projects/Lisa/RNAseq_mESC_degron/Extra/'
raw_data=paste0(work_dir, 'raw_data/')
stats_data=paste0(work_dir, 'stats/')
path_plot=paste0(stats_data, '/plots/')
system(paste0('mkdir -p ', path_plot))

samples_df = read.csv(list.files(path=raw_data, pattern='_data_uniq.txt', full.names=T), header=FALSE, sep='\t')
names(samples_df) = 'Sample'
head(samples_df); dim(samples_df)

# Stats
stats_df = read.csv(paste0(stats_data, 'Stats.txt'), header=FALSE, sep='\t')
names(stats_df) = c('Sample', 'Step', 'Type', 'Counts')
head(stats_df); dim(stats_df)

# concatenated samples only >> 11Wdm4h, 12Wdm4h, etc...
stats_samples = stats_df[stats_df$Sample %in% samples_df$Sample,]
head(stats_samples); dim(stats_samples)

for (s in samples_df$Sample){
    stats_s = read.csv(paste0(stats_data, s, '_stats.dat'), header=FALSE, sep='\t')
    names(stats_s) = c('Sample', 'Step', 'Type', 'Counts')
    head(stats_s); dim(stats_s)
    stats_samples = rbind(stats_samples, stats_s)
    
}
head(stats_samples)
# nums_seq == input
# output == input - too_short
# output == uniq_map + multi_map + unmapped

my_theme=theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
    strip.text.y.right = element_text(angle = 0,  size=5), 
    strip.background = element_rect(colour="black", fill="white"))

colsType = c("#377EB8", "#FFFF33", "#984EA3", "#FF7F00", 'seagreen3', 'palevioletred1', 'grey')
names(colsType) = c('input', 'w_adapters', 'too_short', 'output', 'uniq_map', 'multi_map', 'unmapped')



################################################
# Plot Z
# input
df_Z = stats_samples[stats_samples$Type == 'input',]

plot_Z = ggplot(df_Z, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    geom_text(aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    color="black", size=3, angle=90) +
    scale_fill_manual(values=colsType) +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Raw sequencing reads')

# Plot A
# input = output + too_short
df_A = stats_samples[stats_samples$Type %in% c('output', 'too_short'),]
head(df_A); dim(df_A)

df_A$Type = factor(df_A$Type, levels=c('too_short', 'output'))

plot_A = ggplot(df_A, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsType) +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('TrimGalore sequencing reads')

plot_a = ggplot(df_A, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsType) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('TrimGalore sequencing reads')


# Plot B
# output == uniq_map + multi_map + unmapped
df_B = stats_samples[stats_samples$Type %in% c('uniq_map', 'multi_map', 'unmapped'),]
df_B$Type = factor(df_B$Type, levels=c('unmapped', 'multi_map', 'uniq_map'))
head(df_B); dim(df_B)

plot_B = ggplot(df_B, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsType) +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping')

plot_b = ggplot(df_B, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsType) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping')


# output == uniq_map + multi_map + too_many_loci + too_many_mismatches + too_short + other + chimeric'
kinds = c('uniq_map', 'multi_map', 'too_many_loci', 'too_many_mismatches', 'too_short', 'other', 'chimeric')
df_B = stats_samples[stats_samples$Type %in% kinds,]
df_B$Type = factor(df_B$Type, levels=rev(kinds))

head(df_B); dim(df_B)

colsMapType = c('seagreen3', 'palevioletred1', 'cornflowerblue', 'tan3', 'mediumvioletred', 'darkblue', 'seagreen')
names(colsMapType) = kinds

plot_BB = ggplot(df_B, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=colsMapType) +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping')

plot_bb = ggplot(df_B, aes(x=Sample, y=Counts, fill=Type)) +
    geom_bar(stat='identity', position='fill') +
    scale_y_continuous(labels=seq(0, 100, 25), breaks=seq(0, 1, 0.25)) +
    scale_fill_manual(values=colsMapType) +
    ylab('Percent') +
    my_theme + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    ggtitle('Genome mapping')

pdf(file=paste0(path_plot, '01_Stats.pdf'))
plot_Z
plot_A
plot_a
dev.off()

pdf(file=paste0(path_plot, '02_GenomeMapping.pdf'))
plot_B
plot_b
dev.off()

pdf(file=paste0(path_plot, '02_GenomeMapping_bis.pdf'))
plot_BB
plot_bb
dev.off()



################################################
# mapType

maptype_df=data.frame()
for (s in samples_df$Sample){
    maptype_s = read.csv(paste0(stats_data, s, '_mapType.dat'), header=FALSE, sep='\t')
    names(maptype_s) = c('Type', 'Description', 'Counts')
    maptype_s$Sample = s
    maptype_s = maptype_s[c('Sample', 'Type', 'Description', 'Counts')]
    head(maptype_s); dim(maptype_s)

    maptype_df = rbind(maptype_df, maptype_s)

}
maptype_df <- maptype_df[maptype_df$Counts > 100,]
head(maptype_df); dim(maptype_df)

colsTxType = c(brewer.pal(n = 12, name = "Set3")[c(1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12)], 'blue', 'green', 'red')
names(colsTxType) = sort(unique(maptype_df$Description))

# plot_C
# mapping types
plot_C <- ggplot(maptype_df,
              aes(x=Sample, y=Counts, fill=Description)) +
    geom_bar(stat='identity') +
    #coord_flip() +
    my_theme + theme(axis.title.x=element_blank()) +
    scale_fill_manual('', values=colsTxType) +
    ggtitle('Genome mapping - transcript types')

pdf(file=paste0(path_plot, '03_GenomeMappingTypes.pdf'), width=10)
plot_C
dev.off()



################################################
# read length distribution - trimmed reads
trimmed_readLen_df=data.frame()
for (s in samples_df$Sample){
    trimmed_readLen_s = read.csv(paste0(stats_data, s, '_trimmed_readLenDist.dat'), header=FALSE, sep='\t')
    names(trimmed_readLen_s) = c('ReadLen', 'Counts')
    trimmed_readLen_s$Sample = s
    trimmed_readLen_s = trimmed_readLen_s[c('Sample', 'ReadLen', 'Counts')]
    head(trimmed_readLen_s); dim(trimmed_readLen_s)

    trimmed_readLen_df = rbind(trimmed_readLen_df, trimmed_readLen_s)

}
head(trimmed_readLen_df); dim(trimmed_readLen_df)

# plot_D
# trimmed readLen distribution
plot_D <- ggplot(trimmed_readLen_df,
              aes(x=ReadLen, y=Counts)) +
    geom_bar(stat='identity', fill='royalblue') +
    facet_wrap(~ factor(Sample),
               ncol=ceiling(length(unique(trimmed_readLen_df$Sample))/4)) + 
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - trimmed')

plot_d <- ggplot(trimmed_readLen_df,
              aes(x=ReadLen, y=Counts)) +
    geom_bar(stat='identity', fill='royalblue') +
    facet_wrap(~ factor(Sample),
               ncol=ceiling(length(unique(trimmed_readLen_df$Sample))/4)) + 
    scale_x_continuous(limits=c(140, max(trimmed_readLen_df$ReadLen)),
        breaks=seq(140, max(trimmed_readLen_df$ReadLen), 1)) +
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - trimmed')

pdf(file=paste0(path_plot, '04_ReadLen_trimmed.pdf'), height=9, width=12)
plot_D
plot_d
dev.off()



################################################
# read length distribution - mappedTx
mappedTx_readLen_df=data.frame()
for (s in samples_df$Sample){
    mappedTx_readLen_s = read.csv(paste0(stats_data, s, '_mappedTx_readLenDist.dat'), header=FALSE, sep='\t')
    names(mappedTx_readLen_s) = c('ReadLen', 'Counts')
    mappedTx_readLen_s$Sample = s
    mappedTx_readLen_s = mappedTx_readLen_s[c('Sample', 'ReadLen', 'Counts')]
    head(mappedTx_readLen_s); dim(mappedTx_readLen_s)

    mappedTx_readLen_df = rbind(mappedTx_readLen_df, mappedTx_readLen_s)

}
head(mappedTx_readLen_df); dim(mappedTx_readLen_df)

# plot_E
# trimmed readLen distribution
plot_E <- ggplot(mappedTx_readLen_df,
              aes(x=ReadLen, y=Counts)) +
    geom_bar(stat='identity', fill='royalblue') +
    facet_wrap(~ factor(Sample),
               ncol=ceiling(length(unique(mappedTx_readLen_df$Sample))/4)) + 
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - projected on Transcriptome')


plot_e <- ggplot(mappedTx_readLen_df,
              aes(x=ReadLen, y=Counts)) +
    geom_bar(stat='identity', fill='royalblue') +
    facet_wrap(~ factor(Sample),
               ncol=ceiling(length(unique(mappedTx_readLen_df$Sample))/4)) + 
    scale_x_continuous(limits=c(140, max(mappedTx_readLen_df$ReadLen)),
        breaks=seq(140, max(mappedTx_readLen_df$ReadLen), 1)) +
    my_theme + theme(axis.title.x=element_blank()) +
    ggtitle('Read length distribution - projected on Transcriptome')

pdf(file=paste0(path_plot, '05_ReadLen_onTranscriptome.pdf'), height=9, width=12)
plot_E
plot_e
dev.off()



################################################
### RSEM total counts
RSEM_df=data.frame()
for (s in samples_df$Sample){
    RSEM_s = read.csv(paste0(stats_data, s, '_RSEM.dat'), header=FALSE, sep='\t')
    names(RSEM_s) = c('Sample', 'Step', 'Type', 'Counts')
    head(RSEM_s); dim(RSEM_s)

    RSEM_df = rbind(RSEM_df, RSEM_s)

}
head(RSEM_df); dim(RSEM_df)
RSEM_df = RSEM_df[!is.na(RSEM_df$Counts),]
head(RSEM_df); dim(RSEM_df)


print('06_RSEM_TotalExpectedCounts.pdf')
pdf(file=paste0(path_plot, '06_RSEM_TotalExpectedCounts.pdf'), height=9, width=12)
RSEM_plot = ggplot(RSEM_df, aes(x=Sample, y=Counts, fill=Step)) +
    geom_bar(stat='identity') +
    # geom_text(aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    # color="black", size=3, angle=90) +
    geom_text(data = RSEM_df[RSEM_df$Counts >= 10000, ], aes(y=Counts, label=Counts), vjust=0.5, hjust=1.1, 
    color="black", size=3, angle=90) +
    geom_text(data = RSEM_df[RSEM_df$Counts < 10000, ], aes(y=Counts, label=Counts), vjust=0.5, hjust=-0.1, 
    color="black", size=3, angle=90) +
    my_theme + theme(legend.position="none",
    axis.title.x=element_blank()) +
    ggtitle('RSEM total expected counts')

plot(RSEM_plot)  
dev.off()

