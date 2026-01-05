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

id=$1
work_dir=$2


raw_data=${work_dir}raw_data/
trimmed_data=${work_dir}trimmed_data/
RSEM_data=${work_dir}RSEM/

stats_data=${work_dir}stats/

echo $id

###############################################################################
echo "Calculating pre-processing stats"

stats_out=${stats_data}${id}_stats.dat

# raw_data - done on curnagl
#seqkit stats --tabular --threads 10 ${id}_${run_nb}_${lane_nb}_001.fastq.gz
# >> QC_fastq.sh

# concatenated raw data - done on curnagl
#seqkit stats --tabular --threads 10 ${id}.fastq.gz
# >> QC_fastq_suite.sh
# >> QC_fast_concat.sh (after QC_fastq.sh and QC_fastq_suite.sh)


###############################################################################
trimmed_log=${trimmed_data}${id}_trimming.log
# identical info in ${id}.fastq.gz_trimming_report.txt

# total nb reads
total_reads=$(awk '$0~"Total reads processed:"{print $4}' $trimmed_log | sed 's/,//g')

# nb reads with adapters
reads_w_adapter=$(awk '$0~"Reads with adapters:"{print $4}' $trimmed_log | sed 's/,//g')

# nb reads okay after adapter trimming and quality filtering
reads_okay=$(awk '$0~"Reads written"{print $5}' $trimmed_log | sed 's/,//g')

# nb reads removed because too short
reads_too_short=$(awk '$0~"Sequences removed"{print $14}' $trimmed_log | sed 's/,//g')

# final nb reads okay (present in ${id}_trimmed.fq.gz)
reads_final_ok=$(echo "$reads_okay - $reads_too_short" | bc)


# already created when processing raw_data
echo -e "${id}\ttrim_galore\tinput\t${total_reads}" > $stats_out
echo -e "${id}\ttrim_galore\tw_adapters\t${reads_w_adapter}" >> $stats_out
echo -e "${id}\ttrim_galore\ttoo_short\t${reads_too_short}" >> $stats_out
echo -e "${id}\ttrim_galore\toutput\t${reads_final_ok}" >> $stats_out



###############################################################################
echo "Calculating mapping stats"
mapped_log=${RSEM_data}${id}/mouse_genome.Log.final.out

# input nb reads okay == final nb reads okay after trim_galore
reads_final_ok=$(awk '$0 ~ "Number of input reads"{print $NF}' $mapped_log)

# nb reads uniquely mapped
uniq_map=$(awk '$0 ~ "Uniquely mapped reads number"{print $NF}' $mapped_log)

# nb reads multimapped
multi_map=$(awk '$0 ~ "Number of reads mapped to multiple loci"{print $NF}' $mapped_log)

# nb reads unmapped == too many loci + too many mismatches + too short + other
unmapped=$(echo "$reads_final_ok - $uniq_map - $multi_map" | bc)

#
# nb reads mapped to too many loci
too_many_loci=$(awk '$0 ~ "Number of reads mapped to too many loci"{print $NF}' $mapped_log)

# nb reads unmapped because too many mismatches
too_many_mismatches=$(awk '$0 ~ "Number of reads unmapped: too many mismatches"{print $NF}' $mapped_log)

# too short unmapped reads
too_short=$(awk '$0 ~ "Number of reads unmapped: too short"{print $NF}' $mapped_log)

# 'other' unmapped reads
other=$(awk '$0 ~ "Number of reads unmapped: other"{print $NF}' $mapped_log)

# chimeric reads
chimeric=$(awk '$0 ~ "Number of chimeric reads"{print $NF}' $mapped_log)
#

echo -e "${id}\tSTAR_RSEM\tuniq_map\t${uniq_map}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\tmulti_map\t${multi_map}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\tunmapped\t${unmapped}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\ttoo_many_loci\t${too_many_loci}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\ttoo_short\t${too_short}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\tother\t${other}" >> $stats_out
echo -e "${id}\tSTAR_RSEM\tchimeric\t${chimeric}" >> $stats_out


###############################################################################
echo "Calculation read length distributions"

# concatenated raw_data - done on curnagl
#zcat ${id}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${stats_data}${id}_raw_readLenDist.tmp
#awk '{printf("%s\t%s\n", $2, $1)}' ${stats_data}${id}_raw_readLenDist.tmp > ${stats_data}${id}_raw_readLenDist.dat
#rm ${stats_data}${id}_raw_readLenDist.tmp

# trimmed_data
echo "   trimmed_data"

zcat ${trimmed_data}${id}_trimmed.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${stats_data}${id}_trimmed_readLenDist.tmp
awk '{printf("%s\t%s\n", $2, $1)}' ${stats_data}${id}_trimmed_readLenDist.tmp > ${stats_data}${id}_trimmed_readLenDist.dat
rm ${stats_data}${id}_trimmed_readLenDist.tmp


# BAM
BAM=${work_dir}RSEM/${id}/mouse_genome.Aligned.out.bam
BAM_Tx=${work_dir}RSEM/${id}/mouse_genome.Aligned.toTranscriptome.out.bam
# BAM_Tx=${work_dir}RSEM/${id}/mouse_genome.transcript.bam
# samtools flagstat mouse_genome.Aligned.toTranscriptome.out.bam == samtools flagstat mouse_genome.transcript.bam
# 

#echo "   mapped to mouse_genome"
#readLengthDistribution.py $BAM > ${stats_data}${id}_mapped_readLenDist.dat
# same results as ${stats_data}${id}_trimmed_readLenDist.dat

bioawk=/home/vricci/miniconda3/bin/bioawk
echo "   projected to mouse_Tx"

#readLengthDistribution.py $BAM_Tx > ${stats_data}${id}_mappedTx_readLenDist.dat
samtools view $BAM_Tx | $bioawk -c sam '{hist[length($seq)]++} END {for (l in hist) print l, hist[l]}' | sort -n -k1 > ${stats_data}${id}_mappedTx_readLenDist.dat


###############################################################################
echo "Calculating gene and transcripts map types"

map_out=${stats_data}${id}_mapType.dat
# mapType_Tx.py for BAM with only Tx and not Gene|Tx
GENEBED='/data/databases/mouse/Mmusculus.GRCm39.111.startStop.IDS.bed'

samtools view $BAM_Tx | ${work_dir}mapType_Tx.py -b $GENEBED -c 100000 > $map_out



###############################################################################
echo 'RSEM total counts'
RSEM_out=${stats_data}${id}_RSEM.dat

RSEM_genes=${work_dir}RSEM/${id}/mouse_genome.genes.results
RSEM_isoforms=${work_dir}RSEM/${id}/mouse_genome.isoforms.results

total_expected_count_genes=$(cut -f5 $RSEM_genes | tail -n +2 | awk '{s+=$1}END{print s}')
total_expected_count_isoforms=$(cut -f5 $RSEM_isoforms | tail -n +2 | awk '{s+=$1}END{print s}')

echo -e "${id}\t${bc}\texpected_count\t${total_expected_count_genes}" > $RSEM_out
# echo -e "${id}\t${bc}\texpected_count\t${total_expected_count_isoforms}" >> $RSEM_out
# both are the same value!

# DONE - done to avoid re-running all pipeline!
# cd /data/vricci/projects/Lisa/RNAseq_mESC_degron/Extra/RSEM
# for id in $(ls | grep '[0-9].*$'); do 
#     echo $id
#     RSEM_out=/data/vricci/projects/Lisa/RNAseq_mESC_degron/Extra/stats/${id}_RSEM.dat
#     RSEM_genes=${id}/mouse_genome.genes.results
#     RSEM_isoforms=${id}/mouse_genome.isoforms.results

#     total_expected_count_genes=$(cut -f5 $RSEM_genes | tail -n +2 | awk '{s+=$1}END{print s}')
#     total_expected_count_isoforms=$(cut -f5 $RSEM_isoforms | tail -n +2 | awk '{s+=$1}END{print s}')

#     echo -e "${id}\t${bc}\texpected_count\t${total_expected_count_genes}" > $RSEM_out


# done
