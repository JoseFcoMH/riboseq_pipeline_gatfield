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

# work_dir='/data/vricci/SnakemakeTest/'
# id='PF73_2611_1_001_R1'


raw_data=${work_dir}raw_data/
umi_data=${work_dir}umi_data/
trimmed_data=${work_dir}trimmed_data/
filtered_data=${work_dir}filtered_data/
mapping_data=${workdir}mapping_data/STAR/
RSEM_data=${work_dir}RSEM/

stats_data=${work_dir}stats/
mkdir -p $stats_data

echo $id

###############################################################################
echo "Calculating raw_data stats"

### raw_data
seqkit stats --tabular --threads 10 ${raw_data}${id}.fastq.gz > ${stats_data}/stats_${id}.tmp

num_seqs=$(cut -d$'\t' -f4 ${stats_data}/stats_${id}.tmp | sed -n 2p)
sum_len=$(cut -d$'\t' -f5 ${stats_data}/stats_${id}.tmp | sed -n 2p)
min_len=$(cut -d$'\t' -f6 ${stats_data}/stats_${id}.tmp | sed -n 2p)
avg_len=$(cut -d$'\t' -f7 ${stats_data}/stats_${id}.tmp | sed -n 2p)
max_len=$(cut -d$'\t' -f8 ${stats_data}/stats_${id}.tmp | sed -n 2p)

stats_out=${stats_data}/${id}_stats.txt
stats_len=${stats_data}/${id}_stats_len.txt


echo -e "${id}\traw_fastq\tinput\t${num_seqs}" > $stats_out

echo -e "${id}\traw_fastq\tsum_len\t${sum_len}" > $stats_len
echo -e "${id}\traw_fastq\tmin_len\t${min_len}" >> $stats_len
echo -e "${id}\traw_fastq\tavg_len\t${avg_len}" >> $stats_len
echo -e "${id}\traw_fastq\tmax_len\t${max_len}" >> $stats_len

rm ${stats_data}/stats_${id}.tmp
###


###############################################################################
echo "Calculating barcode stats"

# whitelist.log
whitelist_log=${umi_data}${id}_whitelist.log

# total nb reads == nb trimmed reads!
total_reads=$(awk '$0~"reads matched the barcode pattern"{print $4}' $whitelist_log)

# nb of unique barcodes
nb_uniq_BC=$(awk '$0~"unique cell barcodes"{print $5}' $whitelist_log)

# nb reads matching selected barcodes
w_barcodes=$(awk '$0~"total reads matching the selected cell barcodes"{print $5}' $whitelist_log)

# nb reads that can be corrected
to_correct=$(awk '$0~"total reads which can be error corrected to the selected cell barcodes"{print $5}' $whitelist_log)
# to correct == corrected



###############################################################################
echo "Calculating pre-processing stats"

# stats_out=${stats_data}${id}_stats.dat

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
echo -e "${id}\tcutadapt\tinput\t${total_reads}" >> $stats_out
echo -e "${id}\tcutadapt\tw_adapters\t${reads_w_adapter}" >> $stats_out
echo -e "${id}\tcutadapt\ttoo_short\t${reads_too_short}" >> $stats_out
echo -e "${id}\tcutadapt\toutput\t${reads_final_ok}" >> $stats_out


###############################################################################
# extract.log
extract_log=${umi_data}${id}_extract.log

# input nb reads
input_reads=$(awk '$0~"Input Reads"{print $6}' $extract_log)

# output nb reads
output_reads=$(awk '$0~"Reads output"{print $6}' $extract_log)

# nb reads not correctable
not_corrected=$(awk '$0~"Filtered cell barcode. Not correctable"{print $9}' $extract_log)

# nb reads corrected
corrected=$(awk '$0~"False cell barcode. Error-corrected"{print $8}' $extract_log)

echo -e "${id}\tumi_tools\tinput\t${input_reads}" >> $stats_out
echo -e "${id}\tumi_tools\tuniq_BC\t${nb_uniq_BC}" >> $stats_out
echo -e "${id}\tumi_tools\tnot_corrected\t${not_corrected}" >> $stats_out
echo -e "${id}\tumi_tools\tcorrected\t${corrected}" >> $stats_out
echo -e "${id}\tumi_tools\toutput\t${output_reads}" >> $stats_out



###############################################################################
size_log=${filtered_data}${id}_size.log

filtered_out=$(awk '$0~"Number of sequences filtered out"{print $6}' $size_log)

echo -e "${id}\tconsume\tfiltered_out\t${filtered_out}" >> $stats_out

qual_log=${filtered_data}${id}_qual.log

okay_size=$(awk '$0~"Input"{print $2}' $qual_log)
low_qual=$(awk '$0~"discarded"{print $2}' $qual_log)
reads_to_map=$(awk '$0~"Output"{print $2}' $qual_log)

echo -e "${id}\tfastq_quality_filter\tsize_ok\t${okay_size}" >> $stats_out
echo -e "${id}\tfastq_quality_filter\tlow_qual\t${low_qual}" >> $stats_out
echo -e "${id}\tfastq_quality_filter\tseq_to_map\t${reads_to_map}" >> $stats_out



###############################################################################
echo "Calculating mapping stats"

echo "   mouse_rRNA"
# mouse_rRNA
mapped_log=${mapping_data}${id}/mouse_rRNA.Log.final.out

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

echo -e "${id}\tSTAR_mouse_rRNA\tuniq_map\t${uniq_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\tmulti_map\t${multi_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\tunmapped\t${unmapped}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\ttoo_many_loci\t${too_many_loci}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\ttoo_short\t${too_short}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\tother\t${other}" >> $stats_out
echo -e "${id}\tSTAR_mouse_rRNA\tchimeric\t${chimeric}" >> $stats_out


echo "   human_rRNA"
# human_rRNA
mapped_log=${mapping_data}${id}/human_rRNA.Log.final.out

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

echo -e "${id}\tSTAR_human_rRNA\tuniq_map\t${uniq_map}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\tmulti_map\t${multi_map}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\tunmapped\t${unmapped}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\ttoo_many_loci\t${too_many_loci}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\ttoo_short\t${too_short}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\tother\t${other}" >> $stats_out
echo -e "${id}\tSTAR_human_rRNA\tchimeric\t${chimeric}" >> $stats_out



echo "   mouse_tRNA"
# mouse_tRNA
mapped_log=${mapping_data}${id}/mouse_tRNA.Log.final.out

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

echo -e "${id}\tSTAR_mouse_tRNA\tuniq_map\t${uniq_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\tmulti_map\t${multi_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\tunmapped\t${unmapped}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\ttoo_many_loci\t${too_many_loci}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\ttoo_short\t${too_short}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\tother\t${other}" >> $stats_out
echo -e "${id}\tSTAR_mouse_tRNA\tchimeric\t${chimeric}" >> $stats_out


echo "   mouse_cDNA"
# mouse_cDNA
mapped_log=${mapping_data}${id}/mouse_cDNA.Log.final.out

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

echo -e "${id}\tSTAR_mouse_cDNA\tuniq_map\t${uniq_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\tmulti_map\t${multi_map}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\tunmapped\t${unmapped}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\ttoo_many_loci\t${too_many_loci}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\ttoo_short\t${too_short}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\tother\t${other}" >> $stats_out
echo -e "${id}\tSTAR_mouse_cDNA\tchimeric\t${chimeric}" >> $stats_out



###############################################################################
echo "Calculating group_split stats"

flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")
    # pattern between 'out.' and '.mate1' >> mouse_tRNA.Unmapped.out.AGCTA.mate1
    num_seqs=$(echo $line | cut -d$' ' -f4)
    sum_len=$(echo $line | cut -d$' ' -f5)
    min_len=$(echo $line | cut -d$' ' -f6)
    avg_len=$(echo $line | cut -d$' ' -f7)
    max_len=$(echo $line | cut -d$' ' -f8)

    echo -e "${id}\t${bc}\tsplitseq_to_map\t${num_seqs}" >> $stats_out

    echo -e "${id}\t${bc}\tsum_len\t${sum_len}" >> $stats_len
    echo -e "${id}\t${bc}\tmin_len\t${min_len}" >> $stats_len
    echo -e "${id}\t${bc}\tavg_len\t${avg_len}" >> $stats_len
    echo -e "${id}\t${bc}\tmax_len\t${max_len}" >> $stats_len

done
###


###############################################################################
echo 'Calculating mapping stats STAR_RSEM'
echo "   mouse_genome"

flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")

    # mouse_genome
    mapped_log=${RSEM_data}${id}/${bc}/mouse_genome.Log.final.out

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

    echo -e "${id}\t${bc}\tuniq_map\t${uniq_map}" >> $stats_out
    echo -e "${id}\t${bc}\tmulti_map\t${multi_map}" >> $stats_out
    echo -e "${id}\t${bc}\tunmapped\t${unmapped}" >> $stats_out
    echo -e "${id}\t${bc}\ttoo_many_loci\t${too_many_loci}" >> $stats_out
    echo -e "${id}\t${bc}\ttoo_many_mismatches\t${too_many_mismatches}" >> $stats_out
    echo -e "${id}\t${bc}\ttoo_short\t${too_short}" >> $stats_out
    echo -e "${id}\t${bc}\tother\t${other}" >> $stats_out
    echo -e "${id}\t${bc}\tchimeric\t${chimeric}" >> $stats_out

done



###############################################################################
echo "Calculation read length distributions"

# raw_data 
echo "   raw_data"
zcat ${raw_data}${id}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${stats_data}${id}_raw_readLenDist.tmp
awk '{printf("%s\t%s\n", $2, $1)}' ${stats_data}${id}_raw_readLenDist.tmp > ${stats_data}${id}_raw_readLenDist.dat
rm ${stats_data}${id}_raw_readLenDist.tmp

# trimmed_data
echo "   trimmed_data"
zcat ${trimmed_data}${id}_trimmed.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${stats_data}${id}_trimmed_readLenDist.tmp
awk '{printf("%s\t%s\n", $2, $1)}' ${stats_data}${id}_trimmed_readLenDist.tmp > ${stats_data}${id}_trimmed_readLenDist.dat
rm ${stats_data}${id}_trimmed_readLenDist.tmp

# filtered_data
echo "   filtered_data"
zcat ${filtered_data}${id}_final.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${stats_data}${id}_filtered_readLenDist.tmp
awk '{printf("%s\t%s\n", $2, $1)}' ${stats_data}${id}_filtered_readLenDist.tmp > ${stats_data}${id}_filtered_readLenDist.dat
rm ${stats_data}${id}_filtered_readLenDist.tmp

# BAM - split per barcode
bioawk=/home/vricci/miniconda3/bin/bioawk
echo "   projected to mouse_Tx"

flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")

    BAM=${work_dir}RSEM/${id}/${bc}/mouse_genome.Aligned.out.bam
    BAM_Tx=${work_dir}RSEM/${id}/${bc}/mouse_genome.Aligned.toTranscriptome.out.bam

    if [ ! -f ${stats_data}${id}_${bc}_mappedTx_readLenDist.dat ]; then
        samtools view $BAM_Tx | $bioawk -c sam '{hist[length($seq)]++} END {for (l in hist) print l, hist[l]}' | \
            sort -n -k1 > ${stats_data}${id}_${bc}_mappedTx_readLenDist.dat
    fi
done




###############################################################################
echo "Calculating gene and transcripts map types"

# mapType_Tx.py for BAM with only Tx and not Gene|Tx
GENEBED='/data/databases/mouse/Mmusculus.GRCm39.111.startStop.IDS.bed'
flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")
    map_out=${stats_data}${id}_${bc}_mapType.dat

    BAM=${work_dir}RSEM/${id}/${bc}/mouse_genome.Aligned.out.bam
    BAM_Tx=${work_dir}RSEM/${id}/${bc}/mouse_genome.Aligned.toTranscriptome.out.bam

    samtools view $BAM_Tx | ${work_dir}mapType_Tx.py -b $GENEBED -c 100000 > $map_out
done



###############################################################################
echo "Duplication rate"

dedup_out=${stats_data}${id}_dedup.dat
flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

touch $dedup_out
rm $dedup_out
cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")
    dedup_log=${work_dir}RSEM/${id}/${bc}/mouse_genome.Aligned.toTranscriptome.dedup.out

    # nb of input reads
    input_reads=$(awk '$0 ~ "Reads: Input Reads:"{print $NF}' $dedup_log)

    # nb of deduplicated reads
    output_reads=$(awk '$0 ~ "Number of reads out:"{print $NF}' $dedup_log)

    # duplication rate
    umi_duplication=$(awk -v var1=$input_reads -v var2=$output_reads 'BEGIN { print  ( var1 / var2 ) }')

    echo -e "${id}\t${bc}\tbefore_dedup\t${input_reads}" >> $dedup_out
    echo -e "${id}\t${bc}\tafter_dedup\t${output_reads}" >> $dedup_out
    echo -e "${id}\t${bc}\tduplication_rate\t${umi_duplication}" >> $dedup_out

done


###############################################################################
echo 'RSEM total counts'
RSEM_out=${stats_data}${id}_RSEM.dat
flagstat_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.flagstats
barcodes_split=${mapping_data}${id}/mouse_tRNA.Unmapped.out.mate1.split.barcodes

touch $RSEM_out
rm $RSEM_out
cat $flagstat_split | while read line; do 

    bc=$(echo $line | grep -zoP "(?s)(?<=out[.])(.*)(?=[.]mate1)")
    RSEM_genes=${RSEM_data}${id}/${bc}/${id}_${bc}.genes.results
    RSEM_isoforms=${RSEM_data}${id}/${bc}/${id}_${bc}.isoforms.results

    total_expected_count_genes=$(cut -f5 $RSEM_genes | tail -n +2 | awk '{s+=$1}END{print s}')
    total_expected_count_isoforms=$(cut -f5 $RSEM_isoforms | tail -n +2 | awk '{s+=$1}END{print s}')

    echo -e "${id}\t${bc}\texpected_count\t${total_expected_count_genes}" >> $RSEM_out
    # echo -e "${id}\t${bc}\texpected_count\t${total_expected_count_isoforms}" >> $RSEM_out
    # both are the same value!

done
