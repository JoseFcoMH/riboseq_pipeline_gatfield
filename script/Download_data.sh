path_home='/work/FAC/FBM/CIG/dgatfiel/default/vricci/Lisa/DENR_woCrossLink_Ribo/'
raw_data=${path_home}raw_data/

mkdir -p $raw_data

cd $path_home


for RunID in 498 503 508 549 557; do
        while read -r line; do
                echo $line

                fq=$(echo $line | cut -d'/' -f6) # DENR_L2_R1_001_6HJUzkBTZEHG.fastq.gz
                LibID=$(echo $fq | cut -d'_' -f1) # DENR

                LaneID=$(echo $fq | cut -d'_' -f2 | cut -c2) # 2 of L2
                ReadID=$(echo $fq | cut -d'_' -f3) # R1
                SeqID=$(echo $fq | cut -d'_' -f4) # 001

                FQ=${LibID}_${RunID}_${LaneID}_${SeqID}_${ReadID}.fastq.gz
                echo $FQ

                if [ ! -f "${raw_data}${FQ}" ]; then
                        echo 'Downloading...'
                        wget $line -O ${raw_data}${FQ}
                fi
        done < LIMS_samples.links_${RunID}
done
