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

### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
path_home='/path/to/project/' # don't forget the '/' at the end

### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 
### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT ### TO EDIT 

raw_data=${path_home}raw_data/

mkdir -p $raw_data

cd $path_home

# using UHTS-LIMS of Lausanne Genomic Technologies Facility
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
        done < LIMS_samples.links_${RunID} # including the links for downloads
done
