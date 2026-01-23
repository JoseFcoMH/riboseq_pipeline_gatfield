#!/usr/bin/env python3

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


## This script is used to caclulate the translation start and stop
## positions in a transcript given a GFF file.


import HTSeq
import collections
import argparse

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='Given a GTF file the script '
                                     'outputs the transcripts ID, coordinates '
                                     'of CDS or UTRs, and gene and '
                                     'transcripts biotypes.')
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str,
                        required=True)
    parser.add_argument('-t',
                        help='Coordinates to output. Options are CDS,'
                        ' 5UTR and 3UTR',
                        nargs=1,
                        type=str,
                        required=True)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


args = parserFunction()

debug   = args.debug
gtfFile = args.gtf[0]
tout = args.t[0]

opt = ['5UTR', 'CDS', '3UTR']

if tout not in opt:
    print ('-t', tout, 'not recognised, available options: 5UTR, CDS,',
           '3UTR')
    exit()

## Read GTF file
gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
intervals = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

tPos = dict()
cdsLength = collections.Counter( )
exons = dict()
transcriptLength = collections.Counter( )
cdsStart = dict()
geneId = dict()
geneBiotype = dict()
trnsBiotype = dict()

if debug: print('Reading GTF file')

for feature in gtf_file:
    if len(str(feature.iv.chrom)) > 2:
        continue
    if feature.type == 'transcript':
        transcript_id = feature.attr['transcript_id']
        tPos[transcript_id] = feature.iv
        geneId[transcript_id] = feature.attr['gene_id']
        geneBiotype[transcript_id] = feature.attr['gene_biotype']
        trnsBiotype[transcript_id] = feature.attr['transcript_biotype']
        cdsLength[transcript_id] = 0
        if debug:
            print('Transcript id:', transcript_id)
            print('  GID:', geneId[transcript_id])
            print('  GBT:', geneBiotype[transcript_id])
            print('  TBT:', trnsBiotype[transcript_id])
            print('  CHR:', str(feature.iv.chrom))
    if feature.type == 'exon':
        transcript_id = feature.attr['transcript_id']
        exons[transcript_id] = feature.iv
        transcriptLength[transcript_id] += feature.iv.length
        if debug: print('  Cumulative Transcript length: ', transcriptLength[transcript_id])
    if feature.type == 'CDS':
        transcript_id = feature.attr['transcript_id']
        cdsLength[transcript_id] += feature.iv.length
        if debug:
            print('  CDS length:', feature.iv.length)
            print('  Cumulative CDS length:', cdsLength[transcript_id])
    if feature.type == 'start_codon':
        transcript_id = feature.attr['transcript_id']
        strand = feature.iv.strand
        if transcript_id in cdsStart:
            continue
        else:
            if strand == '+':
                start = feature.iv.start
                tStart = tPos[transcript_id].start
                # cdsStart[transcript_id] = start - tStart
                cdsStart[transcript_id] = transcriptLength[transcript_id] - (exons[transcript_id].end - start - 1)
            else:
                start = feature.iv.end
                tStart = tPos[transcript_id].end
                # cdsStart[transcript_id] = tStart - start
                cdsStart[transcript_id] = transcriptLength[transcript_id] - (start - exons[transcript_id].start - 1)
            if debug:
                print('  Strand: ', strand)
                print('  CDS start: ', start)
                print('  Transcript Start: ', tStart)
                print('  CDS Trns Start: ', cdsStart[transcript_id])



for transcript_id in sorted( geneId ):      
    if transcript_id in cdsStart:
        start = cdsStart[transcript_id] - 1
        end = start + cdsLength[transcript_id]
        if end > transcriptLength[transcript_id] and debug:
            print('  ERROR! CDS length > transcript length')
    else:
        start = 0
        end = 9999999999
    if tout == 'CDS':
        print(geneId[transcript_id] + '|' + transcript_id, start, end,
              geneBiotype[transcript_id] + '|' + trnsBiotype[transcript_id],
              sep='\t')
    elif tout == '5UTR' and start > 0:
        print(geneId[transcript_id] + '|' + transcript_id, 0, start-1,
              geneBiotype[transcript_id] + '|' + trnsBiotype[transcript_id],
              sep='\t')
    elif tout == '3UTR' and end != 9999999999:
        print(geneId[transcript_id] + '|' + transcript_id, end+1,
              transcriptLength[transcript_id], geneBiotype[transcript_id] +
              '|' + trnsBiotype[transcript_id], sep='\t')
