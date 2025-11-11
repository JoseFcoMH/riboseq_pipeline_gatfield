#!/usr/bin/env python3

from __future__ import print_function
import sys
import re
import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio import Entrez
## import mysql.connector

# parse the command line arguments:
def parseArguments():
    parser = argparse.ArgumentParser(description='Given a SAM file of'
                                     ' aligned reads, for each read length '
                                     'it counts the type of transcript '
                                     'it maps to')
    parser.add_argument('-s','--sam',
                        help='SAM file',
                        type=str)
    parser.add_argument('-b','--bed',
                        help='BED file with regions where to count',
                        type=str,
                        required = True)
    parser.add_argument('-c','--count',
                        help='Max number of reads to count',
                        type=int)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return(parser.parse_args())

def FileCheck(fn):
    try:
        open(fn, "r")
        return 1
    except IOError:
        print("Error: file ", fn, " does not appear to exist.", file=sys.stderr)
        sys.exit()
        return 0

def readBed(fname, debug):
    handle = open(fname)
    geneType = dict()
    trnsType = dict()
    if debug: print('##   Reading file')
    for record in handle:
        line = record.rstrip()
        fields = re.split('[ \t]', line)
        chrom = fields[0]
        types = fields[3]
        ids = re.split(r'\|', chrom)
        tid = re.split(r'\|', types)
        geneType[ids[0]] = tid[0]
        trnsType[ids[1]] = tid[1]
    return(geneType, trnsType)

######################################################################
### Main starts here
######################################################################

if __name__ == "__main__":
    args = parseArguments()
    debug = args.debug
    if debug: print('## Debug mode ON')

    if args.sam:
        handle = open(args.sam)
        if debug: print('## SAM file:', args.sam)
    else:
        handle = sys.stdin
        if debug: print('## Get lines from STDIN')
        
    if args.bed:
        bedFile = args.bed
        if debug: print('## BED file:', args.sam)
        geneTypes, trnsTypes = readBed(bedFile, debug)
        
    if args.count:
        maxCount = args.count
        if debug: print('## Counts:', args.count)
    else:
        maxCount = 100000    

    if debug: print("## Reading SAM alignement...")

    # gCounts = dict()
    tCounts = dict()
    readID = list()
    readCount = 0
    
    for record in handle:
        line = record.rstrip()
        if re.match('^@', line):
            if debug: print('## Skipping line', line)
            continue
        fields = re.split('[ \t]', line)
        chrom = fields[2]
        guideID = fields[0]        
        guideFields = re.split(r'\|', guideID)
        matchStrand = fields[1]

        if guideID in readID:
            if debug: print('## Guide already done', guideID)
            continue
        else:
            readID.append(guideID)
        
        if matchStrand == 4:
            if debug: print('Guide', guideID, 'unmapped')
            continue
            
        readCount += 1
        if debug: print('## Counting:', readCount)
        if readCount > maxCount:
            if debug: print('## Reached max counts, stopping')
            break

        # ids = re.split('\|', chrom)
        # gid = ids[0]
        # tid = ids[1]
        tid = chrom
        # if gid in geneTypes:
        #     gType = geneTypes[gid]
        # else:
        #     gType = 'Other'
        if '*' not in tid:
            if tid in trnsTypes:
                tType = trnsTypes[tid]
            else:
                tType = 'Other'
                print(fields)
        # if gType in gCounts:
        #     gCounts[gType] += 1
        # else:
        #     gCounts[gType] = 1
            if tType in tCounts:
                tCounts[tType] += 1
            else:
                tCounts[tType] = 1

        if debug:
            print('Read ID :', guideID)
            # print('  Gene ID      :', gid)
            # print('  Gene type    :', gType)
            # print('  Type Counts  :', gCounts[gType])
            print('  Trns ID      :', tid)
            print('  Trns type    :', tType)
            print('  Type Counts  :', tCounts[tType])


    # for types in gCounts:
    #     print('Gene', types, gCounts[types], sep="\t")

    for types in tCounts:
        print('Transcript', types, tCounts[types], sep="\t")

    
