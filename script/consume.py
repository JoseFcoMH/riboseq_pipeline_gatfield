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

"""
A simple size filter for fastq formatted outputs based read size (MIN, MAX)
"""

import sys
import itertools


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2015-2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

VERSION = __version__


if len(sys.argv) > 1 and sys.argv[1] == '-v':
    sys.stdout.write("consume {}\n".format(VERSION))
    sys.exit(0)
try:
    min_seq_len = int(sys.argv[1])
except IndexError:
    min_seq_len = 21 #for RB, 26
try:
    max_seq_len = int(sys.argv[2])
except IndexError:
    max_seq_len = 60 #for RB, 35

total_count = 0
filtered_count = 0

# to account for newline character
eff_min_seq_len = min_seq_len + 1
eff_max_seq_len = max_seq_len + 1
if len(sys.argv) > 3:
    outfile = open(sys.argv[3], 'w')
else:
    outfile = sys.stdout

for identifier, seq, sep, qual in itertools.zip_longest(*[sys.stdin]*4):
    total_count += 1
    if len(seq) < eff_min_seq_len or len(seq) > eff_max_seq_len:
        filtered_count += 1
    else:
        outfile.write(identifier+seq+sep+qual)
outfile.close()
sys.stderr.write("size_filter.py version {}\n".format(VERSION))
sys.stderr.write("Total number of sequences: {}\n".format(total_count))
sys.stderr.write("Number of sequences filtered out: {}\n".format(filtered_count))
sys.stderr.write("Filter used: [min:{}, max:{}]\n".format(min_seq_len, max_seq_len))
