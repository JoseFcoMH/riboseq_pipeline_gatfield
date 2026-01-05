#!/usr/bin/env python3

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

import sys
TRIM_OFF_VERSIONS = True
if len(sys.argv) < 2 or sys.argv[1] == '-':
    infile = sys.stdin
else:
    infile = open(sys.argv[1])
for line in infile:
    if line[0] == '>':
        parsed = line.strip().split()
        for k,v in [(p+':').split(':')[:2] for p in parsed]:
            if k == 'gene' and v != '':
                if TRIM_OFF_VERSIONS:
                    v_trim = v.split('.')[0]
                    i_trim = parsed[0][1:].split('.')[0]
                    sys.stdout.write(" ".join(['>'+v_trim+'|'+i_trim] +
                                              ['ori_id:'+parsed[0][1:]] +
                                              parsed[1:]) + '\n')
                else:
                    sys.stdout.write(" ".join(['>' + v + '|' + parsed[0][1:]] +
                                              parsed[1:]) + '\n')
    else:
        sys.stdout.write(line)

