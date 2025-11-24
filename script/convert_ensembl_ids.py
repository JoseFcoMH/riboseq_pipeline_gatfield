#!/usr/bin/env python3

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

