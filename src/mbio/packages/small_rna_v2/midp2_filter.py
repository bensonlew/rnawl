# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys

def filter_blast(blast_in, blast_out, num):
    with open(blast_in, 'r') as f, open(blast_out, 'w') as fo:
        align_lines = list()
        reads_now = ""
        for line in f:
            cols = line.strip().split("\t")
            if cols[0] != reads_now:
                if len(align_lines) > 0 and len(align_lines) <= num:
                    fo.write("\n".join(align_lines) + "\n")
                reads_now = cols[0]
                align_lines = list()
            else:
                pass

            [start, end] = cols[2].split("..")
            if int(float(cols[7])) == 1 and int(start) - 1 == 0 and int(cols[1]) - int(end) <=0:
                align_lines.append(line.strip())

        if len(align_lines) > 0 and len(align_lines) <= num:
            fo.write("\n".join(align_lines) + "\n")

if __name__ == "__main__":
    filter_blast(sys.argv[1], sys.argv[2], int(sys.argv[3]))
