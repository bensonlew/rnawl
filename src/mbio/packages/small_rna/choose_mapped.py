# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import math
import collections
from Bio import SeqIO


def choose_fa(uniq, fa_list):
    seq_dict = dict()
    for seq in SeqIO.parse(uniq, "fasta"):
        seq_id = seq.id
        seq_id = seq_id.split("_")[1]
        seq = seq.seq
        seq_dict[seq_id] = seq
    with open(fa_list, 'r') as f, open(fa_list + ".fa", 'w') as fo:
        for line in f:
            s = line.strip()
            s_id = s.split("_")[1]
            fo.write(">{}\n{}\n".format(s, seq_dict[s_id]))


if __name__ == "__main__":
    u = sys.argv[1]
    l = sys.argv[2]
    choose_fa(u, l)
