# -*- coding: utf-8 -*-
import os
import math
import subprocess
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
import glob
import matplotlib
import gzip
import json
import copy
import pandas as pd
from collections import defaultdict
from collections import OrderedDict
import xml.etree.ElementTree as ET

from mbio.api.database.gene_db import genome


class EnsembleKarotype(object):
    def __init__(self):
        self.chr_karyotype = dict()

    def get_species_karyotype(self, seq_region_file, karyotype_file, output):

        with gzip.open(seq_region_file, 'r') as f:
            seq_region_id2name = {l.split("\t")[0]:l.split("\t")[1] for l in f}

        with gzip.open(karyotype_file, 'r') as f, open(output + '.genome2karyotype.txt', 'w') as fo:
            for line in f:
                cols = line.strip("\n").split("\t")
                chr_name = seq_region_id2name[cols[1]]
                if chr_name not in self.chr_karyotype:
                    self.chr_karyotype[chr_name] = [(int(cols[2]), int(cols[3]), cols[4])]
                else:
                    self.chr_karyotype[chr_name].append((int(cols[2]), int(cols[3]), cols[4]))
                fo.write(line.strip("\n") + "\t" + chr_name + "\n")


    def get_species_gene_karyotype(self, gene_pos, output):
        self.gene_pos_dict = dict()
        self.gene_karyotype = OrderedDict()
        with open(gene_pos, 'r') as f:
            f.readline()
            for line in f:
                sg_gene_id = line.split("\t")[0]
                if len(sg_gene_id.split("_")) > 4:
                    chr_name =  "_".join(sg_gene_id.split("_")[:-3])
                    start = sg_gene_id.split("_")[-3]
                    end = sg_gene_id.split("_")[-2]
                    strand = sg_gene_id.split("_")[-1]
                else:
                    chr_name, start, end, strand = sg_gene_id.split("_")
                    
                if chr_name in self.gene_pos_dict:
                    self.gene_pos_dict[chr_name].append((
                        sg_gene_id, (int(start), int(end), strand)
                    ))
                else:
                    self.gene_pos_dict[chr_name] = [(
                        sg_gene_id, (int(start), int(end), strand)
                    )]

        def get_pos(pos, range_list):
            r_len = len(range_list)
            k_types = list()
            for k_start, k_end, k_type in range_list:
                if pos[0] >= k_start and pos[0] <= k_end:
                    k_types.append(k_type)
                if pos[1] >= k_start and pos[1] <= k_end:
                    if k_type not in k_types:
                        k_types.append(k_type)

            return "-".join(k_types)

        for chr_name in self.chr_karyotype:
            karyotype_list = sorted(self.chr_karyotype[chr_name], key=lambda x:x[0])
            gene_list = sorted(self.gene_pos_dict[chr_name],key=lambda x:x[1][0])

            for sg_gene_id, pos in gene_list:
                k_type = get_pos(pos, karyotype_list)
                self.gene_karyotype[sg_gene_id] = k_type

        with open(output + 'gene2karyotype.txt', 'w') as fo:
            for g,k in self.gene_karyotype.items():
                fo.write("{}\t{}\n".format(g,k))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-seq_region', type=str, required=True,
                        help='genome acc')
    parser.add_argument('-karyotype', type=str, required=True,
                        help='karyotype')
    parser.add_argument('-gene_pos', type=str, required=True,
                        help='gene_pos')
    parser.add_argument('-output', type=str, required=True,
                        help='output')
    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    ef = EnsembleKarotype()
    ef.get_species_karyotype(seq_region_file=args.seq_region, karyotype_file=args.karyotype, output=args.output)
    ef.get_species_gene_karyotype(gene_pos=args.gene_pos, output=args.output)
    # ef_j = json(ef.table)
    # print json.dumps(ef.table, indent=4)


