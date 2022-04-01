# -*- coding: utf-8 -*-

import sys
import logging
import os
import numpy as np
import biom
import argparse
import hashlib
from Bio import SeqIO
from biocluster.config import Config
# from q2_types.feature_data import DNAIterator


class FormatReads(object):
    """
    标准化序列和丰度表
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.biom = os.path.join(self.software_dir, "program/Python/bin/biom")

    def run_taxon_anno(self, biom_table, reference, out_file):
        """
        转换表格和序列
        :param biom_table:输入的丰度表；
        :param reference:参考注释信息
        :param out_dir:结果目录
        :return:
        """
        asv_table = os.path.join(out_file, "ASV_table.xls")
        asv_reads = os.path.join(out_file, "ASV_reps.fasta")
        number = 1
        asv_origin_list = []
        asv_dict = {}
        with open(biom_table, 'r') as f, open(asv_table, 'w') as w:
            for line in f:
                if line.startswith("# Constructed from biom file"):
                    pass
                elif line.startswith("#OTU ID"):
                    line = line.strip().split("\t")
                    line[0] = "ASV ID"
                    w.write("\t".join(line)+"\n")
                else:
                    line = line.strip().split("\t")
                    asv_list = [float(x) for x in line[1:]]
                    # asv_name = "ASV" + str(number)
                    asv_name = hashlib.md5(line[0].encode('utf-8')).hexdigest()
                    origin_name = line[0]
                    total_asv_size = int(float(np.sum(asv_list)))
                    if total_asv_size >= 2:
                        number += 1
                        if asv_name not in asv_origin_list:
                            asv_origin_list.append(asv_name)
                            asv_dict[origin_name] = asv_name
                        w.write(asv_name + "\t" + "\t".join(line[1:]) + "\n")

        with open(asv_reads, "w") as wf:
            for seq_record in SeqIO.parse(reference, 'fasta'):
                gene_id = seq_record.id
                gene_seq = str(seq_record.seq)
                if gene_id in asv_dict:
                    new_name = asv_dict[gene_id]
                    wf.write(">{}\n{}\n".format(new_name, gene_seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[asv table]',required=True,help='Input asv table')
    parser.add_argument('-r', metavar='[ref reads]', required=True, help='reference db reads')
    parser.add_argument('-o', metavar='[output dir]',required=True,help='output dir name')
    args = parser.parse_args()
    biom_table = args.i
    reference = args.r
    out_file = args.o
    taxon_anno = FormatReads()
    taxon_anno.run_taxon_anno(biom_table, reference, out_file)