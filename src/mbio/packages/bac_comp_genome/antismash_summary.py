# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
import argparse
from mbio.packages.bac_comp_genome.common_anno import CommonAnno


class antismash(object):
    """
    antismash的丰度统计和gene_list生成
    """
    def antismash_sum(self, dir, out):
        files = os.listdir(dir)
        with open(out, 'wb') as outfile:
            outfile.write("sample\tCluster ID\tType\n")
            for file in files:
                if re.search(r'.antismash_anno.xls', file):
                    sample = file.split(".antismash")[0]
                    with open (dir + "/" + file, "r") as f:
                        lines = f.readlines()
                        for line in lines[1:]:
                            line =line.strip().split("\t")
                            outfile.write(sample + "\t" + line[0] + "\t" + line[2] + '\n')

    def get_sample(self,file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
        return list

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[table_dir]', required=True, help='Input table file')
    parser.add_argument('-o', metavar='[output_dir]', required=True, help='output dir')
    parser.add_argument('-s', metavar='[sample_list]', required=True, help='sample list')
    args = parser.parse_args()
    table_dir = args.i
    out = args.o
    antismash = antismash()
    antismash.antismash_sum(table_dir, out + "/all.antismash_type.xls")
    sample_list = antismash.get_sample(args.s)
    anno = CommonAnno()
    anno.anno_abund(out + "/all.antismash_type.xls", "Type", "Cluster ID", out + "/all.antismash_abund.txt")
    anno.fill_data(out + "/all.antismash_abund.txt", sample_list,
                   out + "/all.antismash_abund.xls", "0", "Type")
    anno.anno_ano_genelist(out + "/all.antismash_type.xls", "Type", "Cluster ID", out + "/all.antismash_genelist.txt")
    anno.fill_data(out + "/all.antismash_genelist.txt", sample_list,
                   out + "/all.antismash_genelist.xls", "-", "Type")