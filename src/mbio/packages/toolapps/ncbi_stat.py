# -*- coding: utf-8 -*-
# __author__ = 'hao.gao' @20210115

import argparse
from biocluster.config import Config
import os,re


class AnnoTaxon(object):
    """
    ncbi下载数据统计
    """
    def run_taxon_anno(self, sample_list, s16, house, out):
        """:
        """
        house_k = {}
        for file in os.listdir(house):
            if re.search('.matches.m8', file):
                name = file.split(".matches.m8")[0]
                num = self.get_num(house + "/"+file)
                house_k[name] = num
        s16_k = {}
        with open(s16, "r") as g:
            lines = g.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                s16_k[lin[0]] = lin[3]
        with open(out, "w") as j:
            for sample in sample_list.split(";"):
                if sample in s16_k.keys():
                    de = s16_k[sample]
                else:
                    de = "-"
                if sample in house_k.keys():
                    des = house_k[sample]
                else:
                    des = "-"
                j.write("{}\t{}\t{}\t{}\n".format(sample, de, des, "Reference"))

    def get_num(self,file):
        with open(file, "r") as f:
            lines = f.readlines()
        return len(lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[sample_list]',required=True,help='sample_list')
    parser.add_argument('-s', metavar='[16s file]', required=True, help='16s file')
    parser.add_argument('-d', metavar='[house_keeping dir]', required=True, help='house_keeping dir')
    parser.add_argument('-o', metavar='[output file]',required=True,help='output file name')
    args = parser.parse_args()
    taxon_anno = AnnoTaxon()
    taxon_anno.run_taxon_anno(args.i, args.s, args.d, args.o)