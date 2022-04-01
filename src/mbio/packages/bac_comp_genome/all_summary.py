# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
import argparse

class all_summary_anno(object):
    """
    注释信息和总览信息表的汇总
    """
    def sum_by_anno(self, anno_sum, prephage, island, antismash, out):
        pre_dict = self.get_dict(prephage, "_prephage_summary.xls", ";", 1, -1)
        is_dict = self.get_dict(island, ".GI_summary.xls", ",", 1, -1)
        ant_dict = self.get_dict(antismash, ".antismash_anno.xls", ",", 0, -4)
        print (ant_dict)
        with open (anno_sum, "r") as f,open (out, "w") as g:
            lines = f.readlines()
            names = lines[0].strip().split("\t")
            g.write("\t".join(names) + "\tPrephage ID\tIsland ID\tAntismash ID\n" )
            for line in lines[1:]:
                lin = line.strip().split("\t")
                des = lin[-1] + "___" + lin[0]
                if des in pre_dict.keys():
                    lin.append(pre_dict[des])
                else:
                    lin.append("-")
                if des in is_dict.keys():
                    lin.append(is_dict[des])
                else:
                    lin.append("-")
                if des in ant_dict.keys():
                    lin.append(ant_dict[des])
                else:
                    lin.append("-")
                g.write("\t".join(lin) + "\n")

    def get_dict(self, dir, des, type, num, num2):
        files = os.listdir(dir)
        dict = {}
        for file in files:
            if re.search(r"{}".format(des), file):
                name = file.split(des)[0]
                with open (dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        lin = line.strip().split("\t")
                        genes = lin[num2].split(type)
                        for gene in genes:
                            dess =name + "___" + gene
                            dict[dess] = lin[num]
        return dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[anno_sum]', required=True, help='Input anno_sum table')
    parser.add_argument('-pre', metavar='[prephage]', help='tcdb function', default="function_database")
    parser.add_argument('-isl', metavar='[island]', help='tcdb function', default="function_database")
    parser.add_argument('-ant', metavar='[antismash]', help='tcdb function', default="function_database")
    parser.add_argument('-o', metavar='[out table]', required=True, help='output file name')
    args = parser.parse_args()
    anno_sum = args.i
    prephage = args.pre
    island = args.isl
    antismash = args.ant
    anno_table = args.o
    all_sum = all_summary_anno()
    all_sum.sum_by_anno(anno_sum, prephage, island, antismash, anno_table)