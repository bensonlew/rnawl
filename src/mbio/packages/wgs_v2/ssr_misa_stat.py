# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190320

"""
对通过misa设计的SSR中的文件进行统计，得到SSR类型统计
输入文件表头：ID	SSR nr.	SSR type	SSR	size	start	end
"""

import os
import argparse


def ssr_types_stat(infile, outfile, sample_name):
    """
    对通过misa设计的SSR中的文件进行统计，得到SSR类型统计
    输入文件表头：ID	SSR nr.	SSR type	SSR	size	start	end
    """
    types_dict = {"c": 0, "c*": 0, "p1": 0, "p2": 0, "p3": 0, "p4": 0, "p5": 0, "p6": 0}
    total_num = 0
    with open(infile, "r") as f:
        for line in f:
            if line.startswith("ID"):
                continue
            item = line.strip().split("\t")
            try:
                types_dict[item[2]] += 1
            except:
                raise Exception("{}不在SSR类型里".format(item[2]))
            total_num += 1
    with open(outfile, "w") as w:
        w.write("Sample ID\tSSR Number\tc\tc*\tp1\tp2\tp3\tp4\tp5\tp6\n")
        w.write(sample_name + "\t" + str(total_num) + "\t" + str(types_dict["c"]) + "\t" + str(types_dict["c*"]) + "\t")
        w.write(str(types_dict["p1"]) + "\t" + str(types_dict["p2"]) + "\t" + str(types_dict["p3"]) + "\t")
        w.write(str(types_dict["p4"]) + "\t" + str(types_dict["p5"]) + "\t" + str(types_dict["p6"]) + "\n")


parser = argparse.ArgumentParser(description="SSR类型统计")
parser.add_argument("-i", "--infile", required=True)
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("-s", "--sample_name", required=True)

args = vars(parser.parse_args())

ssr_types_stat(args["infile"], args["outfile"], args["sample_name"])
