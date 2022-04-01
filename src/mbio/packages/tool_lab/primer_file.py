# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20200429
# 本脚本用于整理引物设计的传入文件为所需格式。
import argparse
import os
parse = argparse.ArgumentParser(description="用于引物设计文件整理")
parse.add_argument("-s", "--seq", help="传入文件序列", required=True)
# parse.add_argument("-t", "--target", help="从何处开始，往后多长的片段是用于扩增的片段", required=True)
parse.add_argument("-mt", "--tm_max", help="输出最大TM值", required=True)
parse.add_argument("-nt", "--tm_min", help="输出最小TM值", required=True)
parse.add_argument("-mg", "--gc_max", help="输出最大GC含量", required=True)
parse.add_argument("-ng", "--gc_min", help="输出最小GC含量", required=True)
parse.add_argument("-n", "--num", help="输出设计引物个数", required=True)
parse.add_argument("-l", "--size", help="输出引物的长度", required=True)
parse.add_argument("-p", "--product", help="输出产物长度范围", required=True)
parse.add_argument("-o", "--output", help="输出文件位置", required=True)
args = vars(parse.parse_args())
result_path = args["output"]
print result_path
with open(os.path.join(result_path, "prime_file.txt"), "w") as w:
    #w.write("PRIMER_SEQUENCE_ID=sca1_459_Indel_A_381.89_9\n")
    w.write("SEQUENCE_TEMPLATE=" + args["seq"] + "\n")
    # w.write("SEQUENCE_TARGET=" + args["target"] + "\n")
    w.write("PRIMER_PRODUCT_SIZE_RANGE=201-300\n")
    w.write("PRIMER_MIN_TM=" + args["tm_min"] + "\n")
    w.write("PRIMER_MAX_TM=" + args["tm_max"] + "\n")
    w.write("PRIMER_MAX_GC=" + args["gc_max"] + "\n")
    w.write("PRIMER_MIN_GC=" + args["gc_min"] + "\n")
    w.write("PRIMER_NUM_RETURN=" + args["num"] + "\n")
    w.write("PRIMER_OPT_SIZE=" + args["size"] + "\n")
    w.write("PRIMER_PRODUCT_SIZE_RANGE=" + args["product"] + "\n")
    w.write("PRIMER_TASK=generic\n")
    # w.write("PRIMER_MAX_END_STABILITY=250\n")
    w.write("=\n")




