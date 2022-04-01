#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 2021/5/8 13:52
# @Author  : U make me wanna surrender my soul
import os


def go_circ(file, num, output, rscript):
    output_pdf = os.path.join(output, 'output/go_circ.pdf')
    script_name = os.path.join(output, 'go_circ.r')
    f = open(script_name, 'w')
    f.write('library(GOplot)\n')
    f.write('circ <- read.table("{}",sep="\t",header=T)\n'.format(file))
    f.write('gghaha <- GOCircle(circ, nsub = {},label.size=3,table.legend = T,rad1=2,rad2=3)\n'.format(num))
    width = 5 + num
    heigh = 5 + num
    f.write('ggsave("{}",gghaha,width = {}, height = {})\n'.format(output_pdf, width, heigh))
    f.close()
    os.system('{} go_circ.r'.format(rscript))



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="请输入GO富集分析文件", type=str, required=True)
    parser.add_argument("-n", help="请输入需要展示的GO Term数量,默认为10", type=int, required=False, default=10)
    parser.add_argument("-o", help="请输入输出图片的名称", type=str, required=True)
    parser.add_argument("-r", help="R软件路径", type=str, required=True)
    args = parser.parse_args()
    go_circ(args.i, args.n, args.o, args.r)
