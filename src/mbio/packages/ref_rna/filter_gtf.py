# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# 此脚本用于过滤gtf

import sys
import os
import shutil

class FilterGtf(object):
    def __init__(self):
        self.chr = dict()
        pass

    def get_chr_range(self, bedfile):
        '''
        根据bed文件获取染色体坐标长度
        '''
        with open(bedfile, 'r') as f:
            for line in f.readlines():
                cols = line.strip().split("\t")
                self.chr.update({cols[0]: (int(cols[1]), int(cols[2]))})
        
    def filter_gtf(self, gtf_in, gtf_out):
        '''
        根据bed坐标范围限制gtf边界
        '''
        with open(gtf_in, 'r') as f, open(gtf_out, 'w') as w:
            for line in f.readlines():
                cols = line.strip().split("\t")
                if len(cols) < 9:
                    pass
                else:
                    start = self.chr[cols[0]][0]
                    end = self.chr[cols[0]][1]
                    if int(cols[3]) > end or int(cols[4]) < start + 1:
                        continue
                    if int(cols[3]) < start + 1:
                        cols[3] = str(start + 1)
                    if int(cols[4]) > end:
                        cols[4] = str(end)
                w.write("\t".join(cols) + "\n")
    
    def filter_gtf_by_bed(self, gtf_in, gtf_out, bed):
        self.get_chr_range(bed)
        self.filter_gtf(gtf_in, gtf_out)


if __name__ == "__main__":
    FilterGtf().filter_gtf_by_bed(sys.argv[1], sys.argv[2], sys.argv[3])
