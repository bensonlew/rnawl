# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(dir, out):
    n = 0
    with open(out, "w") as g:
        for i in os.listdir(dir):
            sample = i.split(".transposon.xls")[0]
            with open(dir + "/" + i, "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    n += 1
                    g.write("{}\n".format("\t".join(['transposon' + str(n), sample, str(lin[1]), str(lin[2]), str(lin[3])])))

def main():
    parser = OptionParser()
    parser.add_option('--d', dest='dir',metavar='[transposon result dir]')
    parser.add_option("--o", dest="out", metavar="[out files]")
    (options,args) = parser.parse_args()
    if not options.dir or not options.out:
        print "python transposon_stat.py --d dir --o out"
        return
    chang_gff(options.dir, options.out)

if __name__=='__main__':
    main()