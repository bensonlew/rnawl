# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
import shutil


def chang_result(db, input, gff, out):
    anno = {}
    with open(db, "r") as f:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            anno[lin[0]] = lin[1]

    gff1 = {}
    with open(gff, "r") as f:
        lines = f.readlines()
        for line in lines:
            if re.search("^#", line):
                continue
            else:
                lin = line.strip().split("\t")
                gff1[lin[1]] = "\t".join([lin[0], lin[3], lin[4], lin[6]])

    with open(input, "r") as f, open(out, "w") as d:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                if lin[0] in gff1.keys() and lin[1] in anno.keys():
                    d.write(lin[0]+"\t"+gff1[lin[0]]+"\t"+anno[lin[1]]+"\n")

def main():
    parser = OptionParser()
    parser.add_option('--in', dest='blast', metavar='[balst file]')
    parser.add_option('--db', dest='db', metavar='[db description]')
    parser.add_option('--gff', dest='gff', metavar='[gff file]')
    parser.add_option("--o", dest="output", metavar="[out file]")
    (options,args) = parser.parse_args()
    if not options.blast or not options.db or not options.output:
        print "python enzyme_type.py --in blast --db db --gff gff --o output"
        return
    chang_result(options.db, options.blast, options.gff, options.output)

if __name__=='__main__':
    main()