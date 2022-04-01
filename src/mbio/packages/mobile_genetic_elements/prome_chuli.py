# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(promoter, pre):
    with open(promoter, "r") as f, open(pre + ".bed", "w") as d, open(pre + ".xls", "w") as g:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                continue
            elif line.startswith(("ID")):
                lin = line.strip().split("\t")
                seqid = lin[1]
            elif line.startswith((">")):
                lin = line.strip().split("\t")
                li = lin[0].split(">")[1].split("..")
                d.write("{}\n".format("\t".join([seqid, li[0], li[1]])))
                g.write("{}\n".format("\t".join([seqid, li[0], li[1], lin[1], lin[2]])))

def main():
    parser = OptionParser()
    parser.add_option('--g', dest ='promoter',metavar='[promoter result file]')
    parser.add_option("--o", dest="prefix", metavar="[prefix of files]")
    (options,args) = parser.parse_args()
    if not options.promoter or not options.prefix:
        print "python prome_chuli.py --g promoter --o prefix"
        return
    chang_gff(options.promoter, options.prefix)


if __name__=='__main__':
    main()