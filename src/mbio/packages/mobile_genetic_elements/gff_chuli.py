# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(gff, out):
    with open(gff, "r") as f, open(out, "w") as d:
        lines = f.readlines()
        n = 1
        for line in lines:
            if line.startswith(("#")):
                continue
            else:
                lin = line.strip().split("\t")
                if lin[2] in ["CDS", "cds", "Cds"]:
                    if lin[6] in ["+"]:
                        d.write("{}\n".format("\t".join([lin[0], lin[3], lin[4]])))
                    elif lin[6] in ["-"]:
                        d.write("{}\n".format("\t".join([lin[0], lin[4], lin[3]])))

def main():
    parser = OptionParser()
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.gff or not options.output:
        print "python gff_chuli.py --g gff --o out"
        return
    chang_gff(options.gff, options.output)


if __name__=='__main__':
    main()