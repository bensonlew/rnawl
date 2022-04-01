# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
import shutil


def chang_result(db, input, gff, out, type):
    anno = {}
    with open(db, "r") as f:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            anno[lin[0]] = lin[1]
    gff2 = {}
    with open(input, "r") as f:
        lines = f.readlines()
        for line in lines:
            if re.search("^#", line):
                continue
            else:
                lin = line.strip().split("\t")
                gff2[lin[0]] = lin[1]

    with open(gff, "r") as f, open(out, "w") as d:
        d.write("Gene_ID\tLocation\tStart\tEnd\tStrand\tType\n")
        lines = f.readlines()
        for line in lines[1:]:
            lin = line.strip().split("\t")
            gene = lin[1]
            if gene in gff2:
                new_gene = gff2[gene]
                if anno[new_gene] in [type] and type != 'all':
                    d.write(gene + "\t" +"\t".join([lin[0], lin[3], lin[4], lin[6]])+"\t"+anno[new_gene]+"\n")
                elif type == 'all':
                    d.write(gene + "\t" +"\t".join([lin[0], lin[3], lin[4], lin[6]])+"\t"+anno[new_gene]+"\n")
            else:
                d.write(gene + "\t" +"\t".join([lin[0], lin[3], lin[4], lin[6]])+"\t"+"protein"+"\n")

def main():
    parser = OptionParser()
    parser.add_option('--in', dest='blast', metavar='[balst file]')
    parser.add_option('--db', dest='db', metavar='[db description]')
    parser.add_option('--gff', dest='gff', metavar='[gff file]')
    parser.add_option("--o", dest="output", metavar="[out file]")
    parser.add_option("--type", dest="type", metavar="[enzyme type]")
    (options,args) = parser.parse_args()
    if not options.blast or not options.db or not options.output:
        print "python enzyme_type.py --in blast --db db --gff gff --o output --type enzyme"
        return
    chang_result(options.db, options.blast, options.gff, options.output, options.type)

if __name__=='__main__':
    main()