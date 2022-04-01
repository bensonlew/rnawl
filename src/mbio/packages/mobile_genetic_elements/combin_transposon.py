# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(transposon, isfile, out):
    gff1 = []
    is_dict = {}
    with open(isfile, "r") as f:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            is_dict[lin[0]] = "\t".join(lin)
            if int(lin[2]) < int(lin[3]):
                gff1.append([lin[1], lin[2], lin[3], lin[0]])
            elif int(lin[3]) < int(lin[2]):
                gff1.append([lin[1], lin[3], lin[2], lin[0]])
    with open(transposon, "r") as f, open(out, "w") as d:
        lines = f.readlines()
        list1 = []
        for line in lines:
            lin = line.strip().split("\t")
            d.write("transposon\t{}\n".format("\t".join(lin)))
            for i in gff1:
                for j in [i[1], i[2]]:
                    if i[0] == lin[1]:
                        if int(lin[2]) <= int(j) and int(lin[3]) >= int(j):
                            list1.append(i[3])
                            print("aaaa"+ str(i))
                            continue
        list1 = set(list1)
        for key in is_dict.keys():
            if key not in list1:
                d.write("is\t{}\n".format(is_dict[key]))

def main():
    parser = OptionParser()
    parser.add_option('--p', dest ='transposon',metavar='[transposon result file]')
    parser.add_option('--i', dest='isfile', metavar='[is result file]')
    parser.add_option('--o', dest='output', metavar='[output file]')
    (options,args) = parser.parse_args()
    if not options.transposon or not options.output or not options.isfile:
        print "python combin_transposon.py --p transposon --i is --o output "
        return
    chang_gff(options.transposon, options.isfile, options.output)


if __name__=='__main__':
    main()