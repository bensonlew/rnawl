# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(transposon, out):
    with open(transposon, "r") as f, open(out, "w") as g:
        lines = f.readlines()
        n =0
        m =1
        start = 0
        end = 0
        for line in lines:
            if re.search("^A transposon", line):
                lin = re.split(r"[ ]+", line)
                start = lin[2].split("..")[0]
                n +=1
                print("aaa" + start)
            elif re.search("^B transposon", line):
                lin = re.split(r"[ ]+", line)
                end = lin[2].split("..")[1]
                n += 1
                print("bbbb" + end)
            if start > 0 and end > 0 and n % 2 == 0:
                print(start + end)
                length = int(end) - int(start)+1
                g.write("{}\n".format("\t".join(['transposon'+str(m), str(start), str(end), str(length)])))
                m +=1
                start =0
                end = 0


def main():
    parser = OptionParser()
    parser.add_option('--g', dest='transposon',metavar='[is result file]')
    parser.add_option("--o", dest="out", metavar="[out files]")
    (options,args) = parser.parse_args()
    if not options.transposon or not options.out:
        print "python transposon_chuli.py --g transposon --o out"
        return
    chang_gff(options.transposon, options.out)


if __name__=='__main__':
    main()