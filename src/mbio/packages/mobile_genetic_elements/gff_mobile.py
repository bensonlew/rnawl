# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(gff, mobile, out, start, end, seqid):
    gff1 = []
    with open(gff, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                continue
            else:
                lin = line.strip().split("\t")
                if int(lin[3]) < int(lin[4]):
                    gff1.append([lin[0], lin[3], lin[4]])
                elif int(lin[4]) < int(lin[3]):
                    gff1.append([lin[0], lin[4], lin[3]])
    with open(mobile, "r") as f, open(out, "w") as d:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            list1 = []
            for j in [lin[int(start)], lin[int(end)]]:
                for i in gff1:
                    if lin[int(seqid)] == i[0]:
                        if int(i[1]) <= int(j) and int(i[2]) >= int(j):
                            list1.append("false")
                            print("aaaa"+ str(i))
                            continue
                        else:
                            list1.append("true")

            list1 = set(list1)
            print(list1)
            if "false" in list1:
                pass
            else:
                d.write("{}\n".format("\t".join(lin)))

def main():
    parser = OptionParser()
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option('--m', dest='mobile', metavar='[mobile file]')
    parser.add_option('--o', dest='output', metavar='[output file]')
    parser.add_option("--s", dest="start", metavar="[column of start]")
    parser.add_option('--e', dest='end', metavar='[column of end]')
    parser.add_option("--i", dest="seqid", metavar="[column of seqid]")
    (options,args) = parser.parse_args()
    if not options.gff or not options.output or not options.mobile or not options.start or not options.end or not options.seqid:
        print "python gff_chuli.py --g gff --m mobile --o out --s 1 --e 2 --i 0"
        return
    chang_gff(options.gff, options.mobile, options.output, options.start, options.end, options.seqid)


if __name__=='__main__':
    main()