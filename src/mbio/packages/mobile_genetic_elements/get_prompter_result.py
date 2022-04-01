# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(prompter, result, out):
    dict_r ={}
    with open(prompter, "r") as f:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            ss = lin[0]+"\t"+lin[1]+"\t"+lin[2]+"\t+\t"+lin[3]+"\t"+lin[4]
            dict_r[lin[0]+"\t"+lin[1]+"\t"+lin[2]] = ss

    with open(result, "r") as f, open(out, "w") as d:
        lines = f.readlines()
        d.write("Location\tStart\tEnd\tStrand\tLength\tseq\n")
        for line in lines:
            lin = line.strip().split("\t")
            if int(lin[3]) == 0:
                des = lin[0] + "\t" + lin[1] + "\t" + lin[2]
                if des in dict_r.keys():
                    d.write("{}\n".format(dict_r[des]))

def main():
    parser = OptionParser()
    parser.add_option('--p', dest='prompter',metavar='[prompter file]')
    parser.add_option('--s', dest='result',metavar='[prompter result file]')
    parser.add_option("--o", dest="output", metavar="[out file]")
    (options,args) = parser.parse_args()
    if not options.prompter or not options.result or not options.output:
        print "python get_prompter_result.py --p prompter --s result --o out"
        return
    chang_gff(options.prompter, options.result, options.output)


if __name__=='__main__':
    main()