# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re, os
import shutil


def chang_gff(isfile, out):
    with open(isfile, "r") as f, open(out, "w") as g:
        lines = f.readlines()
        is_dict = {}
        n =0
        for line in lines:
            if line.startswith(("#")):
                continue
            else:
                lin = line.strip().split("\t")
                if lin[2] == "insertion_sequence":
                    m =re.search("ID=(.*);family=(.*);cluster=(.*)$",lin[8])
                    is_id =m.group(1)
                    is_dict[is_id] = {'location':lin[0],"start":lin[3],"end":lin[4],"strand":lin[6],"family":m.group(2),"cluster":m.group(3)}
                elif lin[2] == "terminal_inverted_repeat" and n == 0:
                    m = re.search("ID=(.*);parent=(.*)$", lin[8])
                    is_id = m.group(2)
                    postion = lin[3]+".."+lin[4]
                    if is_id in is_dict.keys():
                        is_dict[is_id]["IR1"] = postion
                    n +=1
                elif lin[2] == "terminal_inverted_repeat" and n == 1:
                    m = re.search("ID=(.*);parent=(.*)$", lin[8])
                    is_id = m.group(2)
                    postion = lin[3] + ".." + lin[4]
                    if is_id in is_dict.keys():
                        is_dict[is_id]["IR2"] = postion
                    n = 0
        for key,value in is_dict.items():
            print(is_dict[key].keys())
            if "IR1" not in is_dict[key].keys():
                g.write("{}\n".format("\t".join([key, is_dict[key]['location'], is_dict[key]['start'], is_dict[key]['end'],is_dict[key]['strand'], is_dict[key]['family'],is_dict[key]['cluster'],"-",'-'])))
            else:
                g.write("{}\n".format("\t".join([key, is_dict[key]['location'], is_dict[key]['start'], is_dict[key]['end'],is_dict[key]['strand'], is_dict[key]['family'],is_dict[key]['cluster'],is_dict[key]['IR1'],is_dict[key]['IR2']])))

def main():
    parser = OptionParser()
    parser.add_option('--g', dest ='isfile',metavar='[is result file]')
    parser.add_option("--o", dest="out", metavar="[out file]")
    (options,args) = parser.parse_args()
    if not options.isfile or not options.out:
        print "python is_chuli.py --g is --o out"
        return
    chang_gff(options.isfile, options.out)


if __name__=='__main__':
    main()