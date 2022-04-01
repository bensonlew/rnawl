# -*- coding: utf-8 -*-
# gao.hao

import pandas as pd
import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='[blast]', required=True, help='blast table')
parser.add_argument("-d", metavar="[db file]", required=False, help="plasmid db file")
parser.add_argument("-f", metavar="[fasta file]", required=False, help="fasta file")
parser.add_argument("-t", metavar="[table file]", required=False, help="table file")
parser.add_argument('-o', metavar='[output]', required=True, help='outfilename')
args = parser.parse_args()

inputfile = args.i
out = args.o
fasta = args.f
file2 = args.t
db = args.d

def get_gc(file, id):
    gc = 0
    num = 0
    for iterator in SeqIO.parse(file, "fasta"):
        if iterator.id ==id:
            C_num =iterator.seq.count("C")
            G_num = iterator.seq.count("G")
            g_num = iterator.seq.count("g")
            c_num = iterator.seq.count("c")
            gc = float(C_num+G_num+g_num+c_num) / len(iterator.seq)*100
            num = len(iterator.seq)
    return gc,num

def run():
    seq_id = []
    dict2 = {}
    list1 = []
    with open(file2, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            lin = line.strip().split("\t")
            dict2[lin[0]] =lin[-1]
            if lin[3] == "plasmid" or lin[4] == "known-Plasmid":
                seq_id.append(lin[0])
    dict = {}
    with open(inputfile,"r") as f:
        lines = f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            gc,num = get_gc(fasta, lin[0])
            de =lin[1]+"\t"+lin[0]
            if de not in dict:
               dict[de] = [dict2[lin[0]], lin[0],str(gc), str(num), str(lin[2]), str(lin[10])]
            list1.append(lin[0])
    print(seq_id)
    print(list1)
    seq_no = list(set(seq_id) - set(list1))
    print(seq_no)
    with open(db, "r") as g,open(out,"w") as k:
        k.write("location\tseq_id\tgc\tlength\tidentity\tevlue\tacc_num\ttaxon\tplasmid_type\tstrain\tplasmid_name\n")
        if len(seq_no) >0:
            for i in seq_no:
                gc, num = get_gc(fasta, i)
                k.write("{}\n".format(("\t").join([dict2[i],i,str(gc),str(num),"-","-","-","-","-","-","-"])))
        lines = g.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            if len(lin) == 49:
                des = "-"
            elif len(lin) > 49:
                des = lin[49]
            for ss in dict.keys():
                if lin[1] in ss.split("\t"):
                    taxon = ";".join([lin[38], lin[36], lin[34], lin[32], lin[30], lin[28], lin[26]])
                    if re.search("plasmid", lin[2]):
                        train = lin[2].split("plasmid")[0].strip(" ")
                        pla_name = lin[2].split("plasmid")[1].split(",")[0].strip(" ")
                        if pla_name == '':
                            pla_name = '-'
                    else:
                        train = "-"
                        pla_name = '-'
                    k.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\n".format("\t".join(dict[ss]), lin[1], taxon, des, train, pla_name))


if __name__ == '__main__':
    run()