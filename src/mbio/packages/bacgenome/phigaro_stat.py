# -*- coding: utf-8 -*-
# gao.hao

import pandas as pd
import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='[result]', required=True, help='result table')
parser.add_argument("-d", metavar="[anno file]", required=False, help="anno file")
parser.add_argument("-f", metavar="[fasta file]", required=False, help="fasta file")
parser.add_argument('-p', metavar='[output]', required=True, help='outfilename')
args = parser.parse_args()

inputfile = args.i
out = args.p
fasta = args.f
db = args.d

def get_gc(file, id, ids, start, end):
    for iterator in SeqIO.parse(file, "fasta"):
        if iterator.id ==id:
            C_num =iterator.seq.count("C")
            G_num = iterator.seq.count("G")
            g_num = iterator.seq.count("g")
            c_num = iterator.seq.count("c")
            gc = float(C_num+G_num+g_num+c_num) / len(iterator.seq)*100
            seq = iterator.seq[int(start)-1:int(end)]
            num = len(seq)
            iterator.id = ids
            iterator.seq = seq
            return gc,num, iterator

def run():
    list2 = []
    total_len =0
    list1 =[]
    nums = 0
    with open(inputfile,"r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            nums +=1
            lin = line.strip().split("\t")
            phid = "Ph" + '{:d}'.format(nums).zfill(2)
            gc, num, li = get_gc(fasta, lin[0], phid, lin[1],lin[2])
            total_len +=num
            list1.append(li)
            list2.append([phid, lin[0], str(lin[1]), str(lin[2]), str(num), lin[4], str(gc)])
    SeqIO.write(list1, out +"_prephage.fna", "fasta")
    with open(out +".stat.xls","w") as f:
        f.write("{}\t{}\n{}\t{}".format("Ph No.", "Total Length", str(nums), str(total_len)))
    with open(db, "r") as g,open(out +"_prephage_detail.xls","w") as k, open(out +"_prephage_summary.xls","w") as s:
        s.write("Ph No.\tLocation\tPh Start\tPh End\tLength(bp)\tTaxon\tGC(%)\tCDS No.\n")
        k.write("Location\tPh No.\tGene Id\tStrand\tGene Start\tGene End\tGene Name\tNR Description\tCOG Type\n")
        lines = g.readlines()
        for i in list2:
            genes_num =0
            print(i)
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if lin[5] == i[1]:
                    print("aa")
                    if lin[2] >lin[3]:
                        start = lin[3]
                        end = lin[2]
                    else:
                        start = lin[2]
                        end = lin[3]
                    if int(start) >= int(i[2]) and int(end) <= int(i[3]):
                        print(lin[5])
                        genes_num +=1
                        k.write("{}\n".format("\t".join([lin[5], i[0], lin[0], lin[1], lin[2], lin[3], lin[14], lin[6], lin[10]])))
            i.append(str(genes_num))
            s.write("{}\n".format("\t".join(i)))


if __name__ == '__main__':
    run()