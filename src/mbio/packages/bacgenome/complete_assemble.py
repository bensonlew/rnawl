# -*- coding: utf-8 -*-
import pandas as pd
import argparse
from Bio import SeqIO
import re,os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-f', metavar='[fasta]', required=True, help='fasta file')
parser.add_argument('-d', metavar='[plasmid detail]', required=True, help='plasmid detail file')
parser.add_argument('-p', metavar='[prefix]', required=True, help='prefix of file')
parser.add_argument('-dir', metavar='[dir]', required=True, help='output dir')
args = parser.parse_args()

fasta = args.f
detail = args.d
prefix = args.p
dir = args.dir

def run():
    dict ={}
    if os.path.exists(dir +"/seq_dir"):
        shutil.rmtree(dir +"/seq_dir")
    os.mkdir(dir +"/seq_dir")
    with open(detail) as f:
        lines = f.readlines()
        for line in lines[1:]:
            lin = line.strip().split("\t")
            dict[lin[5]] = lin[4]
    with open(prefix+"_assembly_details.xls","w") as f,open(prefix+"_assembly_summary.xls","w") as k, open("plasmid.type.xls","w") as s:
        k.write("Chromosome No.\tPlasmid No.\tGenome Size (bp)\tG+C (%)\n")
        chrr = []
        plsmid = []
        total_len = 0
        gc_total = 0
        for iterator in SeqIO.parse(fasta, "fasta"):
            SeqIO.write(iterator, dir +"/seq_dir/"+iterator.id + ".fasta", "fasta")
            if re.search("Plasmid", iterator.id) or re.search("plasmid", iterator.id):
                plsmid.append(iterator)
                gene = iterator.id[0].lower()+iterator.id[-1]+"_gene"
                print(iterator.id,gene)
                s.write("{}\t{}\n".format(iterator.id, gene))
            elif re.search("Chromosome", iterator.id) or re.search("chromosome", iterator.id):
                chrr.append(iterator)
            a = iterator.seq.upper().count("A")
            t = iterator.seq.upper().count("T")
            c = iterator.seq.upper().count("C")
            g = iterator.seq.upper().count("G")
            hh = g + c
            gc_total +=hh
            total = len(iterator.seq)
            total_len +=total
            gc = float(g+c)/total*100
            f.write("{}\t{}\t{}\t{}\n".format(iterator.id, total, gc, "A:"+str(a)+";T:"+str(t)+";G:"+str(g)+";C:"+str(c)))
        print(len(chrr))
        if len(chrr) >= 1:
            SeqIO.write(chrr, "chromosome.fasta", "fasta")
        if len(plsmid) >= 1:
            SeqIO.write(plsmid, "plasmid.fasta", "fasta")
        gc = float(gc_total)/total_len*100
        dd =len(chrr)
        ff = len(plsmid)
        k.write("{}\t{}\t{}\t{}\n".format(dd, ff, total_len, gc))

if __name__ == '__main__':
    run()