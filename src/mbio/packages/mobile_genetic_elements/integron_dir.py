# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
from collections import defaultdict
import shutil

def change_faa(faa, genome, dict, list,dir):
    for i in list:
        if os.path.exists(dir +"/tmp_" + i):
            shutil.rmtree(dir +"/tmp_" + i)
        os.mkdir(dir +"/tmp_" + i)
        seq_list = []
        outfa = dir +"/tmp_" + i + "/" + i + ".prt"
        out2 = dir +"/tmp_" + i + "/" + i + ".fst"
        for seq_record in SeqIO.parse(genome, "fasta"):
            if seq_record.id == i:
                SeqIO.write(seq_record, out2, "fasta")
        for seq_record in SeqIO.parse(faa, "fasta"):
            if seq_record.id in dict[i]:
                seq_list.append(seq_record)
        SeqIO.write(seq_list, outfa, "fasta")

def chang_gff(gff):
    dict = defaultdict(list)
    with open(gff, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(("#")):
                continue
            else:
                lin =line.strip().split("\t")
                temp =lin[8].split(";")
                id =re.search('ID=(.*)', temp[0])
                dict[lin[0]].append(id.group(1))
    return dict

def get_seqid(genome):
    seqids= []
    for seq_record in SeqIO.parse(genome, "fasta"):
        seqids.append(seq_record.id)
    return seqids

def main():
    parser = OptionParser()
    parser.add_option('--s', dest='name', metavar='[genome or metagenome file name]')
    parser.add_option('--f', dest ='faa',metavar='[gene protein file]')
    parser.add_option('--g', dest ='gff',metavar='[gene gff file]')
    parser.add_option("--o", dest="output", metavar="[out dir]")
    (options,args) = parser.parse_args()
    if not options.faa or not options.name or not options.gff or not options.output:
        print "python Is_gene_dir.py --s name --f faa --g gff --o dir"
        return
    id_dict = chang_gff(options.gff)
    list_id = get_seqid(options.name)
    change_faa(options.faa, options.name, id_dict, list_id, options.output)

if __name__=='__main__':
    main()