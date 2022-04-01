# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# kegg相关函数

import re
import os
import sys
import xml.etree.ElementTree as ET
import lxml.html
import chardet
import random
import requests
from pymongo import MongoClient
from multiprocessing import Pool
import time
from Bio import SeqIO

def get_species_fasta(fasta_file, ko_list, new_get_url):

    gene_ids = list()
    if os.path.exists(fasta_file):
        print "old fasta {}".format(fasta_file)
        for seq in SeqIO.parse(fasta_file, "fasta"):
            gene_ids.append(seq.id)
    with open(ko_list, 'r') as ko_f, open(new_get_url, 'w') as get_f:
        for line in ko_f:
            cols = line.split("\t")
            if cols[0] in gene_ids:
                pass
            else:
                get_f.write("https://www.kegg.jp/dbget-bin/www_bget?-f+-n+a+{}\n".format(cols[0]))


def species_para(species_file):
    with open(species_file, 'r') as s_f:
        for line in s_f:
            cols = line.split("\t")
            old_file = "/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/genes/organisms/{}/{}.pep".format(cols[1], cols[0])
            ko_list = "/mnt/ilustre/users/sanger-dev/app/database/Annotation/download2019/kegg_20190924/ko2gene/{}.ko2gene".format(cols[1])

            new_get_url = "/mnt/ilustre/users/sanger-dev/app/database/Annotation/download2019/kegg_20190924/getfasta/{}.get_fasta.txt".format(cols[1])
            if os.path.exists(ko_list):
                pass
            else:
                print ko_list
                continue
            if os.path.exists(new_get_url):
                pass
            else:
                get_species_fasta(old_file, ko_list, new_get_url)

if __name__ == "__main__":
    species = sys.argv[1]
    species_para(species)
