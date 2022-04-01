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

def merge_species_fasta(old, new, spe_class, spe):
    if os.path.exists(spe_class):
        pass
    else:
        os.system("mkdir -p {}".format(spe_class))

    print "{} {} ".format(old, new)

    if os.path.exists(old) and os.path.exists(new):
        seq_list = list()
        with open("{}/{}.faa".format(spe_class, spe), 'w') as fo:
            for seq in SeqIO.parse(new, "fasta"):
                if seq.id not in seq_list:
                    fo.write(">{}\n{}\n".format(seq.id, str(seq.seq)))
                seq_list.append(seq.id)
            for seq in SeqIO.parse(old, "fasta"):
                if seq.id not in seq_list:
                    fo.write(">{}\n{}\n".format(seq.id, str(seq.seq)))
                seq_list.append(seq.id)
    elif os.path.exists(old):
        os.link(old, "{}/{}.faa".format(spe_class, spe))
    elif os.path.exists(new):
        os.link(new, "{}/{}.faa".format(spe_class, spe))
    else:
        print "{} {} not exists ".format(spe_class, spe)



def species_para(species_file, fasta_dir):
    with open(species_file, 'r') as s_f:
        for line in s_f:
            cols = line.split("\t")
            spe_class = cols[-1].split(";")[1]
            old_file = "/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/genes/organisms/{}/{}.pep".format(cols[1], cols[0])
            new_file = os.path.join(fasta_dir, spe_class.rstrip("s").lower(), cols[1])


            merge_species_fasta(old_file, new_file, spe_class, cols[1])



if __name__ == "__main__":
    species = sys.argv[1]
    new_dir = sys.argv[2]
    species_para(species, new_dir)
