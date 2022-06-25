# -*- coding: utf-8 -*-
# __author__ = 'yuanshaohua'

import argparse
import os
import gzip

from requests import head
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig


parser = argparse.ArgumentParser(description='Get cog annotation from mongo by reading input BLAST tabular file')
parser.add_argument('-i', metavar='[xml_table]', required=True, help='input BLAST tabular file')
parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='Output cog annotation tabular file')
args = parser.parse_args()



func2cate_dict = {
    'A': 'INFORMATION STORAGE AND PROCESSING',
    'B': 'INFORMATION STORAGE AND PROCESSING',
    'C': 'METABOLISM',
    'D': 'CELLULAR PROCESSES AND SIGNALING',
    'E': 'METABOLISM',
    'F': 'METABOLISM',
    'G': 'METABOLISM',
    'H': 'METABOLISM',
    'I': 'METABOLISM',
    'J': 'INFORMATION STORAGE AND PROCESSING',
    'K': 'INFORMATION STORAGE AND PROCESSING',
    'L': 'INFORMATION STORAGE AND PROCESSING',
    'M': 'CELLULAR PROCESSES AND SIGNALING',
    'N': 'CELLULAR PROCESSES AND SIGNALING',
    'O': 'CELLULAR PROCESSES AND SIGNALING',
    'P': 'METABOLISM',
    'Q': 'METABOLISM',
    'R': 'POORLY CHARACTERIZED',
    'S': 'POORLY CHARACTERIZED',
    'T': 'CELLULAR PROCESSES AND SIGNALING',
    'U': 'CELLULAR PROCESSES AND SIGNALING',
    'V': 'CELLULAR PROCESSES AND SIGNALING',
    'W': 'CELLULAR PROCESSES AND SIGNALING',
    'Y': 'CELLULAR PROCESSES AND SIGNALING',
    'Z': 'CELLULAR PROCESSES AND SIGNALING'
}

func2desc_dict = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton'
}


class MetagenomicCogAnnotation(object):
    '''
    Obtain detailed cog annotation information of metagenomic
    '''

    def __init__(self):
        # self.client = Config().get_mongo_client(mtype='metagenomic', ref=True)
        # self.mongodb = self.client[Config().get_mongo_dbname('metagenomic', ref=True)]
        # self.eggNOG_ID = self.mongodb.eggNOG4_seqID
        # self.eggNOG = self.mongodb.eggNOG4
        self.seq_nog = dict()

        self.eggnog_tax = AnnotConfig().get_file_path(
            file ="40674_members.tsv.gz",
            db = "eggnog",
            version = "202006")
        self.eggnog_class = AnnotConfig().get_file_path(
            file ="40674_annotations.tsv.gz",
            db = "eggnog",
            version = "202006")

    def get_eggnog_tax_db(self):
        """
        从参考库EGGNOG_sequence中找到gi号对应的description
        """
        self.gene2eggnog = dict()
        with gzip.open(self.eggnog_tax, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                genes = cols[4].split(",")
                for gene in genes:
                    gene_id = ".".join(gene.split(".")[1:])
                    if gene_id in self.gene2eggnog:
                        self.gene2eggnog[gene_id].add(cols[1])
                    else:
                        self.gene2eggnog[gene_id] = set([cols[1]])


    def get_eggnog_class(self):
        """
        从参考库EGGNOG_sequence中找到gi号对应的description
        """
        nog2detail_dict = dict()
        nog2class_dict = dict()
        with gzip.open(self.eggnog_class, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                try:
                    nog2detail_dict[cols[1]] = cols[3]
                    nog2class_dict[cols[1]] = cols[2]
                except:
                    nog2detail_dict[cols[1]] = ""
                    nog2class_dict[cols[1]] = cols[2]

        return nog2detail_dict, nog2class_dict

    def main(self, align_table, anno_table):
        print 'INFO: start building dictionary'
        self.get_eggnog_tax_db()
        nog2detail_dict, nog2class_dict = self.get_eggnog_class()
        # nog_seq2nog_dict = {document['nog_seq']: document['nog'] for document in self.eggNOG_ID.find({})}
        # nog2detail2_dict = {document['nog']: document for document in self.eggNOG.find({})}
        print 'INFO: start processing {}'.format(align_table)
        with open(align_table) as infile, open(anno_table, 'wb') as outfile:
            header = infile.readline()
            outfile.write('#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len\n')
            NOGseq_list = dict()
            NOG_list = dict()
            for line in infile:
                line = line.strip().split('\t')
                if header.startswith("Score"):
                    seq_id = line[5]
                    gene = line[10].split(".")[1]
                    nogs = ";".join(list(self.gene2eggnog.get(gene, [])))
                    # print nogs
                else:
                    seq_id = line[0]
                    nogs = line[1].strip(";")
                for nog in nogs.split(";"):
                    if nog == '':
                        continue
                    nog_detail = nog2detail_dict.get(nog, "")
                    NOG = nog
                    Function = ";".join(nog2class_dict[nog])

                    Function_des = "; ".join([func2desc_dict[x] for x in Function.split(";")])
                    Category = "; ".join([func2cate_dict[x] for x in Function.split(";")])
                    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(seq_id, NOG, nog_detail, Function, Function_des, Category, 100, 0))


if __name__ == '__main__':
    if args.i and args.o:
        instance = MetagenomicCogAnnotation()
        instance.main(args.i, args.o)
    elif args.h:
        parser.print_help()
