# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import argparse
import os
import gzip
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig


parser = argparse.ArgumentParser(description='Get reactome annotation from mongo by reading input  tabular file')
parser.add_argument('-i', metavar='[xml_table]', required=True, help='input BLAST tabular file')
parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='Output reactome annotation tabular file')
parser.add_argument('-v', metavar='[version]', default="72", help='reactome database version')

args = parser.parse_args()


class ReactomeAnnotation(object):
    '''
    Obtain detailed cog annotation information of metagenomic
    '''

    def __init__(self):
        self.client = Config().get_mongo_client(mtype='metagenomic', ref=True)

        self.gene2path = AnnotConfig().get_file_path(
            file ="Ensembl2Reactome_PE_All_Levels.txt",
            db = "reactome",
            version = "72")
        self.path_name = AnnotConfig().get_file_path(
            file ="ReactomePathways.txt",
            db = "reactome",
            version = "72")
        self.path_relation = AnnotConfig().get_file_path(
            file ="ReactomePathwaysRelation.txt",
            db = "reactome",
            version = "72")

    def get_reactome_class(self):
        """
        从参考库EGGNOG_sequence中找到gi号对应的description
        """
        gene2path_dict = dict()
        gene2name_dict = dict()
        path_dict = dict()
        path_relation_dict = dict()

        with open(self.path_name, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                path_dict[cols[0]] = cols[1]
                
        with open(self.gene2path, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                gene2name_dict[cols[1]] = cols[2].split(" ")[0]
                if cols[3] in path_dict:
                    if cols[1] in gene2path_dict:
                        gene2path_dict[cols[1]].add(cols[3])
                    else:
                        gene2path_dict[cols[1]] = set([cols[3]])



        return gene2name_dict, gene2path_dict, path_dict

    def main(self, align_table, anno_table, version="73"):
        print 'INFO: start building dictionary'
        gene2name_dict, gene2path_dict, path_dict = self.get_reactome_class()
        print gene2name_dict.items()[:5]
        print gene2path_dict.items()[:5]
        print path_dict.items()[:5]
        print 'INFO: start processing {}'.format(align_table)
        with open(align_table) as infile, open(anno_table, 'wb') as outfile:
            # infile.readline()
            outfile.write('#Seq_id\treactome_id\treactome_name\tpathway_id\tpathway_description\n')
            for line in infile:
                line = line.strip().split('\t')
                seq_id = line[0]
                reactomes = line[1]
                reactomes_all = list()
                reactomenames_all = list()
                paths_all = list()
                pathnames_all = list()
                pathlinks_all = list()
                for reactome in reactomes.split(";"):
                    if not reactome in gene2path_dict:
                        continue
                    paths = list(gene2path_dict[reactome])
                    pathnames = [path_dict[p] for p in paths]
                    pathlinks = ["{}&SEL={}".format(p, reactome) for p in paths]
                    reactomes_all.append(reactome)
                    reactomenames_all.append(gene2name_dict[reactome])
                    paths_all.extend(paths)
                    pathnames_all.extend(pathnames)
                    pathlinks_all.extend(pathlinks)

                if len(reactomenames_all) > 0:
                    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(seq_id,
                                                                    ";".join(reactomes_all),
                                                                    ";".join(reactomenames_all),
                                                                    ";".join(paths_all),
                                                                    ";".join(pathnames_all),
                                                                    ";".join(pathlinks_all)))


if __name__ == '__main__':
    if args.i and args.o:
        instance = ReactomeAnnotation()
        instance.main(args.i, args.o)
    elif args.h:
        parser.print_help()
