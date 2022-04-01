# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

import os
import argparse
from pymongo import MongoClient
from biocluster.config import Config


def cog(gene_file, cog_anno, location):
    gene_type = {}
    cog_color = {}
    db_name = 'bacgenome'
    db = Config().get_mongo_client(mtype=db_name, ref=True)[Config().get_mongo_dbname(db_name, ref=True)]
    #db = up['sanger_biodb']
    for i in db.COG_color.find({"version":2}):
        #cog_color[i['_id']] = i['rgb']
        cog_color[i['fid']] = i['rgb']
    with open(cog_anno, 'r') as f2:
        lines = f2.readlines()
        for line in lines[1:]:
            line = line.strip().split('\t')
            gene_type[line[0]] = cog_color[line[3][0]]
    with open(gene_file, 'r') as f1, open('sense_strand_cog.txt', 'w')as o1, open('antisense_strand_cog.txt',
                                                                                  'w') as o2:
        lines = f1.readlines()
        for line in lines[1:]:
            line = line.strip().split('\t')
            if line[4] == "+":
                if line[0] in gene_type:
                    color = gene_type[line[0]]
                else:
                    color = '(169,169,169)'
                o1.write(' '.join([location, line[7], line[8], 'fill_color=' + color]) + '\n')
            if line[4] == "-":
                if line[0] in gene_type:
                    color = gene_type[line[0]]
                else:
                    color = '(169,169,169)'
                o2.write(' '.join([location, line[8], line[7], 'fill_color=' + color]) + '\n')

def _main():
    parser = argparse.ArgumentParser(description='add gene_info in your table ')
    parser.add_argument('-g', '--gene', help="gene.gff")
    parser.add_argument('-c', '--cog', help="anno cog file")
    parser.add_argument('-l', '--location', help="location")
    args = parser.parse_args()
    cog(args.gene, args.cog, args.location)


if __name__ == "__main__":
    _main()