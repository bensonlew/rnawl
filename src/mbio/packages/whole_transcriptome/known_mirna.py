# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging

import pandas as pd
from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]

mature_mirna_col = Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)]['mature_mirna']


def main(args):
    count_matrix = args.i
    detail_table = args.o

    mature_ids = pd.read_table(count_matrix)['miRNA'].tolist()


def get_detail_info(mature_id):
    mature_id

    mature_mirna_col.find_one({'_id': mature_id})


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate known miRNA detail table')
    parser.add_argument('-i', action='store', required=True, help='known miRNA count matrix', metavar='<FILE>')
    parser.add_argument('-o', action='store', required=True, help='known miRNA detail table', metavar='<FILE>')

    args = parser.parse_args()

    main(args)
