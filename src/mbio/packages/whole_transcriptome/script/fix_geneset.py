# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import sys

import pandas as pd

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

db = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]

if __name__ == '__main__':
    logger.debug('Parse sys argv: {}'.format(sys.argv))
    task_id = sys.argv[1]
    gene_type_fp = sys.argv[2]
    gene_diff_dir = sys.argv[3]
    mrna_diff_dir = sys.argv[4]

    category = 'mRNA'
    category_df = pd.read_table(gene_type_fp, names=['seq_id', 'category'], index_col=0, usecols=[0, 2])
    for fname in os.listdir(gene_diff_dir):
        if fname.endswith('.detail.txt'):
            df = pd.read_table(os.path.join(gene_diff_dir, fname), index_col='seq_id')
            df = df.join(category_df)
            df.to_csv(os.path.join(mrna_diff_dir, fname), sep='\t')

    map_dict = dict()
    diff_dir = mrna_diff_dir
    for fname in os.listdir(diff_dir):
        if fname.endswith('.detail.txt'):
            map_dict[fname[:-11]] = os.path.join(diff_dir, fname)

    logger.debug('Map dict: {}'.format(map_dict))

    cursor = db['geneset'].find({'task_id': task_id, 'level': 'G'})
    for main_doc in cursor:
        logger.debug('Find main document: {}'.format(main_doc))
        geneset_id = main_doc['main_id']
        vs_pair = main_doc['name'][7:]
        if vs_pair not in map_dict:
            continue
        logger.debug('Vs pair: {}'.format(vs_pair))
        df = pd.read_table(map_dict[vs_pair])
        df = df[df['significant'] == 'yes']
        category_list = df['category'].tolist()
        logger.debug('Category list length: {}'.format(len(category_list)))
        update_dict = {'seq_list': df['seq_id'].tolist(),
                       'category_list': df['category'].tolist(),
                       'kind_list': df['kind'].tolist(),
                       'regulate_list': df['regulate'].tolist()}
        logger.debug('Geneset id: {}'.format(geneset_id))
        db['geneset_detail'].update_one({'geneset_id': geneset_id}, {'$set': update_dict}, upsert=True)
