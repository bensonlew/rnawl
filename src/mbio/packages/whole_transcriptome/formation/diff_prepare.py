# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os
import shutil
import sys

import pandas as pd
from biocluster.config import Config
from biocluster.file import download

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]


def main(args):
    json_dict = json.load(open(args.json))
    logger.info('succeed in parsing arguments as follow ->\n{}'.format(json_dict))

    database = Config().get_mongo_client(mtype='whole_transcriptome',db_version =args.version)[Config().get_mongo_dbname('whole_transcriptome',db_version =args.version)]

    task_id = json_dict['task_id']
    level = json_dict['level']
    category = json_dict['category']
    kind = json_dict['kind']
    background = json_dict['background']
    work_dir = json_dict['work_dir']
    count_matrix = json_dict['count_matrix']
    kind_table = json_dict['kind_table']

    remote_dir = database['task'].find_one({'task_id': task_id})['output']
    if level == 'G':
        from_file = os.path.join(remote_dir, 'other/count/G.reads.txt')
    elif level == 'T':
        if category in ('mRNA', 'lncRNA'):
            from_file = os.path.join(remote_dir, 'other/count/T.reads.txt')
        elif category == 'miRNA':
            from_file = os.path.join(remote_dir, 'other/count/S.reads.txt')
        elif category == 'circRNA':
            from_file = os.path.join(remote_dir, 'other/count/C.reads.txt')
    to_file = os.path.join(work_dir, os.path.basename(from_file))
    download(from_file, to_file)
    exp_id = database['exp'].find_one({'task_id': task_id, 'level': level})['main_id']
    seq_id = 'transcript_id' if level == 'T' else 'gene_id'
    if category in ('mRNA', 'lncRNA'):
        filter_dict = {'exp_id': exp_id}
        if background == 'mRNA,lncRNA':
            filter_dict['category'] = {'$in': ['mRNA', 'lncRNA']}
        elif background in ('mRNA', 'lncRNA'):
            filter_dict['category'] = category
        if kind != 'all':
            filter_dict['kind'] = kind
        cursor = database['exp_detail'].find(filter_dict)
        seq_ids = [document[seq_id] for document in cursor]
        count_df = pd.read_table(to_file, index_col=0)
        count_df = count_df.query('index in @seq_ids')
        count_df.to_csv(count_matrix, sep='\t')
    else:
        shutil.copy(to_file, count_matrix)

    cursor = database['exp_detail'].find({'exp_id': exp_id})
    data = [(document[seq_id], document['kind']) for document in cursor]
    kind_df = pd.DataFrame(data, columns=['seq_id', 'kind'])
    kind_df.to_csv(kind_table, sep='\t', index=False)

    for name in ('count_matrix', 'kind_table'):
        if os.path.isfile(json_dict[name]):
            logger.info('succeed in exporting {}'.format(json_dict[name]))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Build expression matrix')
    parser.add_argument('-j', action='store', required=True,
                        help='setting file in JSON format', dest='json')
    parser.add_argument('-v', action='store', required=True, type=int, default=1,
                        help='mongo_db_version', dest='version')
    # parser.add_argument('json', action='store', help='setting file in JSON format')

    args = parser.parse_args()

    if hasattr(args, 'json'):
        main(args)
    else:
        parser.print_help()
        sys.exit(-1)
