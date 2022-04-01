# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os
import sys

import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)




def main(args):
    json_dict = json.load(open(args.json))
    logger.info('succeed in parsing arguments as follow ->\n{}'.format(json_dict))

    database = Config().get_mongo_client(mtype='whole_transcriptome',db_version =args.version)[Config().get_mongo_dbname('whole_transcriptome',db_version =args.version)]

    group_dict = json.loads(json_dict['group_dict'])
    seq_id = 'transcript_id' if json_dict['level'] == 'T' else 'gene_id'

    columns = ['seq_id']
    group_lines = ['#sample\tgroup\n']
    if len(group_dict) == 1 and group_dict.keys()[0] == 'all':
        samples = group_dict['all']
        columns.extend(samples)
        group_lines.extend('{}\t{}\n'.format(sample, sample) for sample in samples)
    else:
        for group, samples in group_dict.items():
            columns.extend(samples)
            group_lines.extend('{}\t{}\n'.format(sample, group) for sample in samples)
    open(json_dict['group_table'], 'w').writelines(group_lines)

    main_dict = database['exp'].find_one({'task_id': json_dict['task_id'], 'level': json_dict['level']})
    if 'library' in json_dict:
        task_dict = database['task'].find_one({'task_id': json_dict['task_id']})
        rnas = [rna for rna, lib in task_dict['rna'].items() if lib == json_dict['library']]
        filter_dict = {'exp_id': main_dict['main_id'], 'category': {'$in': rnas}}
    else:
        filter_dict = {'exp_id': main_dict['main_id'], 'category': json_dict['category']}
    if json_dict['kind'] != 'all':
        filter_dict['kind'] = json_dict['kind']
    cursor = database['exp_detail'].find(filter_dict)
    df = pd.DataFrame(list(cursor))
    df = df.rename({seq_id: 'seq_id'}, axis=1)
    df = df.reindex(columns, axis=1)
    df.to_csv(json_dict['exp_matrix'], sep='\t', index=False)

    if json_dict['control_id']:
        document = database['specimen_group_compare'].find_one({'main_id': ObjectId(json_dict['control_id'])})
        compare_names = json.loads(document['compare_names'])
        control_lines = ['#control\tother\n']
        control_lines.extend('{}\t{}\n'.format(*compare_name.split('|')) for compare_name in compare_names)
        open(json_dict['control_table'], 'w').writelines(control_lines)

    for name in ('exp_matrix', 'group_table', 'control_table'):
        if os.path.isfile(json_dict[name]):
            logger.info('succeed in exporting {}'.format(json_dict[name]))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Build expression matrix')
    parser.add_argument('-j', action='store', required=True,
                        help='setting file in JSON format', dest='json')
    parser.add_argument('-v', action='store', required=True, type=int, default=1,
                        help='mongo_db_version', dest='version')

    args = parser.parse_args()

    if hasattr(args, 'json'):
        main(args)
    else:
        parser.print_help()
        sys.exit(-1)
