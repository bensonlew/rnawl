# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import sys

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__ == '__main__':
    task_id = sys.argv[1]
    logger.debug('task id: {}'.format(task_id))
    db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]

    total_num_dict = {'ref': 0, 'new': 0}
    exp_id = db['sg_exp'].find_one({'task_id': task_id})['main_id']
    for doc in db['sg_exp_detail'].find({'exp_id': exp_id}):
        if doc['is_new']:
            total_num_dict['new'] += 1
        else:
            total_num_dict['ref'] += 1
    total_num_dict['all'] = total_num_dict['ref'] + total_num_dict['new']

    for main_doc in db['sg_annotation_stat'].find({'task_id': task_id}):
        stat_id = main_doc['main_id']
        logger.debug('stat id: {}'.format(stat_id))
        cursor = db['sg_annotation_stat_detail'].find({'stat_id': stat_id})
        data = {'ref': {t: set() for t in ('GO', 'Pfam', 'KEGG', 'Swiss-Prot', 'COG', 'NR', 'Total_anno')},
                'new': {t: set() for t in ('GO', 'Pfam', 'KEGG', 'Swiss-Prot', 'COG', 'NR', 'Total_anno')},
                'all': {t: set() for t in ('GO', 'Pfam', 'KEGG', 'Swiss-Prot', 'COG', 'NR', 'Total_anno')}}

        type_tuple = ('GO', 'Pfam', 'KEGG', 'Swiss-Prot', 'COG', 'NR')
        for doc in cursor:
            if doc['type'] in type_tuple:
                transcript_ids = doc['transcript_list'].split(',')
                data[doc['seq_type']][doc['type']].update(transcript_ids)

        with open('tmp.{}.txt'.format(stat_id), 'w') as fw:
            for seq_type in data:
                data[seq_type]['Total_anno'].update(*(data[seq_type][t] for t in type_tuple))
                total_anno_num = len(data[seq_type]['Total_anno'])
                fw.write('{}\t{}\t{}\n'.format(
                    seq_type, total_anno_num, float(total_anno_num) / total_num_dict[seq_type]))
