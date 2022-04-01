# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import sys

from biocluster.config import Config
from bson.son import SON

from mbio.packages.ref_rna_v2.trans_step import step_count

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def fix_transcripts_step(task_id, step2stat):
    logging.info('start searching mongo by task id ({})'.format(task_id))
    db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
    transcripts_id = db['sg_transcripts'].find_one({'task_id': task_id})['main_id']
    data_list = list()
    for step, stat in step2stat.items():
        step_list = []
        fr = open(stat)
        next(fr)
        for line in fr:
            step_dic = dict()
            step_range = line.strip().split('\t')[0]
            num = line.strip().split('\t')[1]
            step_dic[step_range] = num
            step_list.append(step_dic)
        else:
            fr.close()
        data = [
            ('transcripts_id', transcripts_id),
            ('step', int(step)),
            ('step_data', step_list)
        ]
        data = SON(data)
        data_list.append(data)
    else:
        collection = db['sg_transcripts_step']
        collection.delete_many({'transcripts_id': transcripts_id})
        collection.insert_many(data_list)
        logging.info('succeed in inserting documents by id ({})'.format(transcripts_id))


def count_fasta(fasta_file, task_id):
    group_num = 10
    step2stat = dict()
    for step in (200, 300, 600, 1000):
        fasta_to_txt = os.path.join(os.getcwd(), 'tmp.{}.txt'.format(step))
        stat_out = os.path.join(os.getcwd(), 'stat.{}.txt'.format(step))
        step_count(fasta_file, fasta_to_txt, group_num, step, stat_out)
        step2stat.update({step: stat_out})
        logging.info('succeed in exporting {}'.format(stat_out))
    else:
        fix_transcripts_step(task_id, step2stat)


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    task_id = sys.argv[2]
    count_fasta(fasta_file, task_id)
