# -*- coding: utf-8 -*-

import json
import os
import sys

from biocluster.config import Config
from biocluster.file import download

task_id = sys.argv[1]

database = Config().get_mongo_client(mtype='denovo_rna_v2')[Config().get_mongo_dbname('denovo_rna_v2')]

work_dir = os.getcwd()

group_lines = ['#sample\tgroup\n']
group_doc = database['sg_specimen_group'].find_one({'task_id': task_id})
for group, samples in zip(group_doc['category_names'], group_doc['specimen_names']):
    for sample in samples:
        group_lines.append('{}\t{}\n'.format(sample, group))
open(os.path.join(work_dir, 'group.txt'), 'w').writelines(group_lines)

control_lines = ['#control\tother\n']
control_doc = database['sg_specimen_group_compare'].find_one({'task_id': task_id})
for compare_name in json.loads(control_doc['compare_names']):
    control_lines.append('{}\t{}\n'.format(*compare_name.split('|')))
open(os.path.join(work_dir, 'control.txt'), 'w').writelines(control_lines)

count_file_dict = dict()
for document in database['sg_exp'].find({'task_id': task_id}):
    exp_level = document['exp_level']
    from_file = document['count_file']
    count_matrix = os.path.basename(from_file)
    to_file = os.path.join(work_dir, count_matrix)
    download(from_file, to_file)
    for together in (True, False):
        obj = {'count_matrix': count_matrix, 'group_table': 'group.txt', 'control_table': 'control.txt',
               'together': together}
        json.dump(obj, open(os.path.join(work_dir, '{}.{}.json'.format(
            exp_level, 'together' if together else 'split')), 'w'), indent=4)
