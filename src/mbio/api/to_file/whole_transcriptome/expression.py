# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from biocluster.config import Config
import os
import json
from bson.objectid import ObjectId

project_type = 'whole_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_bam_list(data, option_name, dir_path, bind_obj):
    {'task_id': task_id, 'categories': ['mRNA', 'lncRNA', 'circRNA', 'smallRNA']}
    task_id = data['task_id']
    cursor = db['sg_exp'].find({'task_id': task_id})
    for document in cursor:
        if document['category'] in data['categories']:
            exp_id = document['main_id']
            for detail_document in db['sg_exp_detail'].find({'exp_id': exp_id}):
                pass


    output = os.path.join(dir_path, 'bam.list')
    open(output, 'w').writelines(sorted(['{}\n'.format(i['bam_path']) for i in results]))
    return output