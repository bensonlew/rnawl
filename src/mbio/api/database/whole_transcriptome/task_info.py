# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

import datetime
import json
import os

from biocluster.config import Config

from mbio.api.database.whole_transcriptome.api_base import ApiBase
from mbio.packages.rna.annot_config import AnnotConfig


class TaskInfo(ApiBase):
    def __init__(self, bind_object):
        super(TaskInfo, self).__init__(bind_object)

    def add_task_info(self, data_json, annot_group=None):
        main_dict = json.load(open(data_json))
        main_dict['task_id'] = main_dict['id']
        main_dict['created_ts'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        #版本升级v1.1 20201027
        # main_dict['version'] = 'v1'
        main_dict['version'] = 'v1.3'
        if annot_group:
            self.annot_config_dict = AnnotConfig().get_group_option_detail(section=annot_group)
            self.get_db_version(main_dict, annot_group)

        self.create_db_table('task', [main_dict])

    def get_db_version(self, main_dict, annot_group=None):
        annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        if self.annot_config_dict['kegg']['version'] > "2020":
            if False:
                annot_version_dict['kegg'] += "_spe"
        else:
            del annot_version_dict['kegg']
        main_dict.update({
            'database_version': annot_version_dict,
            'annot_group': annot_group
        })

    def add_lib_rna(self, task_id, lib_dict):
        task_doc = self.db['task'].find_one({'task_id': task_id})
        record_id = task_doc['_id']
        lib = [k for k, v in lib_dict.items() if v]
        long_task_id = task_doc['options']['long_task_id']
        long_task_doc = self.db['task'].find_one({'task_id': long_task_id})
        rna = {'mRNA': 'long'}
        if 'lncRNA' in long_task_doc['options']['rna_select']:
            rna['lncRNA'] = 'long'
        if 'circRNA' in long_task_doc['options']['rna_select']:
            rna['circRNA'] = 'long'
        if lib_dict['circle']:
            rna['circRNA'] = 'circle'
        if lib_dict['small']:
            rna['miRNA'] = 'small'
        genome_id = long_task_doc['options']['genome_id']
        lnc_ref = self.add_lnc_ref(genome_id, long_task_doc)
        insert_dict = {'lib': lib, 'rna': rna, 'genome_id': genome_id, 'lnc_ref': lnc_ref}
        self.update_db_record('task', record_id, insert_dict=insert_dict)

    def add_lnc_ref(self, genome_id, long_task_doc):
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        genome_doc = collection.find_one({'genome_id': genome_id})
        lnc_ref = True if 'lnc_dir' in genome_doc else False
        if 'lncRNA' not in long_task_doc['options']['rna_select']:
            lnc_ref = False
        return lnc_ref

    def add_sub_task(self, task_id, long_task_id, small_task_id, circle_task_id):
        task_doc = self.db['task'].find_one({'task_id': task_id})
        record_id = task_doc['_id']
        insert_dict = dict()
        if long_task_id:
            long_task_info = self.db['task'].find_one({'task_id': long_task_id})
            insert_dict['long_task'] = long_task_info
        if small_task_id:
            small_task_info = self.db['task'].find_one({'task_id': small_task_id})
            insert_dict['small_task'] = small_task_info
        if circle_task_id:
            circle_task_info = self.db['task'].find_one({'task_id': circle_task_id})
            insert_dict['circle_task'] = circle_task_info
        if 'sub_output' in task_doc:
            sub_output = task_doc['sub_output']
        else:
            sub_output = {'long': long_task_info['output']}
            if small_task_id:
                sub_output['small'] = small_task_info['output']
            if circle_task_id:
                sub_output['circle'] = circle_task_info['output']
        insert_dict['sub_output'] = sub_output
        self.update_db_record('task', record_id, insert_dict=insert_dict)

    def add_pdf_dir(self, task_id):
        task_doc = self.db['task'].find_one({'task_id': task_id})
        self.db['known_mirna'].update({'task_id': task_id}, {'$set': {
            'pdf_dir': os.path.join(task_doc['output'], 'mirna/03_sRNA_Analysis/01_Known_miRNA/known_pre_structure')
        }}, upsert=True)
        self.db['novel_mirna'].update({'task_id': task_id}, {'$set': {
            'pdf_dir': os.path.join(task_doc['output'], 'mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_pre_structure')
        }}, upsert=True)

    def add_diff_method(self, task_id, long_diff_method):
        task_doc = self.db['task'].find_one({'task_id': task_id})
        record_id = task_doc['_id']
        options = task_doc['options']
        options['long_diff_method'] = long_diff_method
        insert_dict = {'options': options}
        self.update_db_record('task', record_id, insert_dict=insert_dict)
