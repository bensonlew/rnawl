# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig


class AnnoSecAction(MetagenomicController):
    def __init__(self):
        super(AnnoSecAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['submit_location', 'task_type', 'nr_method', 'anno_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, "info": "PARAMETERS MISSING: %s!" , "variables": [argu], "code" : "C2404501"}
                return json.dumps(info)
        task_name = "metagenomic.report.anno_sec"
        module_type = "workflow"
        anno_info = self.metagenomic.get_anno_info(data.anno_id, "sec")
        task_info = self.metagenomic.get_task_info(anno_info['task_id'])
        nr_info = self.metagenomic.get_nr_info(data.nr_method, task_info["task_id"])  # metagenomic需要配置
        nr_table = nr_info['anno_file']
        geneset_id = nr_info["geneset_id"]
        reads_profile = self.metagenomic.get_geneset_info(geneset_id)['reads_num']
        params_json = {
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'nr_method': data.nr_method,
            'anno_id': data.anno_id
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_name = 'SEC_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('anno_sec', mongo_data)
        update_info = {str(main_table_id): 'anno_sec'}
        options = {
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'gene_sec_anno': anno_info['anno_file'],
            'gene_profile': reads_profile,
            'gene_nr_anno': nr_table
        }
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=[])
        task_info = super(AnnoSecAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

