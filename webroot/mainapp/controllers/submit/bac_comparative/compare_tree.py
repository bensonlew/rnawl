# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON

class CompareTreeAction(BacComparativeController):
    def __init__(self):
        super(CompareTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'sample_list', 'submit_location', 'task_type', 'gene_tree_id', 'species_tree_id', 'gene_type', "species_type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.compare_tree'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'gene_tree_id': data.gene_tree_id,
            'species_tree_id': data.species_tree_id,
            'sample_list': data.sample_list
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'TreeComparative' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bac_comparative.insert_main_table('comp_tree', mongo_data)
        update_info[str(main_id)] = 'comp_tree'
        gene_tree = self.bac_comparative.get_panfile(data.task_id, data.gene_tree_id)
        species_tree = self.bac_comparative.get_panfile(data.task_id, data.species_tree_id)
        options = {
            'update_info': json.dumps(update_info),
            'main_id': str(main_id),
            "gene_tree": gene_tree,
            'species_tree': species_tree,
            'sample_list': data.sample_list,
            'gene_type' : data.gene_type,
            'species_type': data.species_type
        }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="CompareTree/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(CompareTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)