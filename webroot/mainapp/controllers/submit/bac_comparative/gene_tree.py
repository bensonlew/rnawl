# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON


class GeneTreeAction(BacComparativeController):
    def __init__(self):
        super(GeneTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'analysis_type', 'submit_location', 'task_type', 'bootstrap', 'merit', "sample_list"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.gene_tree'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if data.analysis_type not in ["hous_keeping", "ortholog"]:
            info = {'success': False, 'info': '%s参数不正确，不是house_keeping,ortholog其中之一!' % data.analysis_type}
            return json.dumps(info)
        if data.analysis_type in ["ortholog"]:
            if not hasattr(data, 'gene_path'):
                info = {'success': False, 'info': '%s参数缺少!' % data.gene_path}
                return json.dumps(info)
            if not hasattr(data, 'homolog_id'):
                info = {'success': False, 'info': '%s参数缺少!' % data.homolog_id}
                return json.dumps(info)
        elif data.analysis_type in ["hous_keeping"]:
            if not hasattr(data, 'core_gene'):
                info = {'success': False, 'info': '%s参数缺少!' % data.core_gene}
                return json.dumps(info)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'analysis_type': data.analysis_type,
            'bootstrap': data.bootstrap,
            'merit': data.merit,
            'method':data.method,
            'sample_list': data.sample_list,
        }
        if data.analysis_type in ["ortholog"]:
            params['type'] = data.type
            params['homolog_id'] = data.homolog_id
            params['clusters'] = data.clusters
            if hasattr(data, 'out_group'):
                params['out_group'] = data.out_group
        elif data.analysis_type in ["hous_keeping"]:
            params['core_names'] = data.core_names
            if hasattr(data, 'out_group'):
                params['out_group'] = data.out_group
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = ''
        if data.analysis_type in ["ortholog"]:
            name = "Ortholog"
        elif data.analysis_type in ["hous_keeping"]:
            name = "House_keeping"
        main_table_name = 'GeneTree_' + name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.bac_comparative.insert_main_table('tree', mongo_data)
        update_info[str(main_id)] = 'tree'
        options = {
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "bootstrap": int(data.bootstrap),
                   'merit': data.merit,
                   'sample_list': data.sample_list
                   }
        if data.analysis_type in ["ortholog"]:
            options['gene_path'] = data.gene_path
            cluster_file = self.bac_comparative.get_panfile(data.task_id, ObjectId(data.homolog_id))
            options['homolog'] = cluster_file
            options['clusters'] = data.clusters
            options['type'] = data.type
            options['analysis_type'] = 'ortholog'
            if hasattr(data, 'out_group'):
                options['out_group'] = data.out_group
        elif data.analysis_type in ["hous_keeping"]:
            options['analysis_type'] = 'house_keeping'
            options['core_names'] = data.core_names
            options['core_gene_blast'] = data.core_gene + "/all.coregene.xls"
            if hasattr(data, 'out_group'):
                options['out_group'] = data.out_group
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="GeneTree/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(GeneTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)