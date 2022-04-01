# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON

class SpeciesTreeAction(BacComparativeController):
    def __init__(self):
        super(SpeciesTreeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'analysis_type', 'submit_location', 'task_type', 'bootstrap', 'sample_list']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.species_tree'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if data.analysis_type not in ["s16", "hous_keeping", "ortholog"]:
            info = {'success': False, 'info': '%s参数不正确，不是s16,hous_keeping,ortholog其中之一!' % data.analysis_type}
            return json.dumps(info)
        if data.analysis_type in ["s16"]:
            if not hasattr(data, 'gene_path'):
                info = {'success': False, 'info': '%s参数缺少!' % data.gene_path}
                return json.dumps(info)
            if not hasattr(data, 'method'):
                info = {'success': False, 'info': '%s参数缺少!' % data.method}
                return json.dumps(info)
            if data.method not in ["NJ","ML"]:
                info = {'success': False, 'info': '%s参数不正确!' % data.method}
                return json.dumps(info)
            else:
                if data.method in ["ML"]:
                    if not hasattr(data, 'merit'):
                        info = {'success': False, 'info': '%s参数缺少!' % data.merit}
                        return json.dumps(info)
            if not hasattr(data, 'sample_list'):
                info = {'success': False, 'info': '%s参数不正确!' % data.sample_list}
                return json.dumps(info)

        elif data.analysis_type in ["ortholog"]:
            if not hasattr(data, 'gene_path'):
                info = {'success': False, 'info': '%s参数缺少!' % data.gene_path}
                return json.dumps(info)
            if not hasattr(data, 'homolog_id'):
                info = {'success': False, 'info': '%s参数缺少!' % data.homolog_id}
                return json.dumps(info)
            if not hasattr(data, 'type'):
                info = {'success': False, 'info': '%s参数缺少!' % data.type}
                return json.dumps(info)
            if not hasattr(data, 'merit'):
                info = {'success': False, 'info': '%s参数缺少!' % data.merit}
                return json.dumps(info)
            if not hasattr(data, 'sample_list'):
                info = {'success': False, 'info': '%s参数不正确!' % data.sample_list}
                return json.dumps(info)
        elif data.analysis_type in ["hous_keeping"]:
            if not hasattr(data, 'core_gene'):
                info = {'success': False, 'info': '%s参数缺少!' % data.core_gene}
                return json.dumps(info)
            if not hasattr(data, 'merit'):
                info = {'success': False, 'info': '%s参数缺少!' % data.merit}
                return json.dumps(info)
            if not hasattr(data, 'sample_list'):
                info = {'success': False, 'info': '%s参数不正确!' % data.sample_list}
                return json.dumps(info)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'analysis_type': data.analysis_type,
            'bootstrap': data.bootstrap,
            'method': data.method,
            'sample_list': data.sample_list
        }
        if data.analysis_type in ["s16"]:
            if hasattr(data, 'out_group'):
                params['out_group'] = data.out_group
            elif data.method in ["ML"]:
                params['merit'] = data.merit
        elif data.analysis_type in ["ortholog"]:
            params['type'] = data.type
            params['homolog_id'] =data.homolog_id
            params['merit'] = data.merit
        elif data.analysis_type in ["hous_keeping"]:
            params['merit'] = data.merit
            if hasattr(data, 'out_group'):
                params['out_group'] = data.out_group
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = ''
        if data.analysis_type in ["s16"]:
            name = "16S"
        elif data.analysis_type in ["ortholog"]:
            name = "Ortholog"
        elif data.analysis_type in ["hous_keeping"]:
            name = "House_keeping"
        main_table_name = 'SpeciesTree_' + name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
                   'sample_list': data.sample_list
                   }
        if data.analysis_type in ["s16"]:
            options['analysis_type'] = "s16"
            options['method'] = data.method
            options['gene_path'] = data.gene_path
            if data.method in ["ML"]:
                options['merit'] = data.merit
            if hasattr(data, 'out_group'):
                options['out_group'] = data.out_group
        elif data.analysis_type in ["ortholog"]:
            options['analysis_type'] = "ortholog"
            options['gene_path'] = data.gene_path
            cluster_file = self.bac_comparative.get_panfile(data.task_id,ObjectId(data.homolog_id))
            options['homolog'] = cluster_file
            options['type'] = data.type
            options['merit'] = data.merit
        elif data.analysis_type in ["hous_keeping"]:
            options['analysis_type'] = "house_keeping"
            options['merit'] = data.merit
            options['core_gene'] = data.core_gene + "/all.cor_gene.fa"
            if hasattr(data, 'out_group'):
                options['out_group'] = data.out_group
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="SpeciesTree/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(SpeciesTreeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)