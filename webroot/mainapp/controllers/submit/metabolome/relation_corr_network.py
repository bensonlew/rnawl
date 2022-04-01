# -*- coding: utf-8 -*-

import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from biocluster.config import Config
import pandas as pd
import os, re
from biocluster.file import download
from bson import ObjectId

class RelationCorrNetworkAction(MetabolomeController):
    def __init__(self):
        super(RelationCorrNetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        global metab_table_path, metab_desc_path
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_table', 'group_id',
                        'metab_set_table', 'group_detail', 'trans_exp_main_id',
                        'trans_geneset_main_id', 'coefficient', 'padjust_method', 'sort', 'top']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_corr_network'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        name = "RelationCorrNetwork_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if not data.coefficient in ["spearman", "pearson", "kendall"]:
            info = {'success': False, 'info': '相关性算法类型错误!', 'code': 'C2300308'}
            return json.dumps(info)
        if not data.padjust_method in ["fdr_bh", "fdr_by", "holm", "bonferroni"]:
            info = {'success': False, 'info': '多重检验矫正方法错误!', 'code': 'C2300308'}
            return json.dumps(info)
        if not data.sort in ["pvalue", "corr"]:
            info = {'success': False, 'info': '选择的展示参数输入错误!', 'code': 'C2300308'}
            return json.dumps(info)
        if project_type == "LC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
        params_json = {
            'metab_table': data.metab_table,
            'metab_set_table': data.metab_set_table,
            "group_detail": json.loads(data.group_detail),
            'group_id': data.group_id,
            'trans_exp_main_id': data.trans_exp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'coefficient': data.coefficient,
            'padjust_method': data.padjust_method,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id,
            'sort': data.sort,
            'top': int(data.top)
        }
        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('relation_corr_network', mongo_data)
        metabolome.insert_main_id('relation_corr_network', main_table_id)
        options = {
            'metab_table': metab_table_path,
            "metab_desc": metab_desc_path,
            'metab_set_table': data.metab_set_table,
            "group_detail": data.group_detail,
            'trans_exp_main_id': data.trans_exp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'coefficient': data.coefficient,
            'padjust_method': data.padjust_method,
            'main_table_id': str(main_table_id),
            "name": name,
            'sort': data.sort,
            'top': int(data.top)
        }
        if table_name == "MetabTable_Origin":
            scale = metabolome.get_scale_type(data.metab_table, task_id)
            if scale and scale == "none":
                options["log10"] = True
            else:
                options["log10"] = False
        elif table_name == "raw":
            options["log10"] = True
        update_info = {str(main_table_id): 'relation_corr_network'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_group_by_detail(group_detail)')
        m_table_name = "Relation/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(RelationCorrNetworkAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)
