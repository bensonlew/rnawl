
# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.meta import Meta
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig



class ProcrustesAction(MetaController):
    def __init__(self):
        super(ProcrustesAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = [ 'submit_location', 'group_detail','task_type',
                        'group_id', 'env_id', 'env_labs', 'otu_id','level_id','sort_method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu], "code" : "C2204201"}
                return json.dumps(info)

        dist_list = ['abund_jaccard', 'binary_chisq', 'binary_chord',
                  'binary_euclidean', 'binary_hamming', 'binary_jaccard',
                  'binary_lennon', 'binary_ochiai',
                  'binary_pearson', 'binary_sorensen_dice',
                  'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
                  'canberra', 'chisq', 'chord', 'euclidean', 'gower',
                  'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                  'pearson', 'soergel', 'spearman_approx', 'specprof',
                  'unifrac','unweighted_unifrac', 'unweighted_unifrac_full_tree',
                  'weighted_normalized_unifrac', 'weighted_unifrac']

        # dist_list = ['bray_curits', 'euclidean']

        if hasattr(data, "meta_method") and data.meta_method not in dist_list:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of  meta_method', "code" : "C2204202"}
            return json.dumps(info)
        if hasattr(data, "spe_method") and data.spe_method not in dist_list:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of spe_mathod ', "code" : "C2204203"}
            return json.dumps(info)

        task_name = 'meta.report.procrustes'
        module_type = 'workflow'
        # task_type = 2
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", "code" : "C2204204"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'sort_method' : data.sort_method,
        }

        if hasattr(data,'spe_method'):
            params_json['spe_method'] = data.spe_method
        if hasattr(data, 'meta_method'):
            params_json['meta_method'] = data.meta_method


        main_table_name = 'Procrustes_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("columns",{'PC1':['p1','c1'],'PC2': ['p2','c2']})
        ]
        main_table_id = self.meta.insert_main_table('sg_procrustes', mongo_data)
        update_info = {str(main_table_id): 'sg_procrustes'}
        options = {
            'otu_table': data.otu_id,
            'otu_id': data.otu_id,
            'level': data.level_id,
            'asso_table': data.env_id,
            'env_labs': data.env_labs,
            #'env_detail' : data.env_labs,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'group_id': data.group_id,
            'update_info': json.dumps(update_info),
            'main_table_id': str(main_table_id),
            'method' : data.sort_method
        }
        if hasattr(data,'asso_dist'):
            options['asso_dist'] = data.meta_method
        if hasattr(data, 'otu_dist'):
            options['otu_dist'] = data.spe_method

        to_file = ["meta.export_otu_table_by_level(otu_table)",
                    "meta.export_group_table_by_detail(group_table)", "env.export_float_env_regression(asso_table)"]

        self.set_sheet_data(name=task_name, options=options, main_table_name="Procrustes/"+main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(ProcrustesAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)