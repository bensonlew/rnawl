# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class MantelTestAction(MetaasvController):
    """
    mantel_test检验接口
    """
    MATRIX = ['abund_jaccard', 'binary_chisq', 'binary_chord', 'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'binary_lennon', 'binary_ochiai', 'binary_otu_gain', 'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran', 'canberra', 'chisq', 'chord', 'euclidean',
              'gower', 'hellinger', 'kulczynski', 'manhattan', 'morisita_horn', 'pearson', 'soergel', 'spearman_approx',
              'specprof', 'unifrac', 'unweighted_unifrac', 'weighted_normalized_unifrac', 'weighted_unifrac']

    MATRIXFACTOR = ['abund_jaccard', 'binary_chisq', 'binary_chord', 'binary_euclidean', 'binary_hamming',
                    'binary_jaccard', 'binary_lennon', 'binary_ochiai', 'binary_otu_gain', 'binary_pearson',
                    'binary_sorensen_dice', 'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran', 'canberra',
                    'chisq', 'chord', 'euclidean', 'gower', 'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                    'pearson', 'soergel', 'spearman_approx', 'specprof']

    def __init__(self):
        super(MantelTestAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'submit_location', "group_id", "env_id", "otu_method", "env_method", "env_labs"]
        if not hasattr(data, 'env_id'):
            info = {'success': False, 'info': '缺少环境因子参数!'}
            return json.dumps(info)
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)

        if data.otu_method not in self.MATRIX:
            variables = []
            variables.append(data.otu_method)
            info = {'success': False, 'info': '选择计算OTU表矩阵的方法不正确%s!' % data.otu_method, 'variables':variables}
            return json.dumps(info)
        if data.env_method not in self.MATRIXFACTOR:
            variables = []
            variables.append(data.otu_method)
            info = {'success': False, 'info': '选择计算环境因子表样本距离矩阵s的方法不正确%s!' % data.otu_method, 'variables':variables}
            return json.dumps(info)
       
        sample_num = 0
        tmp_group_detail = eval(data.group_detail)
        for k in tmp_group_detail.keys():
            sample_num+=len(tmp_group_detail[k])
        if sample_num < 3 :
            info = {'success': False, 'info': "所选样本总数小于3"}
            return json.dumps(info)


        task_name = 'metaasv.report.mantel_test'
        task_type = 'workflow'

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])

        params_json = {
            "asv_id": data.asv_id,
            "level_id": int(data.level_id),
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "group_detail": group_detail_sort(data.group_detail),
            "group_id": data.group_id,
            "env_id": data.env_id,
            "otu_method": data.otu_method,
            "env_method": data.env_method,
            "env_labs": data.env_labs
        }
        if hasattr(data, "units"):
            params_json["units"] = data.units
            main_table_name = "PartialMantelTest_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        else:
            main_table_name = 'MantelTest_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ("env_id", ObjectId(data.env_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('mantel')
        update_info = {str(main_table_id): 'mantel'}

        options = {
            "otu_file": data.asv_id,
            "env_file": data.env_id,
            'update_info': json.dumps(update_info),
            "level": data.level_id,
            "group_detail": data.group_detail,
            "main_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        del params_json["level_id"]
        del params_json["group_detail"]
        options.update(params_json)
        to_file = ['metaasv.export_otu_table_by_detail(otu_file)', "metaasv_env.export_float_env(env_file)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="MantelTest/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(MantelTestAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)