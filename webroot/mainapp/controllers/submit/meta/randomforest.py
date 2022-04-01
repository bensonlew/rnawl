# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
# last_modify = 2018.04.23 by zhujuan 新增十折交叉验证、AUC验证和预测样本分类功能
import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON
import datetime


class RandomforestAction(MetaController):
    def __init__(self):
        super(RandomforestAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'ntree_id', 'method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        if data.method not in ["AUC", "CV"]:
            variables = []
            variables.append(argu)
            info = {'success': False, 'info': '%s参数不正确!' % argu, 'code':'C2203201', 'variables':variables}
            return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) < 2:
            info = {"success": False, "info": "分析只适用于分组方案的分组类别数量大于等于2的情况！", 'code':'C2203202'}
            return json.dumps(info)
        task_name = 'meta.report.randomforest'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2203203'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'ntree_id': data.ntree_id,
            'method': data.method,
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        if hasattr(data, 'group_id2'):
            params_json['group_id2'] = data.group_id2
            if not hasattr(data, 'group_detail2'):
                info = {'success': False, 'info': 'parameters missing:%s' % "group_detail2"}
                return json.dumps(info)
            else:
                params_json['group_detail2'] = group_detail_sort(data.group_detail2)
        if hasattr(data, 'norm_strategy'):
            params_json['norm_strategy'] = data.norm_strategy
        if hasattr(data, 'env_id') and data.env_id:
            params_json['env_id'] = data.env_id
            params_json['env_lab'] = data.env_lab
        main_table_name = 'Randomforest' + \
                          '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            # ('table_type', 'dist'),
            # ('tree_type', 'cluster'),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_none_table('sg_randomforest')
        update_info = {str(main_table_id): 'sg_randomforest'}
        options = {
            'otutable': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'method': data.method,
            'ntree': data.ntree_id,
            'grouptable': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'randomforest_id': str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        if hasattr(data, 'norm_strategy'):
            options['norm_method'] = data.norm_strategy
        if hasattr(data, 'env_id') and data.env_id:
            options['env_id'] = data.env_id
            options['env_labs'] = data.env_lab
        if hasattr(data, 'group_id2'):
            options['predict_sample'] = data.otu_id
            options['group_id2'] = data.group_id2
            options['group_detail2'] = data.group_detail2
            # to_file = ["meta.export_otu_table_by_detail(otutable)", "meta.export_group_table_by_detail(grouptable)",
            #            "meta.export_otu_table_by_detail2(predict_sample)"]
            to_file = ["meta.export_otu_env_by_detail(otutable)", "meta.export_group_table_by_detail(grouptable)",
                       "meta.export_otu_env_table_by_detail2(predict_sample)"]
        else:
            # to_file = ["meta.export_otu_table_by_detail(otutable)", "meta.export_group_table_by_detail(grouptable)"]
            to_file = ["meta.export_otu_env_by_detail(otutable)", "meta.export_group_table_by_detail(grouptable)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Randomforest/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(RandomforestAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        print(self.return_msg)
        return json.dumps(task_info)
