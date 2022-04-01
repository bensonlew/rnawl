# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON
import datetime


class RandomForestAction(MetaasvController):
    """
    Metaasv Randomforest分析
    """
    def __init__(self):
        super(RandomForestAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'ntree_id', 'method',"data_standard"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        if data.method not in ["AUC", "CV"]:
            variables = []
            variables.append(argu)
            info = {'success': False, 'info': '%s参数不正确!' % argu, 'variables':variables}
            return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) < 2:
            info = {"success": False, "info": "分析只适用于分组方案的分组类别数量大于等于2的情况！"}
            return json.dumps(info)
        task_name = 'metaasv.report.randomforest'
        module_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'ntree_id': data.ntree_id,
            'method': data.method,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            "data_standard": data.data_standard,
        }
        if hasattr(data, "env_id"):
            params_json['env_id']=  data.env_id
        if hasattr(data, "env_labs"):
            params_json['env_labs']=  data.env_labs
        if hasattr(data, "predict_detail"):
            params_json['predict_detail']=  json.loads(data.predict_detail)
        # if hasattr(data, 'group_id2'):
        #     params_json['group_id2'] = data.group_id2
        #     if not hasattr(data, 'group_detail2'):
        #         info = {'success': False, 'info': 'parameters missing:%s' % "group_detail2"}
        #         return json.dumps(info)
        #     else:
        #         params_json['group_detail2'] = group_detail_sort(data.group_detail2)
        main_table_name = 'Randomforest' + \
                          '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('randomforest')
        update_info = {str(main_table_id): 'randomforest'}
        options = {
            'otutable': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'method': data.method,
            'ntree': data.ntree_id,
            "data_standard": data.data_standard,
            'grouptable': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        if hasattr(data, "env_id"):
            options['env_id']=  data.env_id
        if hasattr(data, "env_labs"):
            options['env_labs']=  data.env_labs
        if hasattr(data, 'predict_detail'):
            options['predict_sample'] = data.predict_detail
            # options['group_id2'] = data.group_id2
            # options['group_detail2'] = data.group_detail2
            to_file = ["metaasv.export_otu_table_by_detail(otutable)", "metaasv.export_group_table_by_detail(grouptable)",
                       "metaasv.export_otu_table_by_detail2(predict_sample)"]
        else:
            to_file = ["metaasv.export_otu_table_by_detail(otutable)", "metaasv.export_group_table_by_detail(grouptable)"]
        if hasattr(data, 'env_id'):
            options['envtable']=  data.env_id
            to_file.append("metaasv_env.export_float_env_regression(envtable)")
        self.set_sheet_data(name=task_name, options=options, main_table_name="Randomforest/" + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(RandomForestAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        print(self.return_msg)
        return json.dumps(task_info)
