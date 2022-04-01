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



class RocAction(MetaasvController):
    """
    Metaasv Roc分析
    """
    def __init__(self):
        super(RocAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'method', 'confidence_interval', "data_standard"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) != 2:
            info = {"success": False, "info": "分析只适用于分组方案的分组类别数量为2的情况！"}
            return json.dumps(info)
        key1 = list(table_dict.values())
        for i in range(len(key1)):
            if (len(key1[i]) < 3):
                info = {"success": False, "info": "每组内样本数必须大于等于3！"}
                return json.dump(info)
        if float(data.confidence_interval)<=0 or float(data.confidence_interval)>=1:
            info = {"success": False, "info": "置信区间的值必须为属于(0,1)区间的小数"}
            return json.dumps(info)
        task_name = 'metaasv.report.roc'
        task_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'method': data.method,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'confidence_interval':float(data.confidence_interval),
            "data_standard": data.data_standard
        }
        if hasattr(data, "top_n"):
            params_json['top_n']=  int(data.top_n)
        if hasattr(data, "set_id"):
            params_json['set_id']=  data.set_id
        if hasattr(data, "env_id"):
            params_json['env_id']=  data.env_id
        if hasattr(data, "env_labs"):
            params_json['env_labs']=  data.env_labs
        main_table_name = 'Roc' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', 'computing'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_none_table('roc')
        update_info = {str(main_table_id): 'roc'}
        options = {
            'otu_table': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'method': data.method,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'data_standard': data.data_standard,
            'confidence_interval': float(data.confidence_interval),
            'main_table_data': SON(mongo_data)
        }
        if hasattr(data, "top_n"):
            options['top_n']=  int(data.top_n)
        if hasattr(data, "set_id"):
            options['set_id']=  data.set_id
        if hasattr(data, "env_id"):
            options['env_id']=  data.env_id
        if hasattr(data, "env_labs"):
            options['env_labs']=  data.env_labs
        to_file = ["metaasv.export_otu_table_by_detail(otu_table)", "metaasv.export_group_table_by_detail(group_table)"]
        if hasattr(data, 'env_id'):
            options['envtable']=  data.env_id
            to_file.append("metaasv_env.export_float_env_regression(envtable)")
        self.set_sheet_data(name=task_name, options=options, main_table_name="ROC/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(RocAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
