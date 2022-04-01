# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON
import datetime


class RocAction(MetaController):
    def __init__(self):
        super(RocAction, self).__init__(instant=False)
        # super(Roc, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'submit_location',
                        'group_detail', 'group_id', 'method_type',
                        'confidence_interval']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) != 2:
            info = {"success": False, "info": "分析只适用于分组方案的分组类别数量为2的情况！", "code":'C2203501'}
            return json.dumps(info)
        key1 = list(table_dict.values())
        for i in range(len(key1)-1):
            if (len(key1[i]) < 3):
                info = {"success": False, "info": "每组内样本数必须大于等于3！", "code":"C2203503"}
                return json.dumps(info)
        if float(data.confidence_interval)<=0 or float(data.confidence_interval)>=1:
            info = {"success": False, "info": "置信区间的值必须为属于(0,1)区间的小数", 'code': 'C2203504'}
            return json.dumps(info)
        task_name = 'meta.report.roc'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2203502'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'method_type': data.method_type,
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'confidence_interval':float(data.confidence_interval),
        }
        if hasattr(data, 'top_n_id'):
            params_json['top_n_id'] = data.top_n_id
        if hasattr(data, 'otuset_id'):
            params_json['otuset_id'] = data.otuset_id
        if hasattr(data, 'norm_strategy'):
            params_json['norm_strategy'] = data.norm_strategy
        if hasattr(data, 'env_id') and data.env_id:
            params_json['env_id'] = data.env_id
            params_json['env_lab'] = data.env_lab
        main_table_name = 'Roc' + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', 'computing'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_none_table('sg_roc')
        update_info = {str(main_table_id): 'sg_roc'}
        options = {
            'otu_table': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'method': data.method_type,
            'top_n': -1,  # 默认 -1 表示没设置此参数
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            #'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'roc_id': str(main_table_id),
            'confidence_interval': float(data.confidence_interval),
            'main_table_data': SON(mongo_data)
        }
        if hasattr(data, 'top_n_id'):
            options['top_n'] = data.top_n_id
        if hasattr(data, 'otuset_id'):
            options['otuset_id'] = data.otuset_id
        if hasattr(data, 'norm_strategy'):
            options['norm_method'] = data.norm_strategy
        if hasattr(data, 'env_id') and data.env_id:
            options['env_id'] = data.env_id
            options['env_labs'] = data.env_lab
        # to_file = ["meta.export_otu_table_by_detail(otu_table)", "meta.export_group_table_by_detail(group_table)"]
        to_file = ["meta.export_otu_env_by_detail(otu_table)", "meta.export_group_table_by_detail(group_table)"]
        print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name="ROC/" + main_table_name,
                            module_type=task_type, to_file=to_file) # modified by hongdongxuan 20170322 在main_table_name前面加上文件输出的文件夹名
        task_info = super(RocAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        # print(self.return_msg)
        return json.dumps(task_info)
