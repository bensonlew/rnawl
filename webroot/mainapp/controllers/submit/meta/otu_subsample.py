# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import param_pack, group_detail_sort, filter_json_sort
from mainapp.libs.signature import check_sig
from bson import SON
from collections import OrderedDict


class OtuSubsampleAction(MetaController):

    def __init__(self):
        super(OtuSubsampleAction, self).__init__(instant=False)
    
    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ['size', 'submit_location', "otu_id", "task_type", "group_detail", "group_id", "filter_json"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)

        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2202601'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        task_type = 'workflow'
        main_table_name = 'OTUTaxonAnalysis_' + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        my_param = dict()
        my_param["group_id"] = data.group_id
        my_param['otu_id'] = data.otu_id
        my_param["submit_location"] = data.submit_location
        my_param["size"] = data.size
        my_param["filter_json"] = json.loads(data.filter_json,object_pairs_hook=OrderedDict)
        my_param["group_detail"] = group_detail_sort(data.group_detail)
        my_param["task_type"] = data.task_type
        params = param_pack(my_param)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('from_id', str(data.otu_id)),
            ('name', main_table_name),
            ("params", params),
            ('status', 'start'),
            ("level_id", json.dumps([9])),
            ('desc', '正在计算'),
            ("type", "otu_statistic"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]

        main_table_id = self.meta.insert_main_table('sg_otu', mongo_data)
        update_info = {str(main_table_id): 'sg_otu'}
        options = {
            "in_otu_table": data.otu_id,
            "input_otu_id": data.otu_id,
            'group_id': data.group_id,
            "group_detail": data.group_detail,
            "filter_json": data.filter_json,
            "level": "9",
            "size": data.size,
            'update_info': json.dumps(update_info),
            "params": params,
            'main_id': str(main_table_id),
            # 'main_table_data': SON(mongo_data)
        }
        to_file = "meta.export_otu_table_by_level_remove_blank(in_otu_table)" ##fix by qingchen.zhang@20191120
        task_name = 'meta.report.otu_subsample'
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name="OTUTaxonAnalysis/" + main_table_name, module_type=task_type,
                            to_file=to_file)  # modified by hongdongxuan 20170322 前面加上文件输出的文件夹名
        task_info = super(OtuSubsampleAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
