# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import param_pack, group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class EnterotypingAction(MetaController):
    def __init__(self):  # 20170106 2 lines
        super(EnterotypingAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        # return_info = super(Enterotyping, self).POST()
        # if return_info:
        #     return return_info
        data = web.input()
        postArgs = ["otu_id", "level_id", "group_id", "group_detail", "task_type", "submit_location"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '{}parameters missing!'.format(arg)}
                return json.dumps(info)
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if int(data.level_id) in [1]:
            info = {'success': False, 'info': '样本量或物种分类过低，不能进行菌群分型分析！', 'code':'C2201001'}
            return json.dumps(info)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2201002'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        task_name = 'meta.report.enterotyping'
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        main_table_name = 'Enterotyping_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),  # maybe ObjectId(data.otu_id)
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("show", 0),
            ("type", "otu_Enterotyping")
        ]
        main_table_id = self.meta.insert_none_table('sg_enterotyping')
        update_info = {str(main_table_id): 'sg_enterotyping'}
        options = {
            "input_otu_id": data.otu_id,
            "in_otu_table": data.otu_id,
            "group_detail": data.group_detail,
            "group_id": data.group_id,
            "level": str(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        to_file = "meta.export_otu_table_by_detail(in_otu_table)"
        self.set_sheet_data(name=task_name, options=options, main_table_name="Enterotyping/" + main_table_name,
                            module_type='workflow', to_file=to_file)
        task_info = super(EnterotypingAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
        # insert_data = {
        #     "project_sn": project_sn,
        #     'task_id': self.task_id,
        #     'otu_id': from_otu_table,
        #     # 'cluster_name': cluster_name,
        #     # 'spe_name': spe_name,
        #     'name': name,
        #     "params": params,
        #     'status': 'end',
        #     'desc': 'result after Enterotyping analysis',
        #     'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        #     "show": 0,
        #     "type": "otu_Enterotyping"
        # }


        # self.options = {
        #     "input_otu_id": data.otu_id,
        #     "in_otu_table": data.otu_id,
        #     "group_detail": data.group_detail,
        #     "group_id": data.group_id,
        #     "level": str(data.level_id),
        # }

        # self.to_file = "meta.export_otu_table_by_detail(in_otu_table)"
        # self.to_file = "meta.export_otu_table_by_level(in_otu_table)"    # 暂时不改动同样的方式导表
        # my_param = dict()
        # my_param['submit_location'] = data.submit_location
        # my_param['otu_id'] = data.otu_id
        # my_param["level_id"] = int(data.level_id)
        # my_param["group_id"] = data.group_id
        # my_param['group_detail'] = group_detail_sort(data.group_detail)
        # my_param["task_type"] = data.task_type
        # self.params = param_pack(my_param)
        # self.run()
        # return self.returnInfo
