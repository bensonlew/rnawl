# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class PipeAction(MetaController):
    def __init__(self):
        super(PipeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['otu_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'group_info', 'sub_analysis']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "parameters missing:%s" % param}
                return json.dumps(info)
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            info = {'success': False, 'info': '传入的group_detail不是一个字典！', 'code':'C2203001'}
            return json.dumps(info)
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!", 'code':'C2203002'}
            return json.dumps(info)
        task_name = 'meta.report.meta_pipeline'
        task_type = 'workflow'
        task_info = self.meta.get_task_info(otu_info['task_id'])
        main_table_name = 'Pipeline_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        all_analysis = []
        sub_analysis_name = json.loads(data.sub_analysis)
        for key in sub_analysis_name:
            all_analysis.append(str(key))
        if hasattr(data, 'env_id') and hasattr(data, 'env_labs'):
            params_json = {
                'otu_id': data.otu_id,
                'level_id': data.level_id,
                'group_id': data.group_id,
                'group_detail': group_detail_sort(data.group_detail),
                'submit_location': data.submit_location,
                'env_id': data.env_id,
                'env_labs': data.env_labs,
                'group_info': json.loads(data.group_info),
                'filter_json': json.loads(data.filter_json),
                'size': data.size,
                'task_type': data.task_type,
                'first_group_id': data.first_group_id,
                'sub_analysis': json.loads(data.sub_analysis)
            }
        else:
            params_json = {
                'otu_id': data.otu_id,
                'level_id': data.level_id,
                'group_id': data.group_id,
                'group_detail': group_detail_sort(data.group_detail),
                'submit_location': data.submit_location,
                'group_info': json.loads(data.group_info),
                'filter_json': json.loads(data.filter_json),
                'size': data.size,
                'task_type': data.task_type,
                'first_group_id': data.first_group_id,
                'sub_analysis': json.loads(data.sub_analysis)
            }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'pipe_analysis.....'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('submit_location', data.submit_location),
            ('otu_id', ObjectId(data.otu_id)),
            ('analysis_list', str(all_analysis)),
            ('ends_count', 0),
            ('all_count', 0),
            ('pipe_main_id', "none")
        ]
        main_table_id = self.meta.insert_none_table('sg_pipe_batch')
        update_info = {str(main_table_id): 'sg_pipe_batch'}
        options = {
            "data": json.dumps(data),
            "update_info": json.dumps(update_info),
            "pipe_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                            module_type=task_type)
        # self._sheet_data["instant"] = True #投递的接口，将workflow改为
        task_info = super(PipeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }
            }
        return json.dumps(task_info)
