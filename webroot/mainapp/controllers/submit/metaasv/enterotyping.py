# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import param_pack, group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class EnterotypingAction(MetaasvController):
    """
    metaasv 样本菌群分析
    """
    def __init__(self):
        super(EnterotypingAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ["asv_id", "level_id", "group_id", "group_detail", "task_type", "submit_location", "distance_method"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '{}parameters missing!'.format(arg)}
                return json.dumps(info)
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if int(data.level_id) in [1]:
            info = {'success': False, 'info': '样本量或物种分类过低，不能进行菌群分型分析！'}
            return json.dumps(info)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        task_name = 'metaasv.report.enterotyping'
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            "distance_method": data.distance_method
        }
        main_table_name = 'Enterotyping_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]
        main_table_id = self.metaasv.insert_none_table('enterotype')
        update_info = {str(main_table_id): 'enterotype'}
        options = {
            "asv_id": data.asv_id,
            "otu_file": data.asv_id,
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            "level": int(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'main_table_data': SON(mongo_data),
            "distance_method": data.distance_method
        }
        to_file = ["metaasv.export_otu_table_by_detail(otu_file)", "metaasv.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Enterotyping/" + main_table_name,
                            module_type='workflow', to_file=to_file)
        task_info = super(EnterotypingAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
