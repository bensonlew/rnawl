# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId


class MetaSourcetrackerAction(MetaController):
    def __init__(self):
        super(MetaSourcetrackerAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['otu_id', 'level_id', 'group_id', 'group_detail', 'second_group_id', 'second_group_detail', 'add_Algorithm', 's', 'submit_location']
        for param in params_name:
            if not hasattr(data, param):
                info = {"success": False, "info": "parameters missing:%s" % param}
                return json.dumps(info)
        if int(data.level_id) not in range(1, 10):
            variables = []
            variables.append(data.level_id)
            info = {"success": False, "info": "level{}不在规定范围内".format(data.level_id), 'code':'C2202001', 'variables':variables}
            return json.dumps(info)
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            success.append("传入的group_detail不是一个字典")
        second_group_detail = json.loads(data.second_group_detail)
        if not isinstance(second_group_detail, dict):
            success.append("传入的second_group_detail不是一个字典")
        otu_info = Meta().get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2202002'}
            return json.dumps(info)
        task_name = 'meta.report.meta_sourcetracker'
        task_type = 'workflow'
        task_info = self.meta.get_task_info(otu_info['task_id'])
        main_table_name = 'MetaSourcetracker_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            'otu_id': data.otu_id,
            'level_id': data.level_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'second_group_id': data.second_group_id,
            'second_group_detail': group_detail_sort(data.second_group_detail),
            's': data.s,
            'add_Algorithm': data.add_Algorithm,
            'submit_location': data.submit_location,
            'task_type': 'reportTask'
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),  # maybe data.otu_id
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('status', 'start'),
            ('desc', 'meta_sourcetracker分析'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_id = self.meta.insert_main_table('sg_sourcetracker', mongo_data)
        update_info = {str(main_table_id): 'sg_sourcetracker'}
        options = {
            "in_otu_table": data.otu_id,
            "level": data.level_id,
            "map_detail": data.group_id,
            "group_detail": data.group_detail,
            "second_group_detail": data.second_group_detail,
            "add_Algorithm": data.add_Algorithm,
            "s": data.s,
            "meta_sourcetracker_id": str(main_table_id),
            "update_info": json.dumps(update_info),
        }
        to_file = ["meta.export_otu_table_by_level(in_otu_table)", "meta.export_group_table_by_detail_2(map_detail)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="MetaSourcetracker/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(MetaSourcetrackerAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        return json.dumps(task_info)
