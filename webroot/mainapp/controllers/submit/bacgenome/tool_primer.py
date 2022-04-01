# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class ToolPrimerAction(BacgenomeController):
    """
    细菌引物打通接口
    """
    def __init__(self):
        super(ToolPrimerAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'submit_location', 'task_type', 'specimen_id','seq_id','seq_id_new','updownstream', "tool_type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.tool_primer'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'seq_id': data.seq_id,
            'seq_id_new':data.seq_id_new,
            'updownstream':data.updownstream,
            'specimen_id':data.specimen_id,
            'tool_type':data.tool_type
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        update_info = {}
        main_table_name = 'Tool_Primer_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_table_id = self.bacgenome.insert_main_table('sg_tool_lab_primer', mongo_data)
        update_info[str(main_table_id)] = 'sg_tool_lab_primer'
        options = {
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'task_id': data.task_id,
            "specimen_id": str(data.specimen_id),
            "upstream": int(data.updownstream),
            "downstream": int(data.updownstream),
            "relate_name": main_table_name,
            "params": params,
            "tool_type": data.tool_type,
            "seq_id":data.seq_id,
            "seq_id_new":data.seq_id_new
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(ToolPrimerAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
