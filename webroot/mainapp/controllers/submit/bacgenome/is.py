# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class IsAction(BacgenomeController):
    def __init__(self):
        super(IsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'submit_location', 'task_type', 'specimen_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.is_predict'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        update_info = {}
        main_table_name = 'Is_Predict_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('software','ISEScan_v1.7.2.1'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_table_id = self.bacgenome.insert_main_table('is', mongo_data)
        update_info[str(main_table_id)] = 'is'
        options = {
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_table_id),
                   'task_id': data.task_id,
                   "specimen_id": str(data.specimen_id)
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="_".join(main_table_name.strip().split("_")[0:2]) + '/' + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(IsAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
