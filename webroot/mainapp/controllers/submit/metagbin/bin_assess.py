# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metagbin_controller import MetagbinController
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.libs.signature import check_sig
from bson import SON


class BinAssessAction(MetagbinController):
    def __init__(self):
        super(BinAssessAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'bin_id', 'new_bin', 'bin_list', 'submit_location', 'task_type', 'path']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'metagbin.report.bin_assess'
        module_type = 'workflow'
        project_sn = self.metagbin.get_projectsn(data.task_id)
        params = {
            'bin_list': data.bin_list,
            'bin_id': data.bin_id,
            'new_bin': data.new_bin,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'ReBin_' + data.new_bin + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('new_bin',data.new_bin),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.metagbin.insert_main_table('bin',mongo_data)
        Metagbin().insert_main_table_new("bin", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'bin'
        options = {
                   'bin_path': data.path,
                   'update_info': json.dumps(update_info),
                   'bin_list': data.bin_list,
                   'new_bin': data.new_bin,
                   'main_id': str(main_id),
                   }
        self.set_sheet_data(name=task_name,options=options,
                            main_table_name= "ReBin/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(BinAssessAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)