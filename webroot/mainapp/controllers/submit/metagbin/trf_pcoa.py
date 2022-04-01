# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.metagbin_controller import MetagbinController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metagbin import Metagbin
from bson import SON


class TrfPcoaAction(MetagbinController):
    def __init__(self):
        super(TrfPcoaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'bin_id', 'submit_location', 'task_type','path']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'metagbin.report.trf_pcoa'
        module_type = 'workflow'
        project_sn = self.metagbin.get_projectsn(data.task_id)
        params = {
            'bin_id': data.bin_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'TNFanalysis_' + data.bin_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.metagbin.insert_main_table('tnf',mongo_data)
        Metagbin().insert_main_table_new("tnf", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'tnf'
        options = {
                   'input_genome': data.path,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "bin_id": data.bin_id,
                   }
        self.set_sheet_data(name=task_name,options=options,
                            main_table_name="TNFanalysis/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(TrfPcoaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)