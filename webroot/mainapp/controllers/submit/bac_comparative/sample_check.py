# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web,re
import json
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
import datetime

class SampleCheckAction(BacComparativeController):
    def __init__(self):
        super(SampleCheckAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'mapping_file', "main_id", "type", "project_sn", "member_id"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.sample_check'
        module_type = 'workflow'  # 可以不配置
        update_info = {}
        main_table_name = 'SampleCheck_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info[data.main_id] = 'sample'
        options = {
            'task_id': data.task_id,
            'update_info': json.dumps(update_info),
            'main_id': data.main_id,
            "mapping_file": data.mapping_file,
            'type': data.type
                   }
        self.set_sheet_data2(name=task_name,
                            options = options,
                            main_table_name = "SampleCheck/" + main_table_name,
                            module_type = module_type,
                            project_sn = data.project_sn,
                            task_id = data.task_id)
        task_info = super(SampleCheckAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(data.main_id), 'name': "sample_check"}}
        return json.dumps(task_info)