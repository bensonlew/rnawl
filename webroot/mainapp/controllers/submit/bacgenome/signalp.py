# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class SignalpAction(BacgenomeController):
    def __init__(self):
        super(SignalpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'specimen_id', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.signalp'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        mytype = "bac"  # 新增，同时跑两个
        params = {
            'task_id': data.task_id,
            #'type': mytype,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        spe_list = data.specimen_id.split(",")
        update_info = {}
        main_id = []
        main_name = []
        #for i in spe_list:
            #main_table_name = 'SignalP_' + i + '_' + data.type.capitalize() + \
            #                  '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = 'SignalP_'+ datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('type', mytype),
            ('specimen_id', data.specimen_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_name.append(str(main_table_name))
            # main_table_id = self.bacgenome.insert_none_table('signalp')
        main_table_id = self.bacgenome.insert_main_table('anno_signalp', mongo_data)
        update_info[str(main_table_id)] = 'anno_signalp'
        main_id.append(str(main_table_id))
        query = self.bacgenome.get_sample_genefile(data.task_id, data.specimen_id, type='faa')
        options = {'query': query,
                   'type': mytype,
                   'task_id': data.task_id,
                   'specimen_id': data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'main_name': ','.join(main_name),  #
                   'params': params,
                   'main_id': ','.join(main_id)
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="SignalP",
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(SignalpAction, self).POST()
        # return spe_list
        if task_info['success']:
            task_info['content'] = {'ids': {'id': main_id, 'name': main_name}}
        return json.dumps(task_info)
