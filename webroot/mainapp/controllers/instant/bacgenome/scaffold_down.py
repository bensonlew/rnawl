# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web
import json, os, re,datetime
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig


class ScaffoldDownAction(BacgenomeController):
    def __init__(self):
        super(ScaffoldDownAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'type', 'specimen_name', 'seq_list', 'client', 'task_type', 'seq_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.scaffold_down'
        module_type = 'workflow'
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        print data.specimen_name
        params_json = {
            'specimen_name': data.specimen_name,
            'seq_list': data.seq_list,
            'type': data.type,
            'task_type': int(data.task_type),
            'task_id': data.task_id,
            'seq_type': data.seq_type
        }
        main_table_name = 'SeqDown_' + data.type.capitalize() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        file_seq = self.bacgenome.get_assemble_bysample(data.task_id, data.specimen_name, data.type, data.seq_type)
        print file_seq
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('seq_path', file_seq),
            ('task_id', data.task_id),
            ('desc', '正在下载'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_id = self.bacgenome.insert_main_table('seq_down', mongo_data)
        options = {
            'type': data.type,
            'seq_path': str(file_seq),
            'main_id':str(main_id),
            'main_table_name': main_table_name,
            'specimen_name': str(data.specimen_name),
            'client': data.client,
            'seq_list': str(data.seq_list),
            'task_id': data.task_id,
            'seq_type': data.seq_type
        }
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name='Seq_down/' + main_table_name,
                            task_id=data.task_id,
                            project_sn=project_sn, module_type=module_type, params=params_json)
        task_info = super(ScaffoldDownAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)