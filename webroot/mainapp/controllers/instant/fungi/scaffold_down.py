# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web
import json, os, re, time,datetime
from mainapp.controllers.project.fungi_genome_controller import FungiGenomeController
from mainapp.libs.signature import check_sig


class ScaffoldDownAction(FungiGenomeController):
    def __init__(self):
        super(ScaffoldDownAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id','type','specimen_name','seq_list','path','client', 'seq_type','task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'fungi.report.scaffold_down'
        module_type = 'workflow'
        project_sn = self.fungi_genome.get_projectsn(data.task_id)

        params_json = {
            'specimen_name': data.specimen_name,
            'seq_list': data.seq_list,
            'type':data.type,
            'task_type': int(data.task_type),
            'task_id': data.task_id,
            'seq_type':data.seq_type
        }
        main_table_name = 'SeqDown_' + data.type.capitalize() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在下载'),
            ('name', main_table_name),
            ('path',data.path),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_id = self.fungi_genome.insert_main_table('seq_down', mongo_data)
        options = {
            'path': data.path,
            'type': data.type,
            'main_id':str(main_id),
            'main_table_name':main_table_name,
            'specimen_name': str(data.specimen_name),
            'client': data.client,
            'seq_list': str(data.seq_list),
            'task_id': data.task_id,
            'seq_type' : data.seq_type
        }
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name='Seq_down/' + main_table_name,
                            task_id=data.task_id,
                            project_sn=project_sn, module_type=module_type, params=params_json)
        task_info = super(ScaffoldDownAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)

