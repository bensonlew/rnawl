# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON
import web


class CgviewAction(BacgenomeController):
    def __init__(self):
        super(CgviewAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id', 'submit_location', 'task_type','species_name','genome_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.cgview'
        module_type = 'workflow'
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        if not hasattr(data, "seq_type"):
            seq_type = "Circular"
        else:
            seq_type = data.seq_type
        params = {
            'task_id': data.task_id,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'species_name':data.species_name,
            'genome_type':data.genome_type,
            'seq_type':seq_type,
        }
        if data.genome_type == 'scaffold':
            xml_file = self.bacgenome.get_xml_file(data.task_id,data.specimen_id,'scaffold')
        else:
            xml_file = self.bacgenome.get_xml_file(data.task_id, data.specimen_id,data.genome_type)
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if data.genome_type == 'scaffold':
            main_table_name = data.specimen_id + '_' + 'cgview_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        else:
            main_table_name = data.specimen_id + '_' + 'cgview_' + data.genome_type + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('version', "3.0"),
            ("params", params)
        ]
        main_id = self.bacgenome.insert_main_table('cgview', mongo_data)
        update_info[str(main_id)] = 'cgview'
        options = {'task_id': data.task_id,
                   'specimen_id': data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'species_name':data.species_name,
                   'params': params,
                   'main_id': str(main_id),
                   "xml_file": xml_file,
                   'main_table_data': SON(mongo_data),
                   'genome_type':data.genome_type,
                   'seq_type': seq_type
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="cgview/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(CgviewAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)