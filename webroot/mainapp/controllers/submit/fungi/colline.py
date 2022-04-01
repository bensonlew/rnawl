# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.fungi_genome_controller import FungiGenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class CollineAction(FungiGenomeController):
    def __init__(self):
        super(CollineAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'fungi.report.colline'
        module_type = 'workflow'
        project_sn = self.fungi_genome.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'specimen_id': data.specimen_id,
            'ref_id': data.ref_id,
        }
        r_name = data.ref_id
        if data.ref_id == "user_defined_species":
            if hasattr(data,'ref'):
                params['ref'] = data.ref
                r_name = data.ref
            else:
                params['ref_path'] = data.ref_path
                params['file_dir_id'] = data.file_dir_id
                r_name = data.ref_path.split('/')[-1]

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name =  data.specimen_id+ '_vs_' + r_name + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('specimen_id', data.specimen_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
            ]
        update_info = {}
        main_id = self.fungi_genome.insert_main_table('colline', mongo_data)
        update_info[str(main_id)] = 'colline'
        options = {'task_id': data.task_id,
                   'specimen_id': data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   'main_name': main_table_name,
                   'ref_id': data.ref_id,
                   'middle_path' : self.fungi_genome.get_middle_path(data.task_id)
                   }

        if data.ref_id == "user_defined_species":
            if hasattr(data,'ref'):
                options['genbank'] = data.ref
            else:
                options['ref_path'] = data.ref_path

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="synteny/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params
                            )
        task_info = super(CollineAction, self).POST()

        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
