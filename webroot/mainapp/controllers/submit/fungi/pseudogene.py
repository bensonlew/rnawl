# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.fungi_genome_controller import FungiGenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class PseudogeneAction(FungiGenomeController):
    def __init__(self):
        super(PseudogeneAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_list', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'fungi.report.pseudogene'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.fungi_genome.get_projectsn(data.task_id)
        params = {
            'task_id': data.task_id,
            #'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'specimen_list': data.specimen_list
        }
        if hasattr(data, 'ref'):
            params['ref'] = data.ref
            ref_name = data.ref
        if  hasattr(data, 'ref_path'):
            params['ref_path'] = data.ref_path
            ref_name = data.ref_path.split('/')[-1].split('.')[0]
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        spe_list = data.specimen_list.split(",")

        main_id = []
        main_name = []
        update_info = {}
        for i in spe_list:
            main_table_name = i + '_' + ref_name + '_' + 'Pseudogene_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_params = {
                'task_id': data.task_id,
                'specimen_id': i,
                'project_sn': project_sn,
                'submit_location': data.submit_location,
                'task_type': data.task_type
            }
            if hasattr(data, 'ref'):
                mongo_params['ref'] = data.ref
            if  hasattr(data, 'ref_path'):
                mongo_params['ref_path'] = data.ref_path
            mongo_params = json.dumps(mongo_params, sort_keys=True, separators=(',', ':'))
            mongo_data = [
                ('project_sn', project_sn),
                ('status', 'start'),
                ('task_id', data.task_id),
                ('specimen_id', i),
                ('desc', '正在计算'),
                ('name', main_table_name),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", mongo_params)
            ]
            main_name.append(str(main_table_name))
            # main_table_id = self.bacgenome.insert_none_table('signalp')
            main_table_id = self.fungi_genome.insert_main_table('pseudogene', mongo_data)

            update_info[str(main_table_id)] = 'pseudogene'
            main_id.append(str(main_table_id))



        options = {'task_id': data.task_id,
                   'specimen_list': data.specimen_list,
                   #'old_species':data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'params': params,
                   'main_id': ','.join(main_id),
                   'main_name': ','.join(main_name),
                   'middle_path' : self.fungi_genome.get_middle_path(data.task_id)
                   }
        if hasattr(data, 'ref'):
            options['ref'] = data.ref
        if hasattr(data, 'ref_path'):
            options['ref_path'] = data.ref_path


        self.set_sheet_data(name=task_name,
                            options=options,
                            #main_table_name="Pseudogene/" + main_table_name,
                            main_table_name="Pseudogene",
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params
                            #to_file=to_file
                            )
        task_info = super(PseudogeneAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': main_id, 'name': main_name}}
        return json.dumps(task_info)
