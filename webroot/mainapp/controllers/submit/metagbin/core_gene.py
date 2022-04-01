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


class CoreGeneAction(MetagbinController):
    def __init__(self):
        super(CoreGeneAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id',  'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        if not hasattr(data, 'bins') and  hasattr(data, 'bin_list'):
            info = {'success': False, 'info': 'bins和bin_list必须同时存在'}
            return json.dumps(info)
        if not hasattr(data, 'genomes') and  hasattr(data, 'genome_list'):
            info = {'success': False, 'info': 'genomes和genome_list必须同时存在'}
            return json.dumps(info)
        if hasattr(data, 'bins') and not hasattr(data, 'bin_list'):
            info = {'success': False, 'info': 'bins和bin_list必须同时存在'}
            return json.dumps(info)
        if hasattr(data, 'genomes') and not hasattr(data, 'genome_list'):
            info = {'success': False, 'info': 'genomes和genome_list必须同时存在'}
            return json.dumps(info)
        task_name = 'metagbin.report.core_gene'
        module_type = 'workflow'
        project_sn = self.metagbin.get_projectsn(data.task_id)
        params = {
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
        }
        if hasattr(data, 'bins'):
            params['bins']=data.bins
        if hasattr(data, 'genomes'):
            params['genomes'] = data.genomes
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'CoreGene_' + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.metagbin.insert_main_table('core_gene',mongo_data)
        Metagbin().insert_main_table_new("core_gene", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'core_gene'
        options = {
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   }
        if hasattr(data, 'bin_list'):
            options['bin_list']= data.bin_list
        if hasattr(data, 'genome_list'):
            options['genome_list'] = data.genome_list
        self.set_sheet_data(name=task_name,options=options,
                main_table_name="CoreGene/" + main_table_name,
                module_type=module_type,
                project_sn=project_sn,
                task_id=data.task_id,
                params=params)
        task_info = super(CoreGeneAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)