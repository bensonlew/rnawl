# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class ParalogsAction(BacgenomeController):
    def __init__(self):
        super(ParalogsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id',  'submit_location', 'task_type'] #'search_id',
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.paralogs'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        new_species = self.bacgenome.get_specimen_name(data.task_id, data.specimen_id)
        params = {
            'task_id': data.task_id,
            #'search_id': data.search_id,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Paralogs_' + data.specimen_id +  \
                          '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            #('search_id', data.search_id),
            ('specimen_id', data.specimen_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("version",'2.0')
        ]
        main_id = self.bacgenome.insert_main_table('paralogs', mongo_data)
        update_info[str(main_id)] = 'paralogs'
        sequence = self.bacgenome.get_sample_genefile(data.task_id, data.specimen_id,'faa') #guanqing.zou 20180912
        options = {
                    #'gene_file': data.search_id,
                   'sequence': sequence,
                   'task_id': data.task_id,
                   #'search_id': data.search_id,
                   'specimen_id': data.specimen_id,
                   # 'old_species':data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'params': params,
                   'main_id': str(main_id),
                   'main_table_data': SON(mongo_data),
                   }
        #to_file = ['bacgenome.export_gene_faa_by_geneid(gene_file)']   # zouguanqing 20180912
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Paralogs/"+main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            )
        task_info = super(ParalogsAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
