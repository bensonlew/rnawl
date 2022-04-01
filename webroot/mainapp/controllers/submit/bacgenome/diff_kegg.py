# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON
from biocluster.config import Config

class DiffKeggAction(BacgenomeController):
    def __init__(self):
        super(DiffKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id','task_type','submit_location', "pathway_id"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.diff_kegg'
        module_type = 'workflow'
        task_info = self.bacgenome.get_task_info(data.task_id)
        project_sn = task_info['project_sn']
        # new_species = []
        # for s_id in  data.specimen_id.split(","):
        #     new_species.append(self.bacgenome.get_specimen_name(data.task_id, s_id))
        # new_species = ','.join(new_species)
        new_species = data.specimen_id
        params = {
            #'task_id': data.task_id,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'pathway_id':data.pathway_id
        }
        if hasattr(data,'pathway_des'):
            params['pathway_des']=data.pathway_des
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'DiffKegg_' +  datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
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
        main_id = self.bacgenome.insert_main_table('diff_kegg', mongo_data)
        update_info[str(main_id)] = 'diff_kegg'
        pathway_id = data.pathway_id.replace('map','ko')
        ks = self.bacgenome.get_all_common_info('kegg_ko',{'pathway_id':pathway_id},ref=True)
        if ks:
            ks_str = ';'.join([ i['ko_id'] for i in ks])
        else:
            info = {'success': False, 'info': 'K not found' }
            return json.dumps(info)
        options = {
                    'samples': new_species,
                    'update_info': json.dumps(update_info),
                    'main_table_id': str(main_id),
                    'pathway_id': data.pathway_id,
                    'k_list': ks_str
                   }

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="DiffKegg/"+main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            )
        task_info = super(DiffKeggAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
