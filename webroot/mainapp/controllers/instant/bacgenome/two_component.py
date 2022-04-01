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
from biocluster.file import download, exists

class TwoComponentAction(BacgenomeController):
    def __init__(self):
        super(TwoComponentAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id','task_type','submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.two_component'
        module_type = 'workflow'  # 可以不配置
        #project_sn = self.bacgenome.get_projectsn(data.task_id)
        #task_info = self.bacgenome.get_task_info(data.task_id)
        task_info = self.bacgenome.get_one_common_info("sg_task",{"task_id":data.task_id})
        project_sn = task_info['project_sn']
        sanger_prefix = Config().get_project_region_bucket(project_type="bacgenome")

        pfam_anno_list = []
        new_species_list = []
        for s in data.specimen_id.split(','):
            #pfam_anno = sanger_prefix + 'files/'+ str(task_info['member_id']) +  '/' + str(task_info['project_sn']) + '/' + \
            #  data.task_id + '/workflow_results/' + s + '/annotation/Pfam/%s_anno_pfam.xls'%s
            ret = self.bacgenome.get_one_common_info("gene_predict",{"task_id":data.task_id})
            pfam_anno = ret['file_path'][0]+ '/'+ s + '/annotation/Pfam/%s_anno_pfam.xls'%s
            if not exists(pfam_anno):
                pfam_anno = ret['file_path'][0]+ '/'+ s + '/annotation/Pfam/%s_whole_genome_anno_pfam.xls'%s

            pfam_anno_list.append(pfam_anno)
            #new_species = self.bacgenome.get_specimen_name(data.task_id, s)
            new_species = s
            new_species_list.append(new_species)
        pfam_anno = ','.join(pfam_anno_list)
        new_species = ','.join(new_species_list)

        params = {
            'task_id': data.task_id,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'TwoComponent_' + data.specimen_id + '_' +  datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'end'),
            ('task_id', data.task_id),
            ('specimen_id', data.specimen_id),
            ('desc', '完成计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bacgenome.insert_main_table('anno_regulator', mongo_data)
        update_info[str(main_id)] = 'anno_regulator'

        options = {
                    'sample_name': new_species,
                    'update_info': json.dumps(update_info),
                    'main_id': str(main_id),
                    'pfam_anno' : pfam_anno
                   }

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="TwoComponent/"+main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            )
        task_info = super(TwoComponentAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
