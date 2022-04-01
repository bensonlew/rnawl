# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class CompleteKmerPcaAction(BacgenomeController):
    def __init__(self):
        super(CompleteKmerPcaAction, self).__init__(instant=False)
        self.bacgenome._project_type = "bac_assem"

    def _update_status_api(self):
        """
        根据client决定接口api为bacgenome.update_status/bacgenome.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'bac_assem.update_status'
        else:
            return 'bac_assem.tupdate_status'

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'sample', 'genome', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.complete_pca'
        module_type = 'workflow'
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        genome_name,fasta = self.bacgenome.get_genome_info(data.genome)
        params = {
            'task_id': data.task_id,
            'sample': data.sample,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'genome': data.genome
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = data.sample + "_" + genome_name + '_' + 'CompleteKmerPca_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("genome", genome_name),
            ("genome_id", ObjectId(data.genome))
        ]
        main_id = self.bacgenome.insert_main_table('complete_kmer_pca', mongo_data)
        self.bacgenome.insert_main_table_new("complete_kmer_pca", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'complete_kmer_pca'
        options = {
                   'samp': data.sample,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "fa": fasta,
                   'genome_id': str(data.genome)
                   }
        to_file = []
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="CompleteKmerPca/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(CompleteKmerPcaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
