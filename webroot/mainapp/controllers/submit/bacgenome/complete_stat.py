# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class CompleteStatAction(BacgenomeController):
    def __init__(self):
        super(CompleteStatAction, self).__init__(instant=False)
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
        print data
        default_argu = ['task_id', 'sample', 'genome', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.complete_stat'
        module_type = 'workflow'
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        fasta, genome_id = self.bacgenome.get_genome_fill(data.task_id, data.sample, data.genome)
        print(fasta)
        print(genome_id)
        fq1,fq2 = self.bacgenome.get_pe_reads(data.task_id, data.sample)
        params = {
            'task_id': data.task_id,
            'sample': data.sample,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'genome': data.genome
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = data.sample + "_" + data.genome + '_' + 'CompleteStat_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('samp', data.sample),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("genome", data.genome),
            ("genome_id", ObjectId(genome_id))
        ]
        main_id = self.bacgenome.insert_main_table('complete_stat', mongo_data)
        self.bacgenome.insert_main_table_new("complete_stat", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'complete_stat'
        options = {
                   'samp': data.sample,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                    "fa": fasta,
                    'genome_id': str(genome_id),
                    'fq1': fq1,
                    'fq2': fq2
                   }
        to_file = []
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="CompleteStat/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(CompleteStatAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
