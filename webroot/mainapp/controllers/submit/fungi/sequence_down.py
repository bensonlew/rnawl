# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import web,datetime
import json, os, re, time
from mainapp.controllers.project.fungi_genome_controller import FungiGenomeController
from mainapp.libs.signature import check_sig


class SequenceDownAction(FungiGenomeController):
    def __init__(self):
        super(SequenceDownAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'specimen_gene', 'type', 'submit_location', 'specimen_gene_new',
                        'client']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        project_sn = self.fungi_genome.get_projectsn(data.task_id)
        task_name = 'fungi.report.sequence_down'
        module_type = 'workflow'
        params_json = {
            'type': data.type,
            'specimen_gene': data.specimen_gene,
            'task_type': data.task_type,
            'submit_location': data.submit_location,
            'project_sn': project_sn,
            'task_id': data.task_id,
            'specimen_gene_new':data.specimen_gene_new,
        }
        main_table_name = 'SeqDown_' + data.type.capitalize() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        file_seqpath = self.fungi_genome.get_genefile_bysample(data.task_id, "KKKKKKKKKKKKKKKK", type="fnn",
                                                            predict=data.type)
        seq_path = (re.subn("\.fnn$", "", file_seqpath["KKKKKKKKKKKKKKKK"])[0]).split("KKKKKKKKKKKKKKKK")
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在下载'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_id = self.fungi_genome.insert_main_table('seq_down', mongo_data)
        options = {
            'seq_path': str(seq_path),
            'type': data.type,
            'main_id':str(main_id),
            'specimen_gene': str(data.specimen_gene),
            'client': data.client,
            'specimen_gene_new': str(data.specimen_gene_new),
            'task_id': data.task_id,
        }
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name='Seq_down/' + main_table_name,
                            task_id=data.task_id,
                            project_sn=project_sn, module_type=module_type, params=params_json)
        task_info = super(SequenceDownAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
