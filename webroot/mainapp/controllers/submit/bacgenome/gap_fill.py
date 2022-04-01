# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import web, re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON
import os


class GapFillAction(BacgenomeController):
    def __init__(self):
        super(GapFillAction, self).__init__(instant=False)
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
        print (data)
        default_argu = ['task_id', 'sample', 'gap', 'gap_len','scaffold', 'overlap', 'identity', 'submit_location', 'task_type', 'select_seq', 'custom_scaffold']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)

        task_name = 'bacgenome.report.gap_fill'
        module_type = 'workflow'
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        fasta = ''
        fasta2 = ''
        if data.select_seq in ["assemble"]:
            if data.gap not in ["unicycler_auto"]:
                fasta = self.bacgenome.get_gap_fill(data.task_id, data.sample, data.gap)
                fasta2 = self.bacgenome.get_gap_fill(data.task_id, data.sample, data.gap2)
            else:
                fasta = self.bacgenome.get_gap_fill(data.task_id, data.sample, "unicycler")
        elif data.select_seq in ["gap_fill"]:
            fasta = self.bacgenome.get_gap_fill2(ObjectId(data.gap))
            fasta2 = self.bacgenome.get_gap_fill(data.task_id, data.sample, data.gap2)
        params = {
            'task_id': data.task_id,
            'sample': data.sample,
            'select_seq': data.select_seq,
            'gap_len': data.gap_len,
            'gap': data.gap,
            'overlap': data.overlap,
            'identity': data.identity,
            'submit_location': data.submit_location,
            'custom_scaffold':data.custom_scaffold,
            'task_type': int(data.task_type)
        }
        if data.select_seq in ["assemble"]:
            if data.gap not in ["unicycler_auto"]:
                params['gap2'] =data.gap2
                params['gap2_len'] = data.gap2_len
                params['scaffold'] = json.dumps(eval(data.scaffold),sort_keys=True, separators=(',',':'))
            else:
                params['scaffold'] =json.dumps({"Chromosome1": ["1"]},sort_keys=True, separators=(',',':'))
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        main_table_name = data.sample + "_GapFill_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ("sample",data.sample)
        ]
        main_id = self.bacgenome.insert_main_table('gap_fill', mongo_data)
        self.bacgenome.insert_main_table_new("gap_fill", str(main_id), {"main_id": main_id})
        update_info[str(main_id)] = 'gap_fill'
        options = {
                    'samp': data.sample,
                    'update_info': json.dumps(update_info),
                    'main_id': str(main_id),
                    'task_id': data.task_id,
                    "chr_fa": fasta,
                    "chr_len":data.gap_len,
                    "genome_json": data.scaffold,
                    "overlap": int(data.overlap),
                    "select_seq": data.select_seq,
                    "identity": float(data.identity)
                   }
        if data.select_seq in ["assemble"]:
            if data.gap not in ["unicycler_auto"]:
                options['cir_len'] = data.gap2_len
                options['cir_fa'] = fasta2
                options['gap'] = data.gap
        else:
            options['cir_len'] = data.gap2_len
            options['cir_fa'] = fasta2
        if data.select_seq in ["gap_fill"]:
            gap_id = ObjectId(data.gap)
            options["pre_log"] = gap_id
            to_file = ['bac_assem.export_log_by_gapid(pre_log)']
            self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name="GapFill/" + main_table_name,
                                module_type=module_type,
                                project_sn=project_sn,
                                task_id=data.task_id,
                                params=params,
                                to_file=to_file)
        else:
            if data.gap in ["unicycler"]:
                options["pre_log"] = data.gap
                to_file = ['bac_assem.export_log_by_sof(pre_log)']
                self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name="GapFill/" + main_table_name,
                                module_type=module_type,
                                project_sn=project_sn,
                                task_id=data.task_id,
                                params=params,
                                to_file=to_file)
            else:
                self.set_sheet_data(name=task_name,
                                    options=options,
                                    main_table_name="GapFill/" + main_table_name,
                                    module_type=module_type,
                                    project_sn=project_sn,
                                    task_id=data.task_id,
                                    params=params)
        task_info = super(GapFillAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)
