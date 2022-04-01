# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON


class AaiAction(BacComparativeController):
    def __init__(self):
        super(AaiAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['task_id', 'evalue', 'submit_location', 'task_type', 'linkage', 'seq_dir', 'identity', 'aln_len', 'sample_list']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.genome_aai'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if float(data.identity) < 30 or float(data.identity) > 100:
            info = {'success': False, 'info': '%s参数缺少!' % str(data.identity)}
            return json.dumps(info)
        if float(data.evalue) <=0 or float(data.evalue) >=1:
            info = {'success': False, 'info': '%s参数缺少!' % data.evalue}
            return json.dumps(info)
        if data.linkage not in ['average', 'single', 'complete']:
            info = {'success': False, 'info': '%s参数缺少!' % data.linkage}
            return json.dumps(info)
        if float(data.aln_len) <=0 :
            info = {'success': False, 'info': '%s参数缺少!' % data.aln_len}
            return json.dumps(info)
        params = {
            'task_id': data.task_id,
            'evalue': data.evalue,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'linkage': data.linkage,
            'identity': data.identity,
            'sample_list': data.sample_list,
            'aln_len': data.aln_len
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'AAI_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.bac_comparative.insert_main_table('aai', mongo_data)
        update_info[str(main_id)] = 'aai'
        options = {
                   'evalue': data.evalue,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "seq_dir": data.seq_dir,
                   'linkage': data.linkage,
                   'identity': data.identity,
                   'aln_len': data.aln_len,
                   'sample_list':data.sample_list
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="AAI/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(AaiAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)