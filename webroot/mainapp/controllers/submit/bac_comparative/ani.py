# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON


class AniAction(BacComparativeController):
    def __init__(self):
        super(AniAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'method', 'submit_location', 'task_type', 'linkage', 'seq_dir', 'sample_list']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.genome_ani'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if data.method not in ['ANIm', 'ANIb', 'ANIblastall', 'TETRA']:
            info = {'success': False, 'info': '%s参数缺少!' % data.method}
            return json.dumps(info)
        if data.linkage not in ['average', 'single', 'complete']:
            info = {'success': False, 'info': '%s参数缺少!' % data.linkage}
            return json.dumps(info)
        params = {
            'task_id': data.task_id,
            'method': data.method,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'linkage': data.linkage,
            'sample_list': data.sample_list
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'ANI_' + data.method + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_id = self.bac_comparative.insert_main_table('ani', mongo_data)
        update_info[str(main_id)] = 'ani'
        options = {
                   'method': data.method,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "seq_dir": data.seq_dir,
                   'linkage': data.linkage,
                   'sample_list':data.sample_list
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="ANI/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(AniAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)