# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20201013
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bac_comparative import BacComparative
from bson import SON


class PanGenomeAction(BacComparativeController):
    def __init__(self):
        super(PanGenomeAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'pan_id','submit_location', 'task_type', 'cluster_path', "gene_dir"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.pan_genome'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)

        params = {
            'pan_id': data.pan_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            # 'cluster_path': data.cluster_path,
            # 'gene_dir': data.gene_dir,
            'task_id': data.task_id
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = "Pan_Formula_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ("cluster_path", data.cluster_path),
            ("gene_dir", data.gene_dir),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bac_comparative.insert_main_table('pan_formula', mongo_data)
        update_info[str(main_id)] = 'pan_formula'
        options = {'from_main_id': data.pan_id,
                   'update_info': json.dumps(update_info),
                   'main_id': str(main_id),
                   "gene_dir": data.gene_dir,
                   'cluster_path': data.cluster_path,
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Pan_Formula_/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(PanGenomeAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)