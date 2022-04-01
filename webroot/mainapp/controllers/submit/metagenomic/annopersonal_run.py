# -*- coding: utf-8 -*-
import os
import json
import web
import json
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController


anno_col_info = {
    'go': ['anno_go', 'GO_Origin'],  # [col_name, table_name]
    'phi': ['anno_phi', 'PHI_Origin'],
    'mvirdb': ['anno_mvir', 'AnnoMvirdb_Origin'],
    'tcdb': ['anno_tcdb', 'AnnoTcdb_Origin'],
    'qs': ['anno_qs', 'QS_Origin'],
    'pfam': ['anno_pfam', 'PFAM_Origin'],
    'sec': ['anno_sec', 'AnnoSec_Origin'],
    'sec_de': ['anno_sec', 'AnnoSec_Origin_Deunclassified'],
    'sec_lca': ['anno_sec', 'AnnoSec_Origin_LCA'],
    't3ss': ['anno_ttss', 'TTSS_Origin'],
    't3ss_de': ['anno_ttss', 'TTSS_Origin_Deunclassified'],
    't3ss_lca': ['anno_ttss', 'TTSS_Origin_LCA'],
    'probio': ['anno_probio', 'Probio_Origin'],
    'probio_de': ['anno_probio', 'Probio_Origin_Deunclassified'],
    'probio_lca': ['anno_probio', 'Probio_Origin_LCA'],
    'p450': ['anno_cyps', 'CYPS_Origin'],
    'nr_de': ['anno_nr', 'NR_Origin_Deunclassified'],
    'nr_lca': ['anno_nr', 'NR_Origin_LCA'],
}


class AnnopersonalRunAction(MetagenomicController):
    def __init__(self):
        super(AnnopersonalRunAction, self).__init__(instant=False)
        self.my_db = self.metagenomic.db
        self.col = self.my_db['personal_anno']

    def check_run(self, anno):
        anno_info = self.col.find_one(
            {'task_id': self.task_id, 'type': anno}
        )
        if anno_info and anno_info['status'] in ['start', 'end']:
            return False, anno_info['in_overview']
        return True, 0

    def anno_post(self, anno, sheet_data={}):
        sheet_data['name'] = 'metagenomic.report.personal_anno'
        personal_mongo = dict(
            # type=anno, db=self.metagenomic.get_db_index(anno),
            type=anno, db='test',
            name=anno_col_info[anno][1], status='start',
            in_overview=0, desc='注释中', task_id=self.task_id,
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
        self.col.delete_one({'task_id': self.task_id, 'type': anno})
        personal_id = str(self.col.insert_one(personal_mongo).inserted_id)
        sheet_data['options'] = {
            'anno_type': anno,
            'main_col': anno_col_info[anno][0],
            'main_table_name': anno_col_info[anno][1],
            'task_id': self.task_id,
            'update_info': json.dumps({personal_id: 'personal_anno'}),
        }
        self.set_sheet_data(**sheet_data)
        self.col.update_one({'task_id': self.task_id, 'type': anno},
                            {'$set': {'run_id': self.workflow_id}})
        run_info = super(AnnopersonalRunAction, self).POST()
        if not run_info['success']:
            print(run_info)
            return self.sheet_data['id']
        # self.set_sheet_data(name=wf_name, options=options,
        #                     main_table_name= 'Personal_Anno',
        #                     module_type='workflow', project_sn=project_sn,
        #                     task_id=data.task_id, to_file=to_file))
        # pass

    def overview_post(self, running, sheet_data):
        sheet_data['name'] = 'metagenomic.report.personal_overview'
        personal_mongo = dict(
            type='overview', name='overview', status='start', desc='运行中',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            task_id=self.task_id, running_list=running
        )
        overview_info = self.col.find_one({'task_id': self.task_id,
                                           'type': 'overview'})
        if overview_info and overview_info['status'] == 'start':
            running += ',' + overview_info['running_list']
            update = {'running_list': ','.join(set(running.split(',')))}
            self.col.update_one({'task_id': self.task_id, 'type': 'overview'},
                                {'$set': update})
            return
        self.col.delete_one({'task_id': self.task_id, 'type': 'overview'})
        personal_id = str(self.col.insert_one(personal_mongo).inserted_id)
        update_info = json.dumps({personal_id: 'personal_anno'})
        sheet_data['options'] = {'task_id': self.task_id,
                                 'update_info': update_info
                                 }
        self.set_sheet_data(**sheet_data)
        self.col.update_one({'task_id': self.task_id, 'type': 'overview'},
                            {'$set': {'run_id': self.workflow_id}})
        run_info = super(AnnopersonalRunAction, self).POST()
        if not run_info['success']:
            print(run_info)
            return self.workflow_id

    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'anno_list']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, "info": "parameters " + argu + " is missing"}
                return json.dumps(info)
        self.task_id = data.task_id
        task_info = self.my_db['sg_task'].find_one({'task_id': self.task_id})
        print(self.task_id)
        print(task_info)
        sheet_data = {
            'main_table_name': 'Personal_Anno',
            'module_type': 'workflow', 'task_id': self.task_id,
            'project_sn': task_info['project_sn']
        }
        failed_anno = []
        running = []
        for anno in data.anno_list.split(','):
            run, in_overview = self.check_run(anno)
            if run:
                failed_id = self.anno_post(anno, sheet_data)
                if failed_id:
                    failed_anno.append([anno, failed_id])
                else:
                    running.append(anno)
            elif not in_overview:
                running.append(anno)
        post_info = self.overview_post(','.join(set(running)), sheet_data)
        if failed_anno:
            info = {
                "success": False,
                "info": "error in submit {}".format(failed_anno)
            }
        else:
            info = {
                "success": True,
                "info": "submit successed {}".format(running)
            }
        return json.dumps(info)
