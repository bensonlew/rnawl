# -*- coding: utf-8 -*-
# __author__ = 'ysh'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class BlastAction(BacgenomeController):
    def __init__(self):
        super(BlastAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:",data
        default_argu = ['task_id','specimen', 'submit_location', 'task_type','method','evalue','max_target_num','w_size','sequence']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.blast'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        nr_result = self.bacgenome.get_one_common_info('anno_nr', {'task_id':data.task_id,"status" : "end"})
        nr_id = str(nr_result['main_id'])
        if data.method == "blastn":
            if not data.w_size in ["16","20","24","28","32","48","64","128","256"]:
                info = {'success': False, 'info': 'w_size must be 16,20,24,28,32,48,64,128,256 when blastn !'}
                return json.dumps(info)
        if data.method in ["blastp","blastx","tblastn"]:
            if not data.w_size in ["2","3","6"]:
                info = {'success': False, 'info': 'w_size must be 2,3,6 when blastp,blastx,tblastn!'}
                return json.dumps(info)
        if data.method == "tblastx":
            if not data.w_size in ["2","3"]:
                info = {'success': False, 'info': 'w_size must be 2,3 when tblastx!'}
                return json.dumps(info)
        seq = data.sequence.replace('&gt;','>')
        params = {
            'specimen': data.specimen,
            'method': data.method,
            'evalue': data.evalue,
            'max_target_num': int(data.max_target_num),
            'sequence': seq,
            'w_size': int(data.w_size),
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        update_info = {}
        main_table_name = 'blast_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_table_id = self.bacgenome.insert_main_table('blast', mongo_data)
        update_info[str(main_table_id)] = 'blast'
        options = {
                    'specimens': data.specimen,
                    'method': data.method,
                    'evalue': data.evalue,
                    'max_target_num': data.max_target_num,
                    'w_size': data.w_size,
                    'update_info': json.dumps(update_info),
                    'main_id': str(main_table_id),
                    'task_id': data.task_id,
                    'nr_id' : nr_id,
                    'sequence': data.sequence,
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name='GeneBlast/' + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(BlastAction, self).POST()
        # return spe_list
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
