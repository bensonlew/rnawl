# -*- coding: utf-8 -*-
from __future__ import print_function
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from bson import ObjectId
from bson import SON
import datetime

class SsrAction(DenovoRnaV2Controller):
    def __init__(self):
        super(SsrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['rept_1', 'rept_2', 'rept_3', 'rept_4', 'rept_5', 'rept_6', 'ssr_distance','submit_location','task_id', 'task_type']
        for param in params_name:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!!" % param, "code": 'C1602001', "variables": var}
                return json.dumps(info)

        if int(data.rept_1) < 6:
            var = []
            var.append(data.rept_1)
            info = {"success": False, "info": "rept_1的值是%s，不在规定范围内!"%(data.rept_1), "code": 'C1602002', "variables": var}
            return json.dumps(info)
        if int(data.rept_2) < 3:
            var = []
            var.append(data.rept_2)
            info = {"success": False, "info": "rept_2的值是%s，不在规定范围内!" % (data.rept_2), "code": 'C1602003',
                    "variables": var}
            return json.dumps(info)
        if int(data.rept_3) < 3:
            var = []
            var.append(data.rept_3)
            info = {"success": False, "info": "rept_3的值是%s，不在规定范围内!" % (data.rept_3), "code": 'C1602004',
                    "variables": var}
            return json.dumps(info)
        if int(data.rept_4) < 3:
            var = []
            var.append(data.rept_4)
            info = {"success": False, "info": "rept_4的值是%s，不在规定范围内!" % (data.rept_4), "code": 'C1602005',
                    "variables": var}
            return json.dumps(info)
        if int(data.rept_5) < 3:
            var = []
            var.append(data.rept_5)
            info = {"success": False, "info": "rept_5的值是%s，不在规定范围内!" % (data.rept_5), "code": 'C1602006',
                    "variables": var}
            return json.dumps(info)
        if int(data.rept_6) < 3:
            var = []
            var.append(data.rept_6)
            info = {"success": False, "info": "rept_6的值是%s，不在规定范围内!" % (data.rept_6), "code": 'C1602007',
                    "variables": var}
            return json.dumps(info)
        if int(data.ssr_distance) < 0:
            var = []
            var.append(data.ssr_distance)
            info = {"success": False, "info": "ssr_distance的值是%s，不在规定范围内!" % (data.ssr_distance), "code": 'C1602008',
                    "variables": var}
            return json.dumps(info)

        task_name = 'denovo_rna_v2.report.ssr'
        # task_type = 'workflow'
        main_table_name = 'ssr_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        '''params_json里面的内容首先需要主表里面的运行的params，然后还有一些前端关注的参数，最后是比较排序后的字符串'''
        params_json = {
            'rept_1': int(data.rept_1),
            'rept_2': int(data.rept_2),
            'rept_3': int(data.rept_3),
            'rept_4': int(data.rept_4),
            'rept_5': int(data.rept_5),
            'rept_6': int(data.rept_6),
            'ssr_distance': int(data.ssr_distance),
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'task_id' : data.task_id
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

        '''mongo_data里面的内容就是和主表里面的内容保持一致'''
        sg_task_info = self.denovo_rna_v2.get_task_info(data.task_id)
        mongo_data = [
            ('project_sn', sg_task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', 'ssr_结果表'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('name', main_table_name),
        ]
        collection_name = 'sg_ssr'
        main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
        bed = self.denovo_rna_v2.get_bed_path(data.task_id)
        unigene_fa = self.denovo_rna_v2.get_unigenefa_path(data.task_id)
        update_info = {str(main_table_id): 'sg_ssr'}
        '''options里面的参数一部分是用于tool，一部分是用于to_file'''
        options = {
            'rept_1': int(data.rept_1),
            'rept_2': int(data.rept_2),
            'rept_3': int(data.rept_3),
            'rept_4': int(data.rept_4),
            'rept_5': int(data.rept_5),
            'rept_6': int(data.rept_6),
            'ssr_distance': int(data.ssr_distance),
            "update_info": json.dumps(update_info),
            "ssr_id": str(main_table_id),
            'bed': bed,
            'unigene_fa': unigene_fa
        }

        self.set_sheet_data(name=task_name, options=options, main_table_name= main_table_name,
                            task_id=data.task_id, project_sn=sg_task_info['project_sn'],
                             module_type='workflow', params=params, to_file=None)  # modified by hongdongxuan 20170322 在main_table_name前面加上文件输出的文件夹名

        task_info = super(SsrAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        else:
            pass
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)
