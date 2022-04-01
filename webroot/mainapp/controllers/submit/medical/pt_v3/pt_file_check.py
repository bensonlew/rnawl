# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
# last modified: 2017.12.22
import web
import json
import os
import re
from mainapp.controllers.project.pt_controller import PtController
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
from bson import SON
import datetime


class PtFileCheckAction(PtController):
    """
    必需参数：'member_id', 'pepole_name', 'client', 'task_type', 'date'
    task_type参数控制接口调用方式： upload（用户上传）， post（前端传入json字符串）
    task_type为upload时，须提供up_path(用户上传文件的路径)、flowcell（版号）
    task_type为post时，须提供post_data(传入的json字符串)
    """
    def __init__(self):
        super(PtFileCheckAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['member_id', 'pepole_name', 'client', 'task_type', 'date']
        for names in params_name:
            if not (hasattr(data, names)):
                infor = {'success': False, 'info': "缺少参数{}".format(names)}
                return json.dumps(infor)
        if data.client not in ['client01', 'client03']:
            infor = {'success': False, 'info': "client{}类型不合法！只能是client01与client03".format(data.client)}
            return json.dumps(infor)
        # if str(data.client) == 'client01':
        #     base_path = '/mnt/ilustre/data/'
        # else:
        #     base_path = "/mnt/ilustre/tsanger-data/"
        indextype = 'single'  # single or double
        if hasattr(data, 'indextype'):
            indextype = data.indextype
        if data.task_type == "upload":
            params_name = ['up_path', 'flowcell']
            for names in params_name:
                if not (hasattr(data, names)):
                    infor = {'success': False, 'info': "缺少参数{}".format(names)}
                    return json.dumps(infor)

            flowcell = data.flowcell
            if os.path.exists("/mnt/clustre/upload/nextseq1/" + flowcell) or \
                    os.path.exists("/mnt/clustre/upload/nextseq/" + flowcell) or \
                    os.path.exists("/mnt/clustre/users/yixuezhuanhua/raw-data/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/upload/nextseq1/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/upload/nextseq/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/users/yixuezhuanhua/raw-data/" + flowcell) or \
                    re.match(r".*_wailaidata_.*", flowcell):
                pass
            else:
                infor = {'success': False, 'info': "板号{}名字命令错误，请实验人员重新命名上机"
                                                   "表名字！".format(flowcell)}
                return json.dumps(infor)

            infor_table = os.path.basename(data.up_path)
            insert_data = [
                ('created_ts', datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                ('status', "start"),
                ('desc', ""),
                ('pepole_name', data.pepole_name),
                ('member_id', data.member_id),
                ('info_table', infor_table)
            ]
            main_table_id = PT('pt_v3').add_main_id('sg_file_check')
            update_info = json.dumps({str(main_table_id): 'sg_file_check'})
            print data.up_path
            if data.up_path.startswith('/mnt') or re.match('.*://.*', str(data.up_path)):
                up_data_path = data.up_path
            else:
                up_data_path = '://'.join(data.up_path.split(':'))
            options = {
                'up': up_data_path,
                'date': data.date,
                'flowcell': flowcell,
                "update_info": update_info,
                "main_table_data": SON(insert_data),
                "main_id": str(main_table_id),  # 传递到tool中进行导表，代替tool中的update_info字段
                "task_type": "upload",
                'indextype': indextype
            }
            task_name = "medical.paternity_test_v3.report.file_check"
            task_type = "workflow"
            self.set_sheet_data_(name=task_name, options=options, module_type=task_type, db_type='pt_v3')
            task_info = super(PtFileCheckAction, self).POST()
            return json.dumps(task_info)
        elif data.task_type == "post":
            if not (hasattr(data, "json_data")):
                infor = {'success': False, 'info': "缺少参数json_data"}
                return json.dumps(infor)
            if not data.json_data:
                infor = {'success': False, 'info': "json_data为空"}
                return json.dumps(infor)
            flowcell = json.loads(data.json_data)[0]["batch_no"]
            if os.path.exists("/mnt/clustre/upload/nextseq1/" + flowcell) or \
                    os.path.exists("/mnt/clustre/upload/nextseq/" + flowcell) or \
                    os.path.exists("/mnt/clustre/users/yixuezhuanhua/raw-data/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/upload/nextseq1/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/upload/nextseq/" + flowcell) or \
                    os.path.exists("/mnt/ilustre/users/yixuezhuanhua/raw-data/" + flowcell) or\
                    re.match(r".*_wailaidata_.*", flowcell):
                pass
            else:
                infor = {'success': False, 'info': "板号{}名字命令错误，"
                                                   "请实验人员重新命名上机表名字！".format(flowcell)}
                return json.dumps(infor)

            insert_data = [
                ('created_ts', datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                ('status', "start"),
                ('desc', ""),
                ('pepole_name', "系统录入"),
                ('member_id', data.member_id),
                ('info_table', flowcell + ".xlsx")
            ]
            main_table_id = PT('pt_v3').add_main_id('sg_file_check')
            update_info = json.dumps({str(main_table_id): 'sg_file_check'})
            options = {
                'json_data': data.json_data,
                'date': data.date,
                'flowcell': flowcell,
                "update_info": update_info,
                "main_table_data": SON(insert_data),
                "main_id": str(main_table_id),  # 传递到tool中进行导表，代替tool中的update_info字段
                "task_type": "post",
                'indextype': indextype
            }
            task_name = "medical.paternity_test_v3.report.file_check"
            task_type = "workflow"
            self.set_sheet_data_(name=task_name, options=options, module_type=task_type, db_type='pt_v3')
            task_info = super(PtFileCheckAction, self).POST()
            return json.dumps(task_info)
        else:
            infor = {'success': False, 'info': "task_type error(upload|post), '{}' got".format(data.task_type)}
            return json.dumps(infor)
