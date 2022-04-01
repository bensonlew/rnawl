# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
# last modified: 2017.12.22
import web
import json
import os
from mainapp.controllers.project.pt_controller import PtController
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
from bson import SON
import datetime


class PtFileCheckAction(PtController):
    """
    运行逻辑：前端传入参数包含up_path;down_path（和前端确认了dwon+path必须传，如果没有文件夹的时候为''）;batch_path;
    date(前端输入的是当天的时间加上随机数);flowcell（可无）; pepole_name;member_id
    """
    def __init__(self):
        super(PtFileCheckAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['up_path', 'down_path', 'batch_path', 'date', 'member_id', 'pepole_name', 'client']
        for names in params_name:
            if not (hasattr(data, names)):
                infor = {'success': False, 'info': "缺少参数{}".format(names)}
                return json.dumps(infor)
        if data.client not in ['client01', 'client03']:
            infor = {'success': False, 'info': "client{}类型不合法！只能是client01与client03".format(data.client)}
            return json.dumps(infor)
        if str(data.client) == 'client01':
            base_path = '/mnt/ilustre/data/'
        else:
            base_path = "/mnt/ilustre/tsanger-data/"
        flowcell = '_'.join(os.path.basename(str(data.up_path)).strip().split('_')[:-1])
        print flowcell
        if os.path.exists("/mnt/clustre/upload/nextseq1/" + flowcell) or os.path.exists("/mnt/clustre/upl"
                                                                                        "oad/nextseq/" + flowcell):
            pass
        else:
            infor = {'success': False, 'info': "板号{}名字命令错误，请实验人员重新命名上机表名字！".format(flowcell)}
            return json.dumps(infor)
        infor1 = "实验批次表:" + os.path.basename(data.batch_path)
        infor2 = "线上上机表(LIMS系统生成):" + os.path.basename(data.up_path)
        if str(data.down_path):
            infor3 = "线下上机表(返工/测试样品):" + os.path.basename(data.down_path)
        else:
            infor3 = "线下上机表(返工/测试样品): --"
        infor_table = [infor1, infor2, infor3]
        infor_table = "\n".join(infor_table)
        insert_data = [
            ('created_ts', datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            ('status', "start"),
            ('desc', ""),
            ('pepole_name', data.pepole_name),
            ('member_id', data.member_id),
            ('info_table', infor_table)
        ]
        main_table_id = PT('pt_v2').add_main_id('sg_file_check')
        update_info = json.dumps({str(main_table_id): 'sg_file_check'})
        options = {
            'up': os.path.join(base_path, data.up_path),
            'date': data.date,
            'batch': os.path.join(base_path, data.batch_path),
            'flowcell': flowcell,
            "update_info": update_info,
            "main_table_data": SON(insert_data),
            "main_id": str(main_table_id)  # 传递到tool中进行导表，代替tool中的update_info字段
        }
        if hasattr(data, 'down_path') and str(data.down_path):
            options.update({'down': os.path.join(base_path, data.down_path)})
        print options
        task_name = "medical.paternity_test_v2.report.file_check"
        task_type = "workflow"
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, db_type='pt_v2')
        task_info = super(PtFileCheckAction, self).POST()
        return json.dumps(task_info)

