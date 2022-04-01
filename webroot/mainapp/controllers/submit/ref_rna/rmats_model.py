# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 09:34

import re, os, Bio, argparse, sys, fileinput

'''
跑事件模式图所用接口
'''
import json
import random
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from biocluster.config import Config
import types
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
import web
from mainapp.models.mongo.ref_rna import RefRna
import subprocess


class RmatsModelAction(RefRnaController):
    def __init__(self):
        super(RmatsModelAction, self).__init__(instant=False)
    
    def GET(self):
        return 'khl'

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        splicing_rmats_info = self.ref_rna.get_main_info(data.splicing_id, 'sg_splicing_rmats')
        if not splicing_rmats_info:
            info = {"success": False, "info": "splicing_rmats主表信息不存在，请检查参数是否正确！"}
            return json.dumps(info)
        try:
            rmats_out_root_dir = splicing_rmats_info['rmats_out_root_dir']
        except Exception as e:
            info = {"success": False, "info": "splicing_rmats主表信息中不存在结果目录信息！"}
            return json.dumps(info)
        task_info = self.ref_rna.get_task_info(splicing_rmats_info['task_id'])
        
        return_result = self.check_options(data, rmats_out_root_dir)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_type = data.task_type
        task_name = 'ref_rna.report.rmats_model'
        my_param = {'splicing_id': data.splicing_id, 'event_id': data.event_id, 'event_type': data.event_type,
                    'task_type': task_type, 'submit_location': data.submit_location}
        # rmats_out_root_dir = RefRna().get_main_info(ObjectId(data.splicing_id), 'sg_splicing_rmats')[
        #     'rmats_out_root_dir']
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        
        main_table_name = "SplicingRmatsModel_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        task_id = splicing_rmats_info['task_id']
        project_sn = splicing_rmats_info['project_sn']
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':'))),
            ('splicing_id', ObjectId(data.splicing_id))
        ]
        collection_name = "sg_splicing_rmats_model"
        main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}
        update_info = json.dumps(update_info)
        options = {
            "splicing_id": data.splicing_id,
            "update_info": update_info,
            "event_id": data.event_id,
            "event_type": data.event_type,
            "rmats_model_id": str(main_table_id),
            "rmats_out_root_dir": rmats_out_root_dir,
            
        }
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, task_id=task_id,
                            project_sn=project_sn, module_type='workflow')
        
        task_info = super(RmatsModelAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
    
    def check_options(self, data, out_root):
        params_name = ["event_type", "event_id", "splicing_id","task_type"]
        success = []
        for name in params_name:
            if not hasattr(data, name):
                success.append("缺少参数：%s" % name)
        event_type = data.event_type
        if event_type not in ["A3SS", "A5SS", "MXE", "SE", "RI"]:
            success.append("%s分析不存在" % event_type)
        rmats_event_file = os.path.join(out_root,
                                        "MATS_output/" + data.event_type + '.MATS.ReadsOnTargetAndJunctionCounts.alter_id.txt')
        numbers = [e.strip() for e in
                   subprocess.check_output("awk -F '\\t' 'NR>=2{print $1}' %s" % rmats_event_file, shell=True).strip().split("\n")]
        event_id = data.event_id
        if event_id.strip() not in numbers:
            success.append("%s事件的差异记录不存在" % event_id)
        return success
