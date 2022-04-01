# -*- coding: utf-8 -*-
# __author__ = 'wangwenjie'

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.basic_molecular_controller import BasicMolecularController
# from mainapp.controllers.project.bac_identity_controller import BacidentityController


class LmplotAction(BasicMolecularController):
# class BacIdentityAction(BacidentityController):
    """
    标准曲线流程接口
    对接到isanger和拆分系统
    和基础分子MongoDB
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["tid"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        
        # target_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ob_storage"
        # params_json = {
        #     "download_file": data.download_id,
        #     "target_path": target_path
        # }
        # task_id = self.get_new_id(data.tid)
        result = Datasplit("basic_molecular").coll_find_one("sg_lmplot", {"task_id": data.tid})
        iid = result["_id"]
        update_info = {str(iid):"sg_lmplot"}
        options = {"update_info": json.dumps(update_info), "main_id": str(iid)}
        options["list_file"] = str(iid)
        to_file = ["basic_molecular.export_lmplot_params(list_file)"]
        self.set_sheet_data(name="basic_molecular.lmplot", options=options, table_id=data.tid,to_file=to_file, seq_number="lmplot", module_type="workflow")
        task_info = super(LmplotAction, self).POST()
        update_dict = {
                # "statue": split_status,
                "status": "start",
                "desc":"任务开始运行",
                # "cpc": cpc,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        Datasplit("basic_molecular").update_db_record("sg_lmplot", {"_id": ObjectId(iid)}, update_dict)
        print json.dumps(task_info)
        return json.dumps(task_info)
