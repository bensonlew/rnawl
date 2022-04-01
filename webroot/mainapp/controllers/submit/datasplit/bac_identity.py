# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController
# from mainapp.controllers.project.bac_identity_controller import BacidentityController


class BacIdentityAction(DatasplitController):
# class BacIdentityAction(BacidentityController):
    """
    菌鉴流程接口
    对接到tsanger和拆分系统DatasplitController
    和拆分MongoDB
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["tid","list_url","sample_info_url"]
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
        
        mongo_data = [
            ("task_id", data.tid),
            ("status", "start"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Datasplit("datasplit").insert_main_table(collection="bac_identity", data=mongo_data)
        # to_file = ["datasplit.excport_path(download_file)"]
        update_info = {str(main_id):"bac_identity"}
        options = {
            "list_url": data.list_url,
            "sample_info_url": data.sample_info_url,
            "main_id": str(main_id),
            "update_info":json.dumps(update_info)
        }
        self.set_sheet_data(name="tool_lab.bac_identity", options=options, table_id=data.tid, seq_number="bac_identity", module_type="workflow")
        task_info = super(BacIdentityAction, self).POST()
        print json.dumps(task_info)
        return json.dumps(task_info)
