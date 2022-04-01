# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

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


# class BacIdentityAction(DatasplitController):
class BacidentityUploadAction(DatasplitController):
    """
    菌鉴流程接口,测试机投递任务接口，
    对接到BacidentityController和小工具测试tool_lab mongo库
    """
    def __init__(self):
        super(BacidentityUploadAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["bac_id","tid","list_url","sample_info_url","fasta_dir","new_pictrue_name"]
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
        
        # mongo_data = [
        #     ("_id", data.tid),
        #     ("status", "start"),
        #     ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        # ]
        # main_id = Datasplit("tool_lab").insert_main_table(collection="bac_identity", data=mongo_data).
        update_dict ={
            "desc":"菌鉴分析中",
            "work_env":"isanger"
        }
        Datasplit("datasplit").update_db_record("bac_identity",{"_id":ObjectId(data.bac_id)},update_dict)
        # to_file = ["datasplit.excport_path(download_file)"]
        update_info = {data.bac_id:"bac_identity"}
        options = {
            "list_url": data.list_url,
            "sample_info_url": data.sample_info_url,
            "main_id": data.bac_id,
            "update_info":json.dumps(update_info),
            "fasta_dir": data.fasta_dir,
            "new_pictrue_name":data.new_pictrue_name
        }
        self.set_sheet_data(name="datasplit_v2.bac_identity", options=options, table_id=data.tid, seq_number="bac_identity", module_type="workflow")
        task_info = super(BacidentityUploadAction, self).POST()
        print json.dumps(task_info)
        return json.dumps(task_info)
