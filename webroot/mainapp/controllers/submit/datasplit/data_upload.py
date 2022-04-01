# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190327

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class DataUploadAction(DatasplitController):
    """
    数据上传接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["fx_sn", "fx_id", "fx_task_id", "fx_task_result_id", "fx_task_result_server_path"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数: %s" % data.param}
                return json.dumps(info)
        if data.client == "client01":
            s3_path = "s3://datasplit/" + datetime.datetime.now().strftime("%Y") + "/" + str(data.fx_sn) + "/" + str(data.fx_id) + "/" + str(data.fx_task_id)
        else:
            s3_path = "s3://rerewrweset/files/datasplit/" + datetime.datetime.now().strftime("%Y") + "/" + str(data.fx_sn) + "/" + str(data.fx_id) + "/" + str(data.fx_task_id)
        s3_path += str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        result = Datasplit("datasplit").coll_find_one("sg_data_upload", {"s3_path": s3_path, "status": "end"})
        if result:
            info = {"success": False, "info": "对象存储路径:%s已经存在，请检查" % s3_path}
            return json.dumps(info)
        result = Datasplit("datasplit").coll_find_one("sg_data_upload", {"s3_path": s3_path, "status": "start"})
        if result:
            info = {"success": False, "info": "对象存储路径:%s正在上传，请稍等" % s3_path}
            return json.dumps(info)
        if not data.fx_task_result_server_path.startswith("/"):
            info = {"success": False, "info": "fx_task_result_server_path:%s要是服务器路径,请检查" % data.fx_task_result_server_path}
            return json.dumps(info)
        mongo_data = [
            ("status", "start"),
            ("fx_sn", data.fx_sn),
            ("fx_id", data.fx_id),
            ("fx_task_id", data.fx_task_id),
            ("fx_task_result_id", data.fx_task_result_id),
            ("fx_task_result_server_path", data.fx_task_result_server_path),
            ("s3_path", s3_path),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Datasplit("datasplit").insert_main_table(collection="sg_data_upload", data=mongo_data)
        update_info = {str(main_id): "sg_data_upload"}
        options = {
            "upload_path": data.fx_task_result_server_path,
            "target_path": s3_path,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        task_sn = "dataupload_" + str(data.fx_sn)
        self.set_sheet_data(name="datasplit_v2.data_upload", options=options, table_id=task_sn, module_type="tool")
        task_info = super(DataUploadAction, self).POST()
        return json.dumps(task_info)
