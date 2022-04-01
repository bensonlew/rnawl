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


class DatasplitImportAction(DatasplitController):
    """
    线下数据导入对象存储接口
    功能：对数据进检查，将线下文件上传到对象存储，生成新的拆分任务
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["import_id"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数: %s" % data.param}
                return json.dumps(info)
        result = Datasplit("datasplit").coll_find_one("sg_import", {"_id": ObjectId(data.import_id)})
        if not result:
            info = {"success": False, "info": "没有在表sg_import里找到_id: %s，请检查" % data.import_id}
            return json.dumps(info)
        board_number = result["board_number"]
        task_sn = "import_" + result["task_sn"]
        if data.client == "client01":
            s3_path = "s3nb1://datasplit/" + datetime.datetime.now().strftime("%Y") + "/" + board_number + "/import_"
        else:
            s3_path = "s3://rerewrweset/files/datasplit/" + datetime.datetime.now().strftime("%Y") + "/" + board_number + "/import_"
        s3_path += str(data.import_id) + "_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        main_id = data.import_id
        Datasplit("datasplit").update_db_record("sg_import", {"_id": ObjectId(data.import_id)}, {"status": "start"})
        update_info = {str(main_id): "sg_import"}
        options = {
            "target_path": s3_path,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        self.set_sheet_data(name="datasplit_v2.datasplit_import", options=options, table_id=task_sn, module_type="tool")
        task_info = super(DatasplitImportAction, self).POST()
        return json.dumps(task_info)
