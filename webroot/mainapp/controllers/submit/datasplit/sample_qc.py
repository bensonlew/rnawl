# -*- coding: utf-8 -*-
# __author__ = "Xue Qinwen"
# last_modify: 20210914

import os
import web
import json
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class SampleQcAction(DatasplitController):
    """
    释放clean数据拆分的接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["qc_id"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存,请检查".format(name)}
                return json.dumps(info)
        result = Datasplit("datasplit").coll_find_one("sg_qc", {"_id": ObjectId(data.qc_id)})
        iid = result["_id"]
        task_sn = result["task_sn"]
        # split_status = result["split_status"]
        # split_type = result["split_type"]
        update_info = {str(data.qc_id): "sg_qc"}
        options = {"is_upload_cleandata":"true","update_info": json.dumps(update_info), "split_id": data.qc_id}
        task_name = "datasplit_v2.sample_qc"
        options["project_params"] = data.qc_id
        to_file = ["datasplit_v2.export_sample_qc_cleandata_params(project_params)"]
        update_dict = {
                # "statue": split_status,
                "status": "start",
                # "cpc": cpc,
                "desc": "开始进行拆分"
        }
        # elif data.data_source in ["specimen", "meta_raw"] or split_type == "second_split":
        
        Datasplit("datasplit").update_db_record("sg_qc", {"_id": ObjectId(data.qc_id)}, update_dict)
        self.set_sheet_data(name=task_name, options=options, table_id=task_sn, to_file=to_file, seq_number="sample_qc")
        task_info = super(SampleQcAction, self).POST()
        return json.dumps(task_info)

    # def new_split_path(self, split_path):
    #     if os.path.exists(split_path):
    #         return split_path
    #     if "ilustre" in split_path:
    #         split_path1 = split_path.replace("ilustre", "clustre")
    #         if os.path.exists(split_path1):
    #             return split_path1
    #     if "sglustre" in split_path:
    #         split_path1 = split_path.replace("sglustre", "ilustre")
    #         if os.path.exists(split_path1):
    #             return split_path1
    #     return split_path