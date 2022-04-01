# -*- coding: utf-8 -*-
# __author__: zengjing
# last_modify: 20210621

import os
import re
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class DatasplitUploadAction(DatasplitController):
    """
    三代、蛋白、代谢原始文件上传
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["project_type"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        if data.project_type not in ["protein", "metabolome", "pacbio"]:
            info = {"success": False, "info": "参数project_type{}只能是protein/metabolome/pacbio".format(data.project_type)}
            return json.dumps(info)
        if data.project_type == "pacbio":
            params = ["lane_library_ids"]
            for name in params:
                if not hasattr(data, name):
                    info = {"success": False, "info": "参数{}不存在".format(name)}
                    return json.dumps(info)
        else:
            params = ["contract_sn", "task_sn", "sample_name"]
            for name in params:
                if not hasattr(data, name):
                    info = {"success": False, "info": "参数{}不存在".format(name)}
                    return json.dumps(info)
        seq_number = "upload_" + data.project_type
        if data.project_type == "pacbio":
            # lane_library_ids = json.loads(data.lane_library_ids)
            lane_library_ids = data.lane_library_ids.split(",")
            options = {
                "project_table": ";".join(lane_library_ids),
                "lane_library_ids": ";".join(lane_library_ids),
                "is_upload": True,
            }
            to_file = ["datasplit_v2.export_pacbio_upload_info(project_table)"]
            task_id = "upload_pacbio_"+lane_library_ids[0]
            self.set_sheet_data(name="datasplit_v2.upload_pacbio", options=options, table_id=task_id, to_file=to_file, seq_number=seq_number, module_type="tool")
        else:
            # sample_name = json.loads(data.sample_name)
            sample_name = data.sample_name.split(",")
            options = {
                "project_sn": data.contract_sn,
                "task_sn": data.task_sn,
                "samples": ";".join(sample_name),
                "is_upload": True,
            }
            task_id = "upload_protein_"+data.contract_sn
            self.set_sheet_data(name="datasplit_v2.upload_protein_metabolome", options=options, table_id=task_id, to_file="", seq_number=seq_number, module_type="tool")
        task_info = super(DatasplitUploadAction, self).POST()
        return json.dumps(task_info)
