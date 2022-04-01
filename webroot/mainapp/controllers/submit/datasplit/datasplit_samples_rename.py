# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20201105

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class DatasplitSamplesRenameAction(DatasplitController):
    """
    数据拆分样本样本重命名接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["fx_id"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数: %s" % data.param}
                return json.dumps(info)
        result = Datasplit("datasplit").coll_find("sg_split_specimen_rename", {"fx_id": int(data.fx_id), "rename_st": True})
        if not result:
            info = {"success": False, "info": "没有在表sg_split_specimen_rename里找到fx_id: %s有需要重命名的样本,请检查" % data.fx_id}
            return json.dumps(info)
        options = {
            "sample_list": data.fx_id,
            "fx_id": data.fx_id,
            "operation_type": "rename"
        }
        to_file = ["datasplit_v2.export_rename_samples_info(sample_list)"]
        self.set_sheet_data(name="datasplit_v2.samples_merge_rename", options=options, table_id=data.fx_id, to_file=to_file, module_type="workflow")
        task_info = super(DatasplitSamplesRenameAction, self).POST()
        return json.dumps(task_info)
