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


class DatasplitSampleMergeRenameAction(DatasplitController):
    """
    数据拆分样本fastq文件合并或重命名接口
    """
    def __init__(self):
        # super(DatasplitSampleMergeRenameAction, self).__init__(instant=True)
        super(DatasplitSampleMergeRenameAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["fx_id", "coll_id", "operation_type"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数: %s" % data.param}
                return json.dumps(info)
        if data.operation_type not in ["merge", "rename"]:
            info = {"success": False, "info": "operation_type:%s只能是merge/rename,请检查" % data.operation_type}
            return json.dumps(info)
        options = {
            "sample_list": data.coll_id,
            "coll_id": data.coll_id,
        }
        if data.operation_type == "merge":
            try:
                fx_id = int(data.fx_id) #fx_id在释放cleandata时为字符串，需要判定from Qinwen20211019
            except:
                fx_id = data.fx_id
            result = Datasplit("datasplit").coll_find_one("sg_split_specimen_merge", {"_id": ObjectId(data.coll_id), "fx_id": fx_id, "merge_st": True})
            if not result:
                info = {"success": False, "info": "没有在表sg_split_specimen_merge里找到_id: %s需要样本,请检查" % data.coll_id}
                return json.dumps(info)
            update_info = {str(data.coll_id): "sg_split_specimen_merge"}
            options["operation_type"] = "merge"
            options["update_info"] = json.dumps(update_info)
            to_file = ["datasplit_v2.export_merge_sample_info(sample_list)"]
            table_id = "SampleMerge_"+data.fx_id
            seq_number = "SampleMerge"
        elif data.operation_type == "rename":
            result = Datasplit("datasplit").coll_find_one("sg_split_specimen_rename", {"_id": ObjectId(data.coll_id), "fx_id": int(data.fx_id), "rename_st": True})
            if not result:
                info = {"success": False, "info": "没有在表sg_split_specimen_rename里找到_id: %s需要样本,请检查" % data.coll_id}
                return json.dumps(info)
            update_info = {str(data.coll_id): "sg_split_specimen_rename"}
            options["operation_type"] = "rename"
            options["update_info"] = json.dumps(update_info)
            to_file = ["datasplit_v2.export_rename_sample_info(sample_list)"]
            table_id = "SampleRename_"+data.fx_id
            seq_number = "SampleRename"
        self.set_sheet_data(name="datasplit_v2.sample_merge_rename", options=options, table_id=table_id, to_file=to_file, seq_number=seq_number, module_type="module")
        task_info = super(DatasplitSampleMergeRenameAction, self).POST()
        return json.dumps(task_info)
