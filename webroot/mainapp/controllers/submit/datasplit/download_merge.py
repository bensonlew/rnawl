# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

import os
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class DownloadMergeAction(DatasplitController):
    """
    数据拆分接口，进行数据的下载和合并
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["download_id"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        target_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ob_storage"
        params_json = {
            "download_file": data.download_id,
            "target_path": target_path
        }
        options = {
            "download_file": data.download_id,
            "target_path": target_path
        }
        to_file = ["datasplit.excport_path(download_file)"]
        self.set_sheet_data(name="datasplit.download", options=options, table_id=data.download_id, to_file=to_file, seq_number=None, module_type="tool")
        task_info = super(DownloadMergeAction, self).POST()
        return json.dumps(task_info)
