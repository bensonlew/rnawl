# -*- coding: utf-8 -*-
# __author__: zengjing
# last_modify: 20210223

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


class MetabolomePredealAction(DatasplitController):
    """
    代谢预处理
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        # lib_search_dir:搜库文件夹 qcs:保留的qc样本 filter_free:目标代谢物
        # samples:样本名称 olds:需修改的样本名称 news:修改后的样本名称 groups:分组名 controls:对照组 others:实验组
        params = ["project_sn", "fx_sn", "fx_task_id", "lib_search_dir", "qcs", "filter_free", "samples", "news",\
            "groups", "controls", "others","sure_score","kegg_amino_acid"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        lib_search_dir = data.lib_search_dir
        result = Datasplit("datasplit").coll_find_one("sg_specimen_s3_protein", {"search_lib_path": data.lib_search_dir})
        if result:
            if os.path.exists(result["search_lib_ts_path"]):
                lib_search_dir = result["search_lib_ts_path"]
        mongo_data = [
            ("status", "start"),
            ("project_sn", data.project_sn),
            ("fx_sn", data.fx_sn),
            ("fx_task_id", data.fx_task_id),
            ("lib_search_dir", data.lib_search_dir),
            ("qcs", data.qcs),
            ("filter_free", data.filter_free),
            ("samples", data.samples),
            # ("olds", data.olds),
            ("news", data.news),
            ("groups", data.groups),
            ("controls", data.controls),
            ("others", data.others),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Datasplit("datasplit").insert_main_table(collection="sg_predeal", data=mongo_data)
        update_info = {str(main_id): "sg_predeal"}
        options = {
            "lib_search_dir": data.lib_search_dir,
            "qcs": data.qcs,
            "filter_free": data.filter_free,
            "samples": data.samples,
            "news": data.news,
            "groups": data.groups,
            "controls": data.controls,
            "others": data.others,
            "sure_score":data.sure_score,
            "kegg_amino_acid":data.kegg_amino_acid,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        task_sn = "predeal_" + str(data.fx_sn)
        self.set_sheet_data(name="datasplit_v2.metabolome_predeal", options=options, table_id=task_sn, seq_number="metabolome_predeal", module_type="tool")
        task_info = super(MetabolomePredealAction, self).POST()
        return json.dumps(task_info)
