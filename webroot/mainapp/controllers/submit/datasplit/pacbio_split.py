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


class PacbioSplitAction(DatasplitController):
    """
    三代数据拆分的接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["split_id"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存,请检查".format(name)}
                return json.dumps(info)
        result = Datasplit("datasplit").coll_find_one("sg_pacbio", {"_id": ObjectId(data.split_id)})
        iid = result["_id"]
        task_sn = result["task_sn"]
        # split_status = result["split_status"]
        # split_type = result["split_type"]
        update_info = {str(data.split_id): "sg_pacbio"}
        options = {"update_info": json.dumps(update_info), "split_id": data.split_id}
        # params_json = json.loads(result["params"])
        # if data.data_source == "library" or split_type == "first_split":
        # if data.data_source == "library":
        task_name = "datasplit_v2.pacbio_split"
            # params_json = json.loads(result["params"])
            # for lane in params_json["library_split"]:
        split_path = self.new_split_path(result["bam_path"])
        if not os.path.exists(split_path):
            split_path = split_path.replace("ilustre", "clustre")
            print split_path
            if not os.path.exists(split_path):
                update_dict = {"desc": "下机数据路径:%s没有找到,未下机或者下机路径错误"  % split_path}
                Datasplit("datasplit").update_db_record("sg_split", {"_id": ObjectId(data.split_id)}, update_dict)
                info = {"success": False, "info": "下机数据路径:%s 没有找到，请检查" % split_path}
                return json.dumps(info)
        options["pacbio_params"] = data.split_id
            # options["run_type"] = "auto" if split_status["second_split"] == "no" else "manual"
        to_file = ["datasplit_v2.export_pacbio_split(pacbio_params)"]
            # split_status["first_split"] = "start"
            # split_status["cpc"] = "--"
            # cpc = "no"
        # for i in params_json.keys():
        #         if i in ["dna", "rna", "prokaryotic_rna", "lncrna", "microbial_genome"]:
        #             split_status["cpc"] = "no"
        #             cpc = "wait"
        update_dict = {
                # "statue": split_status,
                "status": "start",
                # "cpc": cpc,
                "desc": "开始进行三代拆分"
        }
        # elif data.data_source in ["specimen", "meta_raw"] or split_type == "second_split":
        
        Datasplit("datasplit").update_db_record("sg_pacbio", {"_id": ObjectId(data.split_id)}, update_dict)
        self.set_sheet_data(name=task_name, options=options, table_id=task_sn, to_file=to_file, seq_number="pacbio_split")
        task_info = super(PacbioSplitAction, self).POST()
        return json.dumps(task_info)

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path