# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190311

import os
import web
import json
import time
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class DatasplitV2Action(DatasplitController):
    """
    数据拆分,一次拆分及自动拆分的接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["split_id", "data_source"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存,请检查".format(name)}
                return json.dumps(info)
        result = Datasplit("datasplit").coll_find_one("sg_split", {"_id": ObjectId(data.split_id)})
        task_sn = result["task_sn"]
        seq_number = result["board_number"]
        split_status = result["split_status"]
        split_type = result["split_type"]
        update_info = {str(data.split_id): "sg_split"}
        options = {"update_info": json.dumps(update_info), "split_id": data.split_id}
        params_json = json.loads(result["params"])
        # if data.data_source == "library" or split_type == "first_split":
        if data.data_source == "library":
            task_name = "datasplit_v2.library_split"
            # params_json = json.loads(result["params"])
            for lane in params_json["library_split"]:
                split_path = params_json["library_split"][lane]["split_path"]
                # if not os.path.exists(split_path):
                #     split_path = split_path.replace("ilustre", "clustre")
                print split_path
                split_path = self.new_split_path(split_path)
                if not os.path.exists(split_path):
                    update_dict = {"desc": "下机数据路径:%s没有找到,未下机或者下机路径错误"  % split_path}
                    Datasplit("datasplit").update_db_record("sg_split", {"_id": ObjectId(data.split_id)}, update_dict)
                    info = {"success": False, "info": "下机数据路径:%s 没有找到，请检查" % split_path}
                    return json.dumps(info)
            options["library_params"] = data.split_id
            options["run_type"] = "auto" if split_status["second_split"] == "no" else "manual"
            to_file = ["datasplit_v2.export_library_params(library_params)"]
            split_status["first_split"] = "start"
            split_status["cpc"] = "--"
            cpc = "no"
            for i in params_json.keys():
                if i in ["dna", "rna", "prokaryotic_rna", "lncrna", "microbial_genome"]:
                    split_status["cpc"] = "no"
                    cpc = "wait"
            update_dict = {
                "split_status": split_status,
                "status": "start",
                "cpc": cpc,
                "desc": "开始进行文库拆分"
            }
        # elif data.data_source in ["specimen", "meta_raw"] or split_type == "second_split":
        elif data.data_source in ["specimen", "meta_raw"]:
            task_name = "datasplit_v2.sample_split"
            options["project_params"] = data.split_id
            options["run_type"] = "auto" if split_status["qc"] == "no" else "manual"
            to_file = ["datasplit_v2.export_sample_split_params(project_params)"]
            split_status["second_split"] = "start"
            split_status["cpc"] = "--"
            cpc = "no"
            for i in params_json.keys():
                if i in ["dna", "rna", "prokaryotic_rna", "lncrna", "microbial_genome"]:
                    split_status["cpc"] = "no"
                    cpc = "wait"
            update_dict = {
                "split_status": split_status,
                "status": "start",
                "cpc": cpc,
                "desc": "开始进行样本拆分"
            }
        # elif data.data_source in ["other_raw", "qc"] or split_type == "qc":
        elif data.data_source in ["other_raw", "qc"]:
            task_name = "datasplit_v2.sample_qc"
            options["project_params"] = data.split_id
            to_file = ["datasplit_v2.export_sample_qc_params(project_params)"]
            split_status["qc"] = "start"
            split_status["cpc"] = "--"
            cpc = "no"
            for i in params_json.keys():
                if i in ["dna", "rna", "prokaryotic_rna", "lncrna", "microbial_genome"]:
                    split_status["cpc"] = "no"
                    cpc = "wait"
            update_dict = {
                "split_status": split_status,
                "status": "start",
                "cpc": cpc,
                "desc": "开始进行样本质控"
            }
        elif data.data_source == "cpc":
            task_name = "datasplit_v2.sample_cpc"
            options["sample_list"] = data.split_id
            to_file = ["datasplit_v2.export_sample_cpc_params(sample_list)"]
            split_status["cpc"] = "start"
            update_dict = {
                "split_status": split_status,
                "cpc": "start",
                "desc": "开始进行cpc"
            }
        else:
            info = {"success": False, "info": "参数数据来源:%s不正确,请检查".format(data.data_source)}
            return json.dumps(info)
        Datasplit("datasplit").update_db_record("sg_split", {"_id": ObjectId(data.split_id)}, update_dict)
        self.set_sheet_data(name=task_name, options=options, table_id=task_sn, to_file=to_file, seq_number=seq_number)
        task_info = super(DatasplitV2Action, self).POST()
        # i = 0
        # while i < 3:
        #     i += 1
        #     back_result = json.loads(task_info)
        #     if "success" in back_result.keys():
        #         if back_result["success"] == False:
        #             time.sleep(30*i)
        #             task_info = super(DatasplitV2Action, self).POST()
        back_result = {"success":False,"info": "返回数据不是json或者dict：{}".format(task_info)}
        try:
            back_result = json.loads(task_info)
        except:
            # back_result = task_info
            if isinstance(task_info,dict):
                back_result=task_info
        if "success" in back_result.keys():
            if back_result["success"] == False:
                if data.data_source == "library":
                    split_status["first_split"] = "failed"
                    update_dict = {
                            "split_status": split_status,
                            "status": "failed",
                            "desc": "{}".format(back_result["info"])
                    }
                # elif data.data_source in ["specimen", "meta_raw"] or split_type == "second_split":
                elif data.data_source in ["specimen", "meta_raw"]:
                    split_status["second_split"] = "failed"
                    update_dict = {
                        "split_status": split_status,
                        "status": "failed",
                        "desc": "{}".format(back_result["info"])
                        }
                # elif data.data_source in ["other_raw", "qc"] or split_type == "qc":
                elif data.data_source in ["other_raw", "qc"]:
                    split_status["qc"] = "failed"
                    update_dict = {
                        "split_status": split_status,
                        "status": "failed",
                        "desc": "{}".format(back_result["info"])
                    }
                elif data.data_source == "cpc":
                    split_status["cpc"] = "failed"
                    update_dict = {
                        "split_status": split_status,
                        "cpc": "failed",
                        "desc": "{}".format(back_result["info"])
                }
                Datasplit("datasplit").update_db_record("sg_split", {"_id": ObjectId(data.split_id)}, update_dict)
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
