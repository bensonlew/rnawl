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


class PrimerMismatchAction(DatasplitController):
    """
    多样性结果验证-引物错配
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["split_id", "library_number"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        split_id = data.split_id
        library_number = data.library_number
        query_dict = {"split_id": ObjectId(split_id), "library_number": library_number}
        results = Datasplit("datasplit").coll_find("sg_split_specimen", query_dict)
        if not results:
            info = {"success": False, "info": "没有在表sg_split_specimen里找到library_number: %s,请检查" % library_number}
            return json.dumps(info)
        fq_path = ""
        for result in results:
            if "clean_work_path" in result.keys():
                fq_path = result["clean_work_path"].split(";")[0]
                break
            elif "work_path" in result.keys():
                fq_path = result["work_path"].split(";")[0]
                break
        if fq_path == "":
            info = {"success": False, "info": "没有在表sg_split_specimen里找到work_path,请检查文库:%s是否有拆分结果" % library_number}
            return json.dumps(info)
        qc_dir = "/".join(fq_path.split("/")[:-3]) + "/MetaQc"
        module_info = "/".join(fq_path.split("/")[:-3]) + "/MetaQc/module_workdir.info"
        # if not os.path.exists(qc_dir):
        #     info = {"success": False, "info": "library_number: %s的workspace: %s不存在,不能进行引物错配检查" % (library_number, qc_dir)}
        #     return json.dumps(info)
        # if not os.path.exists(module_info):
        #     info = {"success": False, "info": "library_number: %s的文件: %s不存在,不能进行引物错配检查" % (library_number, module_info)}
        #     return json.dumps(info)
        main_table_name = "primer_mismatch_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            "split_id": split_id,
            "library_number": library_number,
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        data_dict = {
            "name": main_table_name,
            "status": "start",
            "task_sn": main_table_name,
            "split_id": split_id,
            "library_number": library_number,
            "params": params,
            "type": "primer_mismatch",
            "desc": "多样性结果验证-primer错配",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        query_dict = {"split_id": str(split_id), "library_number": library_number}
        result = Datasplit("datasplit").coll_find_one("sg_meta_verify_primer", query_dict)
        if result:
            verify_id = str(result["_id"])
            if result["status"] == "end":
                info = {"success": True, "info": "library_number: %s primer错位运行成功" % library_number}
                return json.dumps(info)
            query_dict = {"_id": ObjectId(verify_id)}
            update_dict = data_dict
            Datasplit("datasplit").update_db_record("sg_meta_verify_primer", query_dict, update_dict)
        else:
            verify_id = Datasplit("datasplit").insert_one("sg_meta_verify_primer", data_dict)
        update_info = {str(verify_id): "sg_meta_verify_primer"}
        options = {
            "lib_list": str(verify_id),
            "update_info": json.dumps(update_info),
            "verify_id": str(verify_id)
        }
        to_file = ["datasplit_v2.export_primer_mismatch(lib_list)"]
        self.set_sheet_data(name="datasplit_v2.primer_mismatch", options=options, table_id="primer_"+library_number, to_file=to_file, seq_number="", module_type="workflow")
        task_info = super(PrimerMismatchAction, self).POST()
        return json.dumps(task_info)
