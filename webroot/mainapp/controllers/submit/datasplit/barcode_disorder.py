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


class BarcodeDisorderAction(DatasplitController):
    """
    多样性结果验证-barcode错乱
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
        result = Datasplit("datasplit").coll_find_one("sg_split_library",
            {"split_id": ObjectId(split_id), "library_number": library_number})
        if not result:
            info = {"success": False, "info": "没有在表sg_split_library里找到library_number: %s,请检查" % library_number}
            return json.dumps(info)
        lib_type = "no_official"
        if re.search("双index官方多样性文库", lib_type):
            lib_type = "official"
        fq_path = result["work_path"].split(";")[0]
        if not os.path.exists(fq_path):
            fq_path = result["path"].split(";")[0]
        if not fq_path:
            info = {"success": False, "info": "library_number: %s的path为空,请检查" % library_number}
            return json.dumps(info)
        lib_insert_size = result["insert_size"]
        main_table_name = "barcode_disorder_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            "split_id": split_id,
            "library_number": library_number,
            "lib_type": lib_type,
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        data_dict = {
            "name": main_table_name,
            "status": "start",
            "task_sn": main_table_name,
            "split_id": split_id,
            "library_number": library_number,
            "params": params,
            "type": "barcode_disorder",
            "desc": "多样性结果验证-barcode错乱",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        query_dict = {"split_id": str(split_id), "library_number": library_number}
        result = Datasplit("datasplit").coll_find_one("sg_meta_verify_barcode", query_dict)
        if result:
            verify_id = str(result["_id"])
            if result["status"] == "end":
                info = {"success": True, "info": "library_number: %s barcode错配运行成功" % library_number}
                return json.dumps(info)
            query_dict = {"_id": ObjectId(verify_id)}
            update_dict = data_dict
            Datasplit("datasplit").update_db_record("sg_meta_verify_barcode", query_dict, update_dict)
        else:
            verify_id = Datasplit("datasplit").insert_one("sg_meta_verify_barcode", data_dict)
        update_info = {str(verify_id): "sg_meta_verify_barcode"}
        options = {
            "fq_dir": os.path.dirname(fq_path) + "/",
            "barcode_primer_info": str(verify_id),
            "library_number": library_number,
            "lib_insert_size": lib_insert_size,
            "lib_type": lib_type,
            "update_info": json.dumps(update_info),
            "verify_id": str(verify_id)
        }
        to_file = ["datasplit_v2.export_barcode_disorder(barcode_primer_info)"]
        self.set_sheet_data(name="datasplit_v2.barcode_disorder", options=options, table_id="barcode_"+library_number, to_file=to_file, seq_number="", module_type="workflow")
        task_info = super(BarcodeDisorderAction, self).POST()
        return json.dumps(task_info)
