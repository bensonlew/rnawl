# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180408

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class CnvCompareAction(DnaController):
    """
    cnv比较分析接口
    """
    def __init__(self):
        super(CnvCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sample1", "sample2", "is_same", "variation_type", "variation_len", "task_id",
                  "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100301", "variables": var}
                return json.dumps(info)
        if data.is_same not in ["true", "false"]:
            var = []
            var.append(data.is_same)
            info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" % data.is_same,
                    "code": "C3100302", "variables": var}
            return json.dumps(info)
        if data.variation_type not in ["", "deletion,duplication", "deletion", "duplication"]:
            var = []
            var.append(data.variation_type)
            info = {"success": False, "info": "变异位点类型%s不合法!" % data.variation_type,
                    "code": "C3100303", "variables": var}
            return json.dumps(info)
        params_json = {
            "sample1": data.sample1,
            "sample2": data.sample2,
            "is_same": data.is_same,
            "variation_type": data.variation_type,
            "variation_len": data.variation_len,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            cnv_anno_path = result['cnv_anno_path']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!",
                    "code": "C3100304", "variables": ""}
            return json.dumps(info)
        if cnv_anno_path.startswith("rerewrweset"):
            base_path = '/mnt/ilustre/data/' if str(data.client) == 'client01' else "/mnt/ilustre/tsanger-data/"
            cnv_anno_path_ = os.path.join(base_path, cnv_anno_path)
        elif cnv_anno_path.startswith("//"):
            region = result['region']
            cnv_anno_path_ = ":".join([region.rstrip(':'), cnv_anno_path])
        elif cnv_anno_path.startswith("/mnt") or re.match(".*://.*", cnv_anno_path):
            cnv_anno_path_ = cnv_anno_path
        else:
            var = []
            var.append(cnv_anno_path)
            info = {"success": False, "info": "存入的文件%s路径格式不正确！" % cnv_anno_path,
                    "code": "C3100305", "variables": var}
            return json.dumps(info)
        sample1 = os.path.join(cnv_anno_path_, data.sample1) + ".cnv.anno.xls"
        sample2 = os.path.join(cnv_anno_path_, data.sample2) + ".cnv.anno.xls"
        # if not os.path.exists(sample1):
        #     info = {"success": False, "info": "{}文件不存在！".format(sample1)}
        #     return json.dumps(info)
        # if not os.path.exists(sample2):
        #     info = {"success": False, "info": "{}文件不存在！".format(sample2)}
        #     return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'CnvCompare_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "cnv差异比较分析！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_cnv_compare", data=mongo_data)
        Dna(db).update_db_record(collection="sg_cnv_compare", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_cnv_compare"}
        options = {
            "sample1": sample1,
            "sample2": sample2,
            "is_same": data.is_same,
            "variation_len": data.variation_len,
            "variation_type": data.variation_type,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        main_table_name = "CnvCompare/CnvCompare_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs.report.cnv_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, params=params, db_type=db)
        task_info = super(CnvCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
