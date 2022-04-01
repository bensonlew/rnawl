# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.09

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SvCompareAction(DnaController):
    """
    SV比较分析接口
    """
    def __init__(self):
        super(SvCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sample1", "sample2", "is_same", "variation_type", "variation_len", "sample1_support",
                  "sample2_support", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101301", "variables": var}
                return json.dumps(info)
        if data.is_same not in ["true", "false"]:
            var = []
            var.append(data.is_same)
            info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" % data.is_same,
                    "code": "C3101302", "variables": var}
            return json.dumps(info)
        for t in data.variation_type.split(","):
            if t not in ["DEL", "INV", "ITX", "CTX", "INS"]:
                var = []
                var.append(data.variation_type)
                info = {"success": False, "info": "变异位点类型%s不合法!" % data.variation_type,
                        "code": "C3101303", "variables": var}
                return json.dumps(info)
        m = re.match(r"(.*)-(.*)", data.variation_len)
        if not m:
            var = []
            var.append(data.variation_len)
            info = {"success": False, "info": "变异区域长度%s不合法!" % data.variation_len,
                    "code": "C3101304", "variables": var}
            return json.dumps(info)
        elif data.variation_len == "0-0":
            variation_len = None
        else:
            variation_len = m.group(1)
            if m.group(2):
                variation_len += ":" + m.group(2)
        n = re.match(r"(.*)-(.*)", data.sample1_support)
        if not n:
            var = []
            var.append(data.sample1_support)
            info = {"success": False, "info": "Reads支持度%s不合法!" % data.sample1_support,
                    "code": "C3101305", "variables": var}
            return json.dumps(info)
        elif data.sample1_support == "0-0":
            sample1_support = None
        else:
            sample1_support = n.group(1)
            if n.group(2):
                sample1_support += ":" + n.group(2)
        l = re.search(r"(.*)-(.*)", data.sample2_support)
        if not l:
            var = []
            var.append(data.sample2_support)
            info = {"success": False, "info": "Reads支持度%s不合法!" % data.sample2_support,
                    "code": "C3101306", "variables": var}
            return json.dumps(info)
        elif data.sample2_support == "0-0":
            sample2_support = None
        else:
            sample2_support = l.group(1)
            if l.group(2):
                sample2_support += ":" + l.group(2)
        params_json = {
            "sample1": data.sample1,
            "sample2": data.sample2,
            "is_same": data.is_same,
            "variation_type": data.variation_type,
            "variation_len": data.variation_len,
            "sample1_support": data.sample1_support,
            "sample2_support": data.sample2_support,
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
        try:
            sv_anno_path = result['sv_anno_path']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!",
                    "code": "C3101307", "variables": ""}
            return json.dumps(info)
        # base_path = "s3://"
        # if sv_anno_path.startswith("rerewrweset"):
        #     sv_anno_path = os.path.join(base_path, sv_anno_path)
        sv_anno_path = Dna(db).set_file_path(data.task_id, sv_anno_path, data.client)
        sample1_sv_anno = os.path.join(sv_anno_path, data.sample1) + ".sv.anno.xls"
        sample2_sv_anno = os.path.join(sv_anno_path, data.sample2) + ".sv.anno.xls"
        # if not os.path.exists(sample1_sv_anno):
        #     info = {"success": False, "info": "{}文件不存在！".format(sample1_sv_anno)}
        #     return json.dumps(info)
        # if not os.path.exists(sample2_sv_anno):
        #     info = {"success": False, "info": "{}文件不存在！".format(sample2_sv_anno)}
        #     return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "SvCompare_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "sv差异比较分析！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_sv_compare", data=mongo_data)
        Dna(db).update_db_record(collection="sg_sv_compare", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_sv_compare"}
        options = {
            "sample1_sv_anno": sample1_sv_anno,
            "sample2_sv_anno": sample2_sv_anno,
            "is_same": data.is_same,
            "variation_len": variation_len,
            "variation_type": data.variation_type,
            "sample1_support": sample1_support,
            "sample2_support": sample2_support,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs.sv_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name="sv_compare", options=options, module_type="tool", params=params, db_type=db)
        task_info = super(SvCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
