# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.23

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bsa import Bsa
from mainapp.controllers.project.bsa_controller import BsaController


class SlidingwinAnalysisAction(BsaController):
    """
    BSA标记筛选和分析接口
    lasted modified by hd 20180612
    """
    def __init__(self):
        super(SlidingwinAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sliding_type", "sliding_strategy", "s_deep", "b_deep", "task_id"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1200301", "variables": var }
                return json.dumps(info)
        if data.sliding_type not in ["distance", "variant_num"]:
            var =[]
            var.append(data.sliding_type)
            info = {"success": False, "info": "滑窗方式%s需为distance/variant_num!" % data.sliding_type,
                    "code": "C1200302", "variables": var}
            return json.dumps(info)
        if data.sliding_type == "distance":
            sliding_type = "bp"
            m = re.match(r"(.*)-(.*)", data.sliding_strategy)   # modified by hd 20180423 解决传进来的是小数
            if m:
                win = float(m.group(1)) * 1000000
                step = float(m.group(2)) * 1000
                if win < step:
                    var =[]
                    var.append(win)
                    var.append(step)
                    info = {"success": False, "info": "滑窗策略M前数值%s需大于K前数值%s!" % (win, step),
                            "code": "C1200303", "variables": var}
                    return json.dumps(info)
                win = str(win)
                step = str(step)
            else:
                var = []
                var.append(data.sliding_strategy)
                info = {"success": False, "info": "滑窗策略%s必须是'数值M-数值K'的形式!" % data.sliding_strategy,
                        "code": "C1200304", "variables": var}
                return json.dumps(info)
        else:
            sliding_type = "num"
            m = re.match(r"(\d+)-(\d+)", data.sliding_strategy)
            if m:
                win = int(m.group(1))
                step = int(m.group(2))
                if win < step:
                    var = []
                    var.append(win)
                    var.append(step)
                    info = {"success": False, "info": "滑窗策略中划线前数值%s需大于中划线后数值%s!" % (win, step),
                            "code": "C1200305", "variables": var}
                    return json.dumps(info)
                win = str(win)
                step = str(step)
            else:
                var = []
                var.append(data.sliding_strategy)
                info = {"success": False, "info": "滑窗策略%s必须是'正整数-正整数'的形式!" % data.sliding_strategy,
                        "code": "C1200306", "variables": var}
                return json.dumps(info)
        query_dic = {"task_id": data.task_id}
        result = Bsa().find_one_record(collection="sg_task", query_dic=query_dic)
        # noinspection PyBroadException
        try:
            pop_vcf_path_ = result["pop_vcf_path"]
            pop_summary_path_ = result["pop_summary_path"]
            ref_chrlist_ = result["ref_chrlist"]
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有pop_vcf_path、pop_summary_path信息，请检查!",
                    "code": "C1200307", "variables": ""}
            return json.dumps(info)
        query_dic = {"task_id": data.task_id}
        result = Bsa().find_one_record(collection="sg_slidingwin", query_dic=query_dic)
        # noinspection PyBroadException
        try:
            mb = result["mb"]
            wb = result["wb"]
            mp = result["mp"]
            wp = result["wp"]
            variant_type = result["variant_type"]
        except:
            info = {"success": False, "info": "sg_slidingwin表里没有mb、wb、mp、wp、variant_type信息，请检查!",
                    "code": "C1200308", "variables": ""}
            return json.dumps(info)
        pop_summary_path = Bsa().set_file_path(data.task_id, pop_summary_path_, data.client)
        pop_vcf_path = Bsa().set_file_path(data.task_id, pop_vcf_path_, data.client)
        ref_chrlist = Bsa().set_file_path(data.task_id, ref_chrlist_, data.client)
        params_json = {
            "s_deep": int(data.s_deep),
            "b_deep": int(data.b_deep),
            "sliding_type": data.sliding_type,
            "sliding_strategy": data.sliding_strategy,
            "submit_location": "slidingwin",
            "task_type": 2
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "SlidingwinAnalysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("member_id", member_id),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "标记筛选和分析主表"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("mb", mb),
            ("wb", wb),
            ("mp", mp),
            ("wp", wp),
            ("variant_type", variant_type)
        ]
        main_id = Bsa().insert_main_table(collection="sg_slidingwin", data=mongo_data)
        query_dict = {"_id": main_id}
        update_dict = {"main_id": main_id}
        Bsa().update_db_record(collection="sg_slidingwin", query_dict=query_dict, update_dict=update_dict)
        update_info = {str(main_id): "sg_slidingwin"}
        options = {
            "pop_vcf": pop_vcf_path,
            "pop_summary": pop_summary_path,
            "ref_chrlist": ref_chrlist,
            "mb": mb,
            "wb": wb,
            "mp": mp,
            "wp": wp,
            "pdep": int(data.s_deep),
            "bdep": int(data.b_deep),
            "method": sliding_type,
            "variant_type": variant_type,
            "win": win,
            "step": step,
            "update_info": json.dumps(update_info),
            "slidingwin_id": str(main_id)
        }
        main_table_name = "SlidingwinAnalysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        sheet_data = self.set_sheet_data(name="bsa.report.slidingwin_analysis", member_id=member_id,
                                         project_sn=project_sn, task_id=data.task_id, main_table_name=main_table_name,
                                         options=options, module_type="workflow", params=params)
        print sheet_data
        task_info = super(SlidingwinAnalysisAction, self).POST()
        return json.dumps(task_info)
