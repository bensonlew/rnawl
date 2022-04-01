# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.24

import re
import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bsa import Bsa
from mainapp.controllers.project.bsa_controller import BsaController


class SlidingwinFilterAnalysisAction(BsaController):
    """
    BSA关联区域定位接口
    lasted modified by hd 20180612
    """
    def __init__(self):
        super(SlidingwinFilterAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["slidingwin_id", "threshold", "threshold_value", "task_id"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1200401", "variables": var}
                return json.dumps(info)
        if data.threshold not in ["quantile", "index"]:
            var = []
            var.append(data.threshold)
            info = {"success": False, "info": "阈值确定类型：%s不对,只能为quantile/index!" % data.threshold,
                    "code": "C1200402", "variables": var}
            return json.dumps(info)
        if float(data.threshold_value) < 0 or float(data.threshold_value) > 1:
            var = []
            var.append(data.threshold_value)
            info = {"success": False, "info": "阈值确定范围在0-1之间:%s!" % data.threshold_value,
                    "code": "C1200403", "variables": var}
            return json.dumps(info)
        query_dic = {"task_id": data.task_id, "main_id": ObjectId(data.slidingwin_id)}
        result = Bsa().find_one_record(collection="sg_slidingwin", query_dic=query_dic)
        # noinspection PyBroadException
        try:
            slidingwin_file_ = result["slidingwin_result_path"]
            i_c_result_ = result["calc_index_path"]
            slid_file_ = result["slid_result_path"]
            mb = result["mb"]
            wb = result["wb"]
            mp = result["mp"]
            wp = result["wp"]
        except:
            info = {"success": False, "info": "sg_slidingwin表里信息不全，请检查!", "code": "C1200404", "variables": ""}
            return json.dumps(info)
        query_dic = {"task_id": data.task_id}
        result = Bsa().find_one_record(collection="sg_task", query_dic=query_dic)
        # noinspection PyBroadException
        try:
            pop_summary_ = result["pop_summary_path"]
            p_f_vcf_ = result["pop_vcf_path"]
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有pop_summary_path、pop_vcf_path信息，请检查!",
                    "code": "C1200405", "variables": ""}
            return json.dumps(info)
        pop_summary = Bsa().set_file_path(data.task_id, pop_summary_, data.client)
        p_f_vcf = Bsa().set_file_path(data.task_id, p_f_vcf_, data.client)
        slidingwin_file = Bsa().set_file_path(data.task_id, slidingwin_file_, data.client)
        i_c_result = Bsa().set_file_path(data.task_id, i_c_result_, data.client)
        slid_file = Bsa().set_file_path(data.task_id, slid_file_, data.client)
        # if re.match(r"^/mnt/ilustre/.*", pop_summary_):
        #     slidingwin_file = slidingwin_file_
        #     i_c_result = i_c_result_
        #     slid_file = slid_file_
        #     pop_summary = pop_summary_
        #     p_f_vcf = p_f_vcf_
        # else:
        #     sanger_path = 's3://'
        #     slidingwin_file = os.path.join(sanger_path, slidingwin_file_)
        #     i_c_result = os.path.join(sanger_path, i_c_result_)
        #     slid_file = os.path.join(sanger_path, slid_file_)
        #     pop_summary = os.path.join(sanger_path, pop_summary_)
        #     p_f_vcf = os.path.join(sanger_path, p_f_vcf_)
        params_json = {
            "slidingwin_id": data.slidingwin_id,
            "threshold": data.threshold,
            "threshold_value": data.threshold_value,
            "submit_location": "region",
            "task_type": 2
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "FilterAnalysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("slidingwin_id", ObjectId(data.slidingwin_id)),
            ("params", params),
            ("desc", "关联区域过滤主表"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("mb", mb),
            ("wb", wb),
            ("mp", mp),
            ("wp", wp)
        ]
        main_id = Bsa().insert_main_table(collection="sg_region", data=mongo_data)
        query_dict = {"_id": main_id}
        update_dict = {"main_id": main_id}
        Bsa().update_db_record(collection="sg_region", query_dict=query_dict, update_dict=update_dict)
        update_info = {str(main_id): "sg_region"}
        options = {
            "slidingwin_file": slidingwin_file,
            "slid_file": slid_file,
            "region_type": data.threshold,
            "region_value": data.threshold_value,
            "i_c_result": i_c_result,
            "pop_summary": pop_summary,
            "p_f_vcf": p_f_vcf,
            "mb": mb,
            "wb": wb,
            "mp": mp,
            "wp": wp,
            # "step": step,
            "update_info": json.dumps(update_info),
            "region_id": str(main_id)
        }
        sheet_data = self.set_sheet_data(name="bsa.report.slidingwin_filter_analysis", member_id=member_id,
                                         project_sn=project_sn, task_id=data.task_id,
                                         main_table_name="slidingwin_filter_analysis", options=options,
                                         module_type="workflow", params=params)
        print sheet_data
        task_info = super(SlidingwinFilterAnalysisAction, self).POST()
        return json.dumps(task_info)
