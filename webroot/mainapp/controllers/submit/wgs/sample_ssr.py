# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

import re
import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SampleSsrAction(DnaController):
    """
    样本基因组SSR分析接口
    """
    def __init__(self):
        super(SampleSsrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["group_id", "group_dict", "tm1", "tm2", "product_size", "primer_num", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101001", "variables": var}
                return json.dumps(info)
        if float(data.tm1) >= float(data.tm2):
            var = []
            var.append(data.tm1)
            var.append(data.tm2)
            info = {"success": False, "info": "tm1:%s不能大于等于tm2:%s,请检查!" % (data.tm1, data.tm2),
                    "code": "C3101002", "variables": var}
            return json.dumps(info)
        if int(data.primer_num) > 5 or int(data.primer_num) < 1:
            var = []
            var.append(data.primer_num)
            info = {"success": False, "info": "primer_num:%s范围只能是[2,5],请检查!" % data.primer_num,
                    "code": "C3101003", "variables": var}
            return json.dumps(info)
        specimen_names = json.loads(data.group_dict)
        specimen_ids = []
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        for group in specimen_names.keys():
            sample_names = specimen_names[group]
            for sample_name in sample_names:
                specimen_ids.append(sample_name)
                result = Dna(db).find_one(collection="sg_specimen", query_dic={"task_id": data.task_id, "old_name": sample_name})
                if not result:
                    var = []
                    var.append(data.primer_num)
                    info = {"success": False, "info": "sg_specimen表里没有找到样本:%s,请检查!" % sample_name,
                            "code": "C3101004", "variables": var}
                    return json.dumps(info)
        specimen_ids = list(set(specimen_ids))
        print specimen_ids
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            clean_fastq_path = result['clean_fastq_path']
        except:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到task_id:%s的project_sn, member_id信息，请检查!"
                                              % data.task_id, "code": "C3101005", "variables": var}
            return json.dumps(info)
        clean_fastq_path = Dna(db).set_file_path(data.task_id, clean_fastq_path, data.client)
        params_json = {
            "group_id": data.group_id,
            "group_dict": json.loads(data.group_dict),
            "tm1": data.tm1,
            "tm2": data.tm2,
            "product_size": data.product_size,
            "primer_num": int(data.primer_num),
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "Ssr_Analysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "样本基因组SSR分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        ssr_id = Dna(db).insert_main_table(collection="sg_ssr_specimen", data=mongo_data)
        Dna(db).update_db_record(collection="sg_ssr_specimen", query_dict={"_id": ssr_id},
                                 update_dict={"main_id": ssr_id})
        update_info = {str(ssr_id): "sg_ssr_specimen"}
        is_bucket = "false"
        if "region" in result.keys() and result['region']:
            is_bucket = "true"
        options = {
            "specimen_names": ','.join(specimen_ids),
            "tm1": data.tm1,
            "tm2": data.tm2,
            "product_size": data.product_size,
            "primer_num": data.primer_num,
            "update_info": json.dumps(update_info),
            "main_id": str(ssr_id),
            "task_id": data.task_id,
            "clean_fastq_path": clean_fastq_path,
            "is_bucket": is_bucket
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "SampleSsr/Ssr_Analysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.report.sample_ssr", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params, db_type=db)
        task_info = super(SampleSsrAction, self).POST()
        task_info['id'] = str(ssr_id)
        return json.dumps(task_info)
