# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao
# modified 20180906
import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class ChromosomeWindowAction(DnaController):
    """
    群体进化 基因注释接口
    """

    def __init__(self):
        super(ChromosomeWindowAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "step_num", "variant_type", "file_id", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3200301", "variables": var
}
                return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id,
                    "code": "C3200302", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        file_table = Dna("dna_evolution").find_one(collection="sg_variant_compare",
                                                 query_dic={"_id": ObjectId(data.file_id)})
        if not file_table:
            var = []
            var.append(data.file_id)
            info = {"success": False,
                    "info": "sg_variant_compare表里没有_id: %s对应的信息，请检查!" % data.file_id,
                    "code": "C3200303", "variables": var}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                file_path = file_table['pop_table_path']
            except:
                info = {"success": False,
                        "info": "sg_variant_compare中没有pop_table_path字段，请核查！", "code": "C3200304",
                        "variables": ""}
                return json.dumps(info)
        file_path = Dna("dna_evolution").set_file_path(data.task_id, file_path, data.client)
        params_json = {
            "task_type": int(data.task_type),
            "step_num": int(data.step_num),
            "variant_type": data.variant_type,
            "submit_location": data.submit_location,
            "file_id": data.file_id
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        main_table_name = 'Chromosome_Window_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        # main_table_name = Dna("dna_evolution").set_main_table_name("chromosome_window", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "基因组覆盖分布"),
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_chromosome_window", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_chromosome_window", query_dict={"_id": main_id},
                                              update_dict={"main_id": main_id, "compare_id": ObjectId(data.file_id)})
        update_info = {str(main_id): "sg_chromosome_window"}

        options = {
            "project_sn": project_sn,
            "main_id": str(main_id),
            "task_id": data.task_id,
            "step_num": data.step_num,
            "pop_table": file_path,
            "variant_type": data.variant_type,
            "update_info": json.dumps(update_info),    # 给work_flow用
        }
        self.set_sheet_data(name="dna_evolution.chromosome_window", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="ChromosomeWindows/" + main_table_name,
                            options=options, module_type="tool", params=params, db_type="dna_evolution",
                            analysis_name="ChromosomeWindows")
        task_info = super(ChromosomeWindowAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
