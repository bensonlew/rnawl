# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.26

import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class VennAnalysisAction(DnaController):
    """
    venn分析接口
    """
    def __init__(self):
        super(VennAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["set_ids", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101501", "variables": var}
                return json.dumps(info)
        print data.set_ids
        set_ids = json.loads(data.set_ids)
        print set_ids
        if len(set_ids) < 2 or len(set_ids) > 6:
            info = {"success": False, "info": "只能选择2-6个比较分析结果进行venn分析!",
                    "code": "C3101502", "variables": ""}
            return json.dumps(info)
        venn_type = []
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        for set_id in set_ids:
            result = Dna(db).find_one(collection="sg_site_set", query_dic={"task_id": data.task_id, "set_id": set_id})
            # noinspection PyBroadException
            try:
                type = result["type"]
                if type.startswith("snp"):
                    venn_type.append("snp")
                else:
                    venn_type.append("indel")
            except:
                var = []
                var.append(set_id)
                info = {"success": False, "info": "sg_site_set表里没有找到set_id:%s的信息，请检查!" % set_id,
                        "code": "C3101503", "variables": var}
                return json.dumps(info)
        venn_type = list(set(venn_type))
        print venn_type
        if len(venn_type) != 1:
            info = {"success": False, "info": "必须选择同属于snp或者同属于indel的比较结果进行venn分析，请检查!",
                    "code": "C3101504", "variables": ""}
            return json.dumps(info)
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到task_id:%s的project_sn, member_id信息，请检查!"
                                              % data.task_id, "code": "C3101505", "variables": var}
            return json.dumps(info)
        is_bucket = "true"
        if "region" in result.keys() and result['region']:
            rerewrwese_path = ""
        else:
            rerewrwese_path = "http://bcl.sanger.com/data/" if data.client == 'client01' \
                else "http://bcl.tsanger.com/data/"
        params_json = {
            "set_ids": set_ids,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "Venn_Analysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Venn分析主表"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_venn_analysis", data=mongo_data)
        Dna(db).update_db_record(collection="sg_venn_analysis", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_venn_analysis"}
        options = {
            "set_ids": ",".join(set_ids),
            "compare_type": venn_type[0],
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "is_bucket": is_bucket,
            "rerewrwese_path": rerewrwese_path
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "VennAnalysis/Venn_Analysis_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.venn_analysis", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="tool", params=params, db_type=db)
        task_info = super(VennAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
