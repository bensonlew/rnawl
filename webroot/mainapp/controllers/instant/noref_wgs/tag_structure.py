# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.07

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class TagStructureAction(DnaController):
    """
    Tag结构图
    """
    def __init__(self):
        super(TagStructureAction, self).__init__(instant=True)

    def POST(self):
        data = web.input()
        print data
        params = ["tag_id", "task_id", "project_type", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        if not data.tag_id:
            info = {"success": False, "info": "参数tag_id不能为空!"}
            return json.dumps(info)
        db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            member_id = result["member_id"]
            project_sn = result["project_sn"]
            populations_tag = result['populations_tag']
        except:
            info = {"success": False, "info": "sg_task表里没有populations_tag信息，请检查!"}
            return json.dumps(info)
        params_json = {
            "tag_id": data.tag_id,
            "task_id": data.task_id,
            "task_type": int(data.task_type),
            "project_type": data.project_type,
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "tag_structure_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("status", "start"),
            ("params", params),
            ("tag_id", data.tag_id),
            ("desc", "tag结构图"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_tag", data=mongo_data)
        Dna(db).update_db_record(collection="sg_tag", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_tag"}
        options = {
            "populations_tag": populations_tag,
            "tag_id": data.tag_id,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        main_table_name = "tag_structure_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="noref_wgs.tag_structure", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="tool", params=params, db_type=db)
        task_info = super(TagStructureAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
