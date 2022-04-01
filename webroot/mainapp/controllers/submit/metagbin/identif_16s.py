# -*- coding: utf-8 -*-
#__author__ = 'haidong.gu'
import web
import re
import os
import json
import types
import datetime
from bson import SON
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.controllers.project.metagbin_controller import MetagbinController


class Identif16sAction(MetagbinController):
    def __init__(self):
        super(Identif16sAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        """
        设置接口及其参数
        :return:
        """
        data = web.input()
        print data
        default_argu = ['task_id', 'query', 'submit_location', 'task_type', 'db']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s'%argu}
                return json.dumps(info)
        project_sn = self.metagbin.get_projectsn(data.task_id)
        task_name = 'metagbin.report.taxon_identify'
        module_type = 'workflow'
        name = "16S_" + data.query +"_"+ datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            'query': data.query,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "db": data.db
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        update_info = {}

        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_table_id = Metagbin().insert_main_table("identif_16s", mongo_data)
        Metagbin().insert_main_table_new("identif_16s", str(main_table_id), {"main_id": main_table_id})
        update_info[str(main_table_id)] = 'identif_16s'
        options = {
            "update_info": json.dumps(update_info),
            "main_id": str(main_table_id),
            "bin_name": data.query,
            "method": "blasr",
            "query": Metagbin().get_sample_genefile(data.task_id, data.query, predict="rrna"),
            "blasr_db": data.db,
            "task_id": data.task_id
        }
        self.set_sheet_data(
            name=task_name,
            options=options,
            main_table_name= name.strip().split("_")[0] + '/' + name,
            module_type=module_type,
            project_sn=project_sn,
            task_id=data.task_id,
            params=params_json
        )
        task_info = super(Identif16sAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)