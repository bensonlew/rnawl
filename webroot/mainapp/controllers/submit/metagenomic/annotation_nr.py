# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
import types
from mainapp.models.mongo.metagenomic import Metagenomic
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from bson import SON
import os
from mainapp.libs.signature import check_sig
from mbio.packages.metagenomic.id_convert import id2name
from mbio.packages.metagenomic.id_convert import name2id
from mainapp.controllers.project.metagbin_controller import MetagbinController


class AnnotationNrAction(MetagenomicController):
    def __init__(self):
        super(AnnotationNrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print (data)
        default_argu = ["query", "upload", "task_id", "project_sn"]
        #, "identity", "align_length", "task_type",

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables": [argu], "code" : "C2403801"}
                return json.dumps(info)
        task_name = "metagenomic.report.nr_annotation"#调用的workflow
        module_type = "workflow"
        project_sn = data.project_sn
        task_id = data.task_id
        #metagenomic = Metagenomic()
        name = "NR_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "query": data.query,
            "upload": data.upload,
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        options = {
            "query": data.query,
            "upload": data.upload,
        }
        print options
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=task_id,
                            params=params)
        task_info = super(AnnotationNrAction, self).POST()
        if task_info['success']:
            print  task_info
            return json.dumps(task_info)
