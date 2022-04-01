# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import web
import json
import datetime
from bson import SON
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metagbin import Metagbin
from mainapp.controllers.project.metagbin_controller import MetagbinController


class AssemblyPredictAction(MetagbinController):
    def __init__(self):
        super(AssemblyPredictAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        """
        设置接口及其参数
        :return:
        """
        data = web.input()
        print data
        default_argu = ['task_id', 'genome_id', 'taxon', 'submit_location', 'task_type', 'bam_file', 'ref_dir', 'insert', 'max']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)
        project_sn = self.metagbin.get_projectsn(data.task_id)
        task_name = 'metagbin.report.assembly_predict'
        module_type = 'workflow'
        assembly_name = "Assembly_" + data.genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            'genome_id': data.genome_id,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "soft": data.soft
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', assembly_name),
            ('genome_id', data.genome_id),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.metagbin.insert_main_table('assembly', mongo_data)
        update_info[str(main_id)] = 'assembly'

        options = {
            "project_sn": project_sn,
            "task_id": data.task_id,
            "insert_size": int(data.insert),
            "max_length": int(data.max),
            "bin_id": data.genome_id,
            "main_id": str(main_id),
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "ref_fa": data.ref_dir,
            'update_info': json.dumps(update_info),
            "bam_file": data.bam_file,
        }
        if hasattr(data, "taxon"):
            options["taxon"] = str(data.taxon)
        self.set_sheet_data(
                        name=task_name,
                        options=options,
                        #main_table_name= assembly_name.strip().split("_")[0] + '/' + assembly_name,
                        main_table_name="Assembly_predict" + '/' + assembly_name,
                        module_type=module_type,
                        project_sn=project_sn,
                        task_id=data.task_id,
                        params=params)
        task_info = super(AssemblyPredictAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': assembly_name}}
        return json.dumps(task_info)