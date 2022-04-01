# -*- coding: utf-8 -*-

import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime


class BugbaseContributionAction(MetaasvController):

    def __init__(self):
        super(BugbaseContributionAction, self).__init__(instant=False)


    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['asv_id', 'bugbase_id', 'submit_location', 'group_detail', 'group_id']

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu], "code" : "C2204501"}
                return json.dumps(info)

        task_name = 'metaasv.report.bugbase_contribution'
        task_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'bugbase_id': data.bugbase_id,
            'group_id': data.group_id,
            'method': data.method,
            'top': data.top,
            'species_level': data.species_level,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }

        main_table_name = 'BugBaseContribution_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('bugbase_id', ObjectId(data.bugbase_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_main_table('bugbase_contribution', mongo_data)

        update_info = {str(main_table_id): 'bugbase_contribution'}
        options = {
            'otu_table': data.asv_id,
            'asv_id': data.asv_id,
            'bugbase_id': data.bugbase_id,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'bugbase_table': data.bugbase_id,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'tax_file': data.asv_id,
            'species_level': data.species_level,
            'method': data.method,
            'top': data.top,
            'normalized_table': data.bugbase_id
        }

        to_file = ["metaasv.export_otu_table_without_tax(otu_table)","metaasv.export_group_table_by_detail_unsort(group_table)","metaasv.export_bugbase_contribution_table(bugbase_table)","metaasv.export_tax_table_by_asv_id(tax_file)", "metaasv.export_normalized_table_by_bugbase_id(normalized_table)"]

        self.set_sheet_data(name=task_name, options=options, main_table_name="BugBaseContribution/"+main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(BugbaseContributionAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
