# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class RmatsDiffcompAction(RefRnaV2Controller):
    def __init__(self):
        super(RmatsDiffcompAction, self).__init__(instant=False)
        self.task_name = 'ref_rna_v2.report.rmats_diffcomp'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['delta_PSI']
        self.expected_args += ['significant_value']
        self.expected_args += ['significant_diff']
        self.expected_args += ['main_id']
        # self.expected_args +=['rmats_detail']
        # self.expected_args +=['compare_plan']
        # self.expected_args += ['last_id']
        self.input_data = web.input()
        self.collection = "sg_splicing_rmats_diffcomp"

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2902301',
                        'variables': variables}
                return json.dumps(info)
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        result_dir_path = []
        compare_plan = []
        for item in input_data_dict["main_id"].split(","):
            main_info = self.ref_rna_v2.get_main_info(item, "sg_splicing_rmats", self.input_data.task_id)
            result_dir_path.append(main_info["result_dir"] + "/all_events_detail_big_table.txt")
            compare_plan.append(main_info["compare_plan"])
        for each in self.expected_args:
            if each == "task_type":
                params_dict[each] = int(input_data_dict[each])
            # elif each == "rmats_detail":
            # params_dict[each] = ",".join(result_dir_path)
            # elif each == "last_id":
            #    last_info = self.ref_rna_v2.get_main_info_by_record("sg_rmats_diffcomp",task_id=self.input_data.task_id)
            #   if last_info is None:
            #       params_dict[each] = ""
            #    else:
            #       params_dict[each] = last_info["_id"]
            # elif each == "compare_plan":
            # params_dict[each] = ",".join(compare_plan)
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':')), result_dir_path, compare_plan

    def create_main_table(self, packed_params, collection_name):
        compare_plan_list = []
        result_info = self.ref_rna_v2.get_task_info(self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        name = 'splicing_diffcomp_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        for item in self.input_data.main_id.split(","):
            main_info = self.ref_rna_v2.get_main_info(item, "sg_splicing_rmats", self.input_data.task_id)
            compare_plan_list.append(main_info["compare_plan"])
        diffgroups = ",".join(compare_plan_list)
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='splicing_diffcomp_table',
            params=packed_params,
            status="start",
            diffgroups=diffgroups,
        )
        main_id = self.ref_rna_v2.insert_main_table(collection_name, main_info)
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        last_info = self.ref_rna_v2.get_main_info_by_record("sg_splicing_rmats_diffcomp",
                                                            task_id=self.input_data.task_id)
        if last_info is None:
            last_id = ""
        elif last_info["status"] == "end":
            last_id = last_info["_id"]
        else:
            last_id = ""
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            if each not in input_data_dict.keys():
                pass
            else:
                options[each] = input_data_dict[each]
        new_task_id = self.ref_rna_v2.get_new_id(self.input_data.task_id)
        main_table_data = {'run_id': new_task_id}
        options.update({
            "main_id": str(main_id),
            "main_table_data": main_table_data,
            "update_info": json.dumps({str(main_id): collection_name}),
            "last_id": str(last_id),
            "rmats_detail": ",".join(self.pack_params()[1]),
            "compare_plan": ",".join(self.pack_params()[2])
        })
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,
            module_type="workflow",
            project_sn=self.project_sn,
            new_task_id=new_task_id,
            task_id=self.input_data.task_id
        )

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()[0]
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(RmatsDiffcompAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.ref_rna_v2.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.ref_rna_v2.update_group_compare_is_use(self.input_data.task_id, self.input_data.control_id)
        return json.dumps(task_info)
