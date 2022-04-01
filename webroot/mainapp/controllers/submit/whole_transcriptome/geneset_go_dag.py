# -*- coding: utf-8 -*-

import datetime
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.ref_rna import *


class GenesetGoDagAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetGoDagAction, self).__init__(instant=False)
        self.task_name = 'whole_transcriptome.report.geneset_go_dag'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        # self.expected_args += ['go_list']
        # self.expected_args += ['top_num']
        # self.expected_args += ['significant_diff']
        # self.expected_args += ['significant_value']
        self.expected_args += ['go_enrich_id']
        self.input_data = web.input()
        self.collection = "geneset_go_dag"

    def check_params(self):
        # 获取geneset_id
        go_enrich_id = self.input_data.go_enrich_id
        result = self.db['geneset_go_enrich'].find_one({'main_id': ObjectId(go_enrich_id)})
        if result:
            params = json.loads(result['params'])
            self.geneset_id = params['geneset_id']
        if hasattr(self.input_data, "go_list"):
            self.expected_args += ['go_list']
        elif hasattr(self.input_data, "significant_diff") and hasattr(self.input_data, "significant_value") and hasattr(
                self.input_data, "top_num"):
            self.expected_args += ['top_num']
            self.expected_args += ['significant_diff']
            self.expected_args += ['significant_value']
        else:
            info = {"success": False, "info": u"GO列表和筛选阈值至少填入一个", 'code': 'C2901401', 'variables': ''}
            return info
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901402', 'variables': variables}
                return json.dumps(info)
            else:
                return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        # result_dir_path = []
        # main_info = self.ref_rna_v2.get_main_info(input_data_dict["main_id"], "geneset_go_enrich", self.input_data.task_id)
        for each in self.expected_args:
            if each == "task_type":
                params_dict[each] = int(input_data_dict[each])
            elif each == "significant_value":
                params_dict[each] = float(input_data_dict[each])
            elif each == "top_num":
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self, packed_params, collection_name):
        result_info = self.whole_transcriptome.get_task_info(self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        name = 'GO_DAG_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            # project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset_go_dag_table',
            params=packed_params,
            status="start",
            version='v1'
        )
        main_id = self.whole_transcriptome.insert_main_table("geneset_go_dag", main_info)
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        main_info = self.whole_transcriptome.get_main_info(input_data_dict["go_enrich_id"], "geneset_go_enrich",
                                                           self.input_data.task_id)
        result_dir = main_info["result_dir"]
        options = dict()
        for each in self.expected_args:
            if input_data_dict[each] == "":
                options[each] = None
            else:
                options[each] = input_data_dict[each]
        options.update({
            "go_enrich_id": str(main_id),
            "update_info": json.dumps({str(main_id): collection_name}),
            "go_enrich_detail": result_dir
        })
        task_info = self.whole_transcriptome.get_task_info(self.input_data.task_id)
        if "database_version" in task_info["long_task"]:
            go_version = task_info["long_task"]["database_version"].get("go", "20200628").split("_")[0]
            options.update({"go_version": go_version})

        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,
            module_type="workflow",
            project_sn=self.project_sn,
            task_id=self.input_data.task_id
        )

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        self.whole_transcriptome.delete_main_table("geneset_go_dag", self.input_data.task_id)
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(GenesetGoDagAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.whole_transcriptome.update_group_compare_is_use(self.input_data.task_id,
                                                                     self.input_data.control_id)
        self.whole_transcriptome.insert_geneset_info(self.geneset_id, 'geneset_go_dag', str(main_id))
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        import os
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/whole_transcriptome/geneset_go_dag "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="whole_transcriptome",
            submit_location="genesetgo_acyclic",
            task_type="2",
            go_enrich_id="5dc8f63517b2bf5b6cf33055",
            significant_diff="pvalue",
            significant_value="0.05",
            top_num="20",
            # go_list="GO:0043631,GO,0043487,GO:0050906,GO:0031572,GO:0090329,GO:0007157,GO:0006941,GO:0009124,GO:0022402,GO:0051606,GO:0009064"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
