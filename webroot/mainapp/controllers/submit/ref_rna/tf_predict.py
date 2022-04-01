# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class TfPredictAction(RefRnaController):
    """ TfPredict controller"""
    def __init__(self):
        super(TfPredictAction, self).__init__(instant=False)
        # specify where is the single workflow for the current task
        self.task_name = 'ref_rna.report.tf_predict'
        # specify all params needed, they will be saved as value of 'params' in the main table
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['s']
        self.expected_args += ['organism']
        self.expected_args += ['blast_all']
        # self.expected_args += ['seqfile']
        self.expected_args += ['E']
        self.expected_args += ['evalue']
        # receive the params from web
        self.input_data = web.input()
        # in which the main table will be created
        self.collection = "sg_tf_predict"

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(TfPredictAction,self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        # do some other check
        '''
        The following codes shows how to correct file path received from web
        target_dir = 'sanger' if self.input_data.client == "client01" else 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)
        trait_path = base_path + self.input_data.trait_path
        '''
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == "task_type":
                # maybe we do not need to do so
                params_dict[each] = int(input_data_dict[each])
            elif each == "group_dict":
                # deal with dumped dict param from web
                params_dict[each] = json.loads(input_data_dict[each], object_pairs_hook=OrderedDict)
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self, packed_params, collection_name):
        task_info = self.ref_rna.get_table_info_by_task_id(self.input_data.task_id, 'sg_task')
        self.project_sn = task_info["project_sn"]
        name = "TfPredict_{}_".format(self.input_data.E)
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            # type=self.input_data.?,  # You may need to save some other info for convenience
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='TfPredict main table',
            params=packed_params,
            status="start"
        )
        main_id = self.ref_rna.insert_main_table(collection_name, main_info)
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        sg_task_info = self.ref_rna.get_task_info(self.input_data.task_id)
        if 'is_demo' in sg_task_info and int(sg_task_info['is_demo']) != 0:
            task_id = sg_task_info['demo_id']
            result = self.ref_rna.get_table_info_by_task_id(task_id, 'sg_express_class_code')
        else:
            result = self.ref_rna.get_table_info_by_task_id(self.input_data.task_id, 'sg_express_class_code')
        if not result:
            raise Exception("query sg_express_class_code failed")
        seqdb_path = result['refrna_seqdb']
        class_code_id = str(result['_id'])
        options.update({
            "main_id": str(main_id),
            "update_info": json.dumps({str(main_id): collection_name}),  # to update sg_status
            "seq_db": seqdb_path+','+self.input_data.task_id,
            "class_code_id": class_code_id,
        })
        # prepare to file
        to_files = ["ref_rna.get_all_pep_seq(seq_db)", "ref_rna.get_gene_detail(class_code_id)"]
        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            to_file=to_files,
            project_sn=self.project_sn,
            task_id=self.input_data.task_id
        )


class TestFunction(unittest.TestCase):
    """
    This is test for the web_api func. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna/tf_predict "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_29268",
            task_type="2",  # maybe you need to change it
            submit_location="TfPredict",
            s="animal",
            organism="unknown",
            blast_all="yes",
            E="0.001",
            # domE="0.0001",
            evalue="0.0001",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

