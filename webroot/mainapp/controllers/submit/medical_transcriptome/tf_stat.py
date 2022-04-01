# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class TfStatAction(MedicalTranscriptomeController):
    """ TfStat controller"""
    def __init__(self):
        super(TfStatAction, self).__init__(instant=False)
        # specify where is the single workflow for the current task
        self.task_name = 'medical_transcriptome.report.tf_stat'
        # specify all params needed, they will be saved as value of 'params' in the main table
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['geneset_id']
        self.expected_args += ['families']
        self.expected_args += ['level']

        self.expected_args += ['tf_predict_id']
        self.expected_args += ['e_value']
        # receive the params from web
        self.input_data = web.input()
        # in which the main table will be created
        self.collection = "sg_tf_stat"

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(TfStatAction,self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.medical_transcriptome.update_group_compare_is_use(self.input_data.task_id, self.input_data.control_id)
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903101', 'variables': variables}
                return json.dumps(info)
        # do some other check
        main_info = self.medical_transcriptome.get_main_info(self.input_data.tf_predict_id, 'sg_tf_predict', self.input_data.task_id)
        if 'E' in json.loads(main_info['params']):
            pre_evalue = float(json.loads(main_info['params'])['E'])
        else:
            pre_evalue = float(json.loads(main_info['params'])['evalue'])
        if float(pre_evalue) < float(self.input_data.e_value):
            variables = []
            variables.append(pre_evalue)
            info = {'success': False, 'info': "evalue阈值不能大于上次TF预测时使用的阈值%s" % pre_evalue, 'code': 'C2903102', 'variables': variables}
            return json.dumps(info)
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
        return json.dumps(params_dict, sort_keys=True, separators=(',',':'))

    def create_main_table(self, packed_params, collection_name):
        result_info = self.medical_transcriptome.get_task_info(self.input_data.task_id)
        print(result_info)
        self.project_sn = result_info["project_sn"]
        name = "TfStat_"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            # type=self.input_data.?,  # You may need to save some other info for convenience
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='TfStat main table',
            params=packed_params,
            status="start",
            tf_predict_id=ObjectId(self.input_data.tf_predict_id)
        )
        main_id = self.medical_transcriptome.insert_main_table(collection_name, main_info)
        if self.input_data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.medical_transcriptome.insert_geneset_info(self.input_data.geneset_id, collection_name, str(main_id))
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        new_task_id = self.medical_transcriptome.get_new_id(self.input_data.task_id)
        main_table_data = {'run_id': new_task_id}
        options.update({
            "main_id": str(main_id),
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): collection_name})  # to update sg_status
            # "tofile_Param": self.input_data.geneset_id,
        })
        # prepare to file
        # to_files = ["ref_rna_v2.YOUR_tofile_Func(tofile_Param)", ]
        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            # to_file=to_files,
            project_sn=self.project_sn,
            new_task_id=new_task_id,
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
        cmd += "s/medical_transcriptome/tf_stat "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",  # maybe you need to change it
            submit_location="TfStat",
            geneset_id="5f45c5b917b2bf77c0eb227d",
            families="ETS,THAP,AF-4,AP-2",
            tf_predict_id="5f9a174e17b2bf72fea69d3a",
            e_value="0.00001",
            level='G'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
