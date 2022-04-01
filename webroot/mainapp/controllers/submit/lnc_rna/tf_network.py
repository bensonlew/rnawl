# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class TfNetworkAction(LncRnaController):
    """ TfNetwork controller"""
    def __init__(self):
        super(TfNetworkAction, self).__init__(instant=False)
        # specify where is the single workflow for the current task
        self.task_name = 'lnc_rna.report.tf_network'
        # specify all params needed, they will be saved as value of 'params' in the main table
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['tfbs_predict_id']
        self.expected_args += ['thresh']
        self.expected_args += ['qv_thresh']
        self.expected_args += ['gene_type']
        self.expected_args += ['corr_pvalue']
        self.expected_args += ['corr_cutoff']
        self.expected_args += ['corr_use_pdajust']
        # receive the params from web
        self.input_data = web.input()
        # in which the main table will be created
        self.collection = "sg_tf_network"

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2902901', 'variables': variables}
                return json.dumps(info)
        # do some other check
        main_info = self.lnc_rna.get_main_info(self.input_data.tfbs_predict_id, 'sg_tfbs_predict', self.input_data.task_id)
        tfbs_params = json.loads(main_info['params'])
        pre_corr_cutoff = float(tfbs_params['corr_cutoff'])
        pre_corr_pvalue = float(tfbs_params['corr_pvalue'])
        pre_thresh = float(tfbs_params['thresh'])
        if pre_thresh < float(self.input_data.thresh):
            variables = []
            variables.append(pre_thresh)
            info = {'success': False, 'info': "FIMO对应的pvalue阈值不能大于上次TFBS预测时使用的阈值%s" % pre_thresh, 'code': 'C2902902', 'variables': variables}
            return json.dumps(info)
        if pre_corr_cutoff > float(self.input_data.corr_cutoff):
            variables = []
            variables.append(pre_corr_cutoff)
            info = {'success': False, 'info': "相关系数的阈值不能小于上次TFBS预测时使用的阈值%s" % pre_corr_cutoff, 'code': 'C2902903', 'variables': variables}
            return json.dumps(info)
        if pre_corr_pvalue < float(self.input_data.corr_pvalue):
            variables = []
            variables.append(pre_corr_pvalue)
            info = {'success': False, 'info': "相关系数的pvalue阈值不能大于上次TFBS预测时使用的阈值%s" % pre_corr_pvalue, 'code': 'C2902904', 'variables': variables}
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
        result_info = self.lnc_rna.get_task_info(self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        name = "TfNetwork"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            # type=self.input_data.?,  # You may need to save some other info for convenience
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='TfNetwork main table',
            params=packed_params,
            status="start",
            tfbs_predict_id=ObjectId(self.input_data.tfbs_predict_id)
        )
        main_id = self.lnc_rna.insert_main_table(collection_name, main_info)
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        options.update({
            "main_id": str(main_id),
            "update_info": json.dumps({str(main_id): collection_name}),  # to update sg_status
            "tfbs_predict": self.input_data.tfbs_predict_id,
        })
        # prepare to file
        to_files = ["lnc_rna.export_predict_result(tfbs_predict)"]
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

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(TfNetworkAction,self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.lnc_rna.update_group_compare_is_use(self.input_data.task_id, self.input_data.control_id)
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the web_api func. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/lnc_rna/tf_network "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="lnc_rna",
            task_type="2",  # maybe you need to change it
            submit_location="tfnetwork",
            tfbs_predict_id="5c76357b17b2bf3ec4f2f71e",
            thresh="0.0001",
            qv_thresh="0",
            corr_pvalue="0.05",
            corr_cutoff="0.5",
            corr_use_pdajust="0",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
