# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from collections import OrderedDict
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.controllers.project.lnc_rna_controller import LncRnaController



class TargetCistransAction(LncRnaController):
    """
    lnc_rna cis靶基因预测接口
    """
    def __init__(self):
        super(TargetCistransAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        ret = self.check_params()
        if ret is not True:
            return ret
        self.pack_params()
        self.create_main_table()
        self.prepare_workflow()
        task_info = super(TargetCistransAction, self).POST()
        task_info['content'] = {'ids': {'id': str(self.main_id), 'name': self.name}}
        return json.dumps(task_info)

    def check_params(self):
        self.input_data = web.input()
        print "input_data is "
        print self.input_data
        self.expected_args = ['up_dis', 'down_dis', 'task_type', 'submit_location', 'task_id']
        self.expected_args2 = ["submit_location", 'group_dict', 'group_id', 'corr_cutoff', 'corr_way', 'padjust_way', "pvalue_type", "pvalue_cutoff"]

        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        task_info = self.lnc_rna.get_task_info(self.input_data.task_id)
        # exp_g_dict =
        if not task_info:
            info = {'success': False, 'info': '无法找到任务信息: {}'.format(task_id)}
            return json.dumps(info)
        self.task_dict = task_info
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args + self.expected_args2:
            if each == 'task_type':
                params_dict[each] = int(input_data_dict[each])
            elif each == "group_dict":
                params_dict[each] = json.loads(self.input_data.group_dict, object_pairs_hook=OrderedDict)
            else:
                params_dict[each] = input_data_dict[each]
        self.params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))


    def create_main_table(self):
        result_info = self.lnc_rna.get_task_info(self.input_data.task_id)

        self.project_sn = result_info['project_sn']
        time_now = datetime.datetime.now()
        self.name = 'Target_cistrans_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Geneset go Dag main table built at controller'
        main_info = {
            'task_id': self.input_data.task_id,
            'project_sn': self.project_sn,
            'created_ts': created_ts,
            'name': self.name,
            "version": "v1.1",
            'desc': desc,
            'params': self.params,
            'status': 'start',
        }
        self.main_id = self.lnc_rna.insert_main_table('sg_target_cistrans', main_info)

    def prepare_workflow(self):
        input_data_dict = dict(self.input_data)
        group_dict = json.loads(self.input_data.group_dict, object_pairs_hook=OrderedDict)
        if str(self.input_data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])

        # main_info = self.lnc_rna.get_main_info(self.input_data.go_enrich_id, 'sg_geneset_go_enrich', self.input_data.task_id)
        exp_dict_G = self.lnc_rna.get_main_info_by_record("sg_exp", exp_level="G", task_id = input_data_dict['task_id'])
        exp_dict_T = self.lnc_rna.get_main_info_by_record("sg_exp", exp_level="T", task_id = input_data_dict['task_id'])

        options = dict()
        for each in self.expected_args:
            if input_data_dict[each] == '':
                options[each] = None
            else:
                options[each] = input_data_dict[each]

        annot_info = self.lnc_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(self.input_data.task_id), status="end",  type="origin")
        if annot_info and "result_dir" in annot_info:
            annot_file = os.path.join(annot_info["result_dir"], "allannot_class/all_annot.xls")
        else:
            annot_file = self.task_dict["annot"]

        update_info = {str(self.main_id): "sg_target_cistrans"}

        options.update({
            "target_cis_id": str(self.main_id),
            "novol": self.task_dict["novel_lnc_gtf"],
            "known": self.task_dict["known_lnc_gtf"],
            'update_info': json.dumps(update_info),
            "mrna_gtf": self.task_dict["mrna_gtf"],
            "group_id": input_data_dict["group_id"],
            "exp_matrix_lnc": str(exp_dict_T["main_id"])  + ";" + ";" + "All",
            "exp_matrix_target": str(exp_dict_G["main_id"])  + ";All" + ";",
            "annotation": annot_file,
            "group_dict": json.dumps(group_dict),
        })

        if input_data_dict['pvalue_type'] == "pvalue":
            options.update({
                "pvalue_cutoff": input_data_dict['pvalue_cutoff']
            })
        else:
            options.update({
                "qvalue_cutoff": input_data_dict['pvalue_cutoff']
            })
        options.update({
            "cor_cutoff": input_data_dict['corr_cutoff']
        })



        to_file = ["lnc_rna.export_geneset_lncset_exp_matrix(exp_matrix_lnc)",
                   "lnc_rna.export_geneset_lncset_exp_matrix(exp_matrix_target)"]
        self.set_sheet_data(name='lnc_rna.report.target_cistrans',
                            options=options,
                            main_table_name=self.name,
                            module_type='workflow',
                            to_file=to_file,
                            project_sn=self.project_sn,
                            task_id=self.input_data.task_id)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_default(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/lnc_rna/target_cistrans '
        cmd += '-b http://192.168.12.101:51232 '
        args = {
            'task_id': 'tsg_35567',
            'group_dict':  r'{"CL": ["CL1", "CL2", "CL5"], "HFL": ["HFL3", "HFL4", "HFL6"], "HGL": ["HGL1", "HGL3", "HGL4"]}'.replace(
                '"', '\\"'),
            'group_id':"5d7f031717b2bf538f78daac",
            'corr_cutoff': "0.8",
            'corr_way': "pearson",
            'padjust_way': "BH",
            "pvalue_type": "qvalue",
            "pvalue_cutoff": "0.05",
            'submit_location': 'target_cistrans',
            'task_type': '2',
            'up_dis': '20',
            'down_dis': '20',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
