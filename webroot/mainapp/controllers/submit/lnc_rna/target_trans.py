# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.controllers.project.lnc_rna_controller import LncRnaController

class TargetTransAction(LncRnaController):
    """
    lnc_rna trans靶基因预测接口
    """
    def __init__(self):
        super(TargetTransAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        ret = self.check_params()
        if ret is not True:
            return ret
        self.pack_params()
        self.create_main_table()
        self.prepare_workflow()
        task_info = super(TargetTransAction, self).POST()
        task_info['content'] = {'ids': {'id': str(self.main_id), 'name': self.name}}
        return json.dumps(task_info)

    def check_params(self):
        self.input_data = web.input()
        self.expected_args = ['lnc_geneset_id', 'm_geneset_id', 'task_type', 'method', 'submit_location', 'task_id']
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        task_info = self.lnc_rna.get_task_info(self.input_data.task_id)
        self.lncset_dict = self.lnc_rna.get_main_info(self.input_data.lnc_geneset_id, "sg_geneset", self.input_data.task_id)
        self.mset_dict = self.lnc_rna.get_main_info(self.input_data.m_geneset_id, "sg_geneset", self.input_data.task_id)
        if not task_info:
            info = {'success': False, 'info': '无法找到任务信息: {}'.format(task_id)}
            return json.dumps(info)
        if self.lncset_dict["gene_length"] * self.mset_dict["gene_length"] > 10000:
            if self.input_data.method in ["intarna", "lnctar"]:
                info = {'success': False, 'info': '基因数量过多，该方法运行过慢，请缩减基因集'}
                return json.dumps(info)
            else:
                pass
        self.task_dict = task_info
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == 'task_type':
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        self.params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self):
        result_info = self.lnc_rna.get_task_info(self.input_data.task_id)

        self.project_sn = result_info['project_sn']
        time_now = datetime.datetime.now()
        self.name = 'Target_trans_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Geneset go Dag main table built at controller'
        main_info = {
            'task_id': self.input_data.task_id,
            'project_sn': self.project_sn,
            'created_ts': created_ts,
            'name': self.name,
            'desc': desc,
            'params': self.params,
            'status': 'start',
        }
        self.main_id = self.lnc_rna.insert_main_table('sg_target_trans', main_info)

    def prepare_workflow(self):
        input_data_dict = dict(self.input_data)
        # main_info = self.lnc_rna.get_main_info(self.input_data.go_enrich_id, 'sg_geneset_go_enrich', self.input_data.task_id)
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

        latest_info_target = self.lnc_rna.get_main_info_by_record("sg_target", task_id=str(self.input_data.task_id), status="end",  type="latest")
        if latest_info_target:
            latest_main_id_target = latest_info_target["_id"]
        else:
            latest_main_id_target = ""
        update_info = {str(self.main_id): "sg_target_trans"}

        options.update({
            "target_trans_id": str(self.main_id),
            "known": self.task_dict["known_lnc_fa"],
            "novol": self.task_dict["novel_lnc_fa"],
            "mrna": self.task_dict["assemble_fa"],
            "novol_gtf": self.task_dict["novel_lnc_gtf"],
            "known_gtf": self.task_dict["known_lnc_gtf"],
            "m_geneset_type": self.mset_dict['type'],
            "lnc_geneset_type": self.lncset_dict['type'],
            "mrna_gtf": self.task_dict["mrna_gtf"],
            "last_id_target": str(latest_main_id_target),
            'update_info': json.dumps(update_info),
            "annotation": annot_file
        })
        to_file = ['lnc_rna.export_geneset_list(m_geneset_id)', 'lnc_rna.export_geneset_list(lnc_geneset_id)']
        self.set_sheet_data(name='lnc_rna.report.target_trans',
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
        cmd += 's/lnc_rna/target_trans '
        cmd += '-b http://192.168.12.101:9090 '
        args = {
            'task_id': 'lnc_rna',
            'submit_location': 'target_trans',
            'task_type': '2',
            'm_geneset_id': '5c90894a17b2bf407167bc7f',
            'lnc_geneset_id': '5c90890017b2bf407167bc7d',
            'method': 'rnaplex'
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
