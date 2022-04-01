# -*- coding: utf-8 -*-
"""
@time    : 2018/10/26 9:25
@file    : diff_gene_extr.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import os
import unittest

# from webroot.mainapp.controllers.project.prok_rna_controller import ProkRNAController

class DiffGeneExtrAction(ProkRNAController):
    def __init__(self):
        super(DiffGeneExtrAction, self).__init__(instant=False)
        self.web_data = web.input()
        self.task_id = self.web_data.task_id
        self.workflow_path = 'prok_rna.report.diff_gene_anno'
        self.expected_args = ['task_id']

    @check_sig
    def POST(self):
        # 参数检测
        check_stat = self.check_params(self.web_data, self.expected_args)
        if check_stat is not True:
            return check_stat
        # 创建主表
        main_table_name = 'diff_gene_anno_extr'
        main_id, task_name, project_sn = self.create_main_table(main_table_name)
        # 准备参数
        self.prepare_workflow_params(project_sn, main_id, main_table_name)
        # 运行workflow 并传回参数
        task_info = super(DiffGeneExtrAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': task_name
            }
        }
        # # task_info['group_dict'] = group_dict
        # if 'group_id' in self.web_data and str(self.web_data.group_id).lower() != 'all':
        #     _ = self.prok_rna.update_group_is_use(self.web_data.task_id, self.web_data.group_id)
        # if 'control_id' in self.web_data:
        #     _ = self.prok_rna.update_group_compare_is_use(self.task_id, self.web_data.control_id)
        return json.dumps(task_info)

    def check_params(self, data, args):
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        return True

    def create_main_table(self, collection_name):
        task_id = self.task_id
        info = self.prok_rna.get_table_info_by_task_id(task_id, table_name='sg_diff')
        print(info)
        project_sn = info['project_sn']
        exp_level = info['exp_level']
        now_time = datetime.datetime.now()
        task_name = '_'.join(('DiffGeneExtract', exp_level, now_time.strftime('%Y%m%d_%H%M%S')))
        main_info = dict(
            name=task_name,
            task_id=task_id,
            project_sn=project_sn,
            desc='differentially expressed gene annotation',
            status="start",
            created_ts=now_time.strftime('%Y-%m-%d %H:%M:%S'),
            params=json.dumps({  # json format
                'task_id': self.task_id,
                'exp_level': exp_level,
                'submit_location': self.web_data.submit_location
            })
        )
        # 创建主表
        main_id = self.prok_rna.insert_main_table(collection_name, main_info)
        return main_id, task_name, project_sn

    def prepare_workflow_params(self, project_sn, main_id, main_table_name):
        options = {
            'gene_list': {'task_id': self.web_data['task_id'], 'main_table_name': 'sg_diff'},
            'anno_matrix': {'task_id': self.web_data['task_id'], "type" : "origin"},
            'task_id': self.task_id,
            'diff_main_id': str(main_id),
            'update_info': json.dumps({'main_id': str(main_id)})
        }

        to_file = [
            'prok_rna.export_diff_genes(gene_list)',
            'prok_rna.get_anno_file_path(anno_matrix)'
        ]

        self.set_sheet_data(
            name=self.workflow_path,
            options=options,
            main_table_name=main_table_name,
            task_id=self.task_id,
            project_sn=project_sn,
            module_type='workflow',
            to_file=to_file
        )


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            # cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/prok_rna/diff_gene_extr "
            cmd += "-b http://192.168.12.102:9090 "
            args = dict(
                task_id="tsg_32038",
                task_type="2",
                submit_location='diff_gene_extr'
            )
            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
            print(cmd)
            print('=============================================================' * 2)
            os.system(cmd)

    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
