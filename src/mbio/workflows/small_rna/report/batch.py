# -*- coding: utf-8 -*-
import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class BatchWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BatchWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            # {'name': 'count', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'batch_method', 'type': 'string'},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'exp_id', 'type': 'string'},
            {'name': 'sg_exp_batch_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'count_batch', 'type': 'outfile', 'format': 'denovo_rna_v2.common'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'exp_other', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'params', 'type': 'string'},
            {'name': 'ellipse', 'type': 'string'}
            # {'name': 'group_dict', 'type': 'string'}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Express/04 Exp_Batch')
        self.inter_dirs = []

        # self.task_id = self._sheet.id
    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(BatchWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        if self.option('has_batch') == True:
            self.run_batch()
        else:
            self.run_SVA()
        super(BatchWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_exp_batch", main_id=self.option('sg_exp_batch_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_batch(self):
        self.batch = self.add_tool('denovo_rna_v3.batch.batch')
        self.batch.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'batch_matrix': self.option('batch_matrix'),
            'batch_method': self.option('batch_method')
        })
        self.batch.on('end',self.run_exp_pca)
        self.batch.run()



    def run_SVA(self):
        self.SVA = self.add_tool('denovo_rna_v3.batch.sva')
        self.SVA.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
        })

        self.SVA.on('end', self.run_exp_pca)
        self.SVA.run()

    def run_exp_pca(self):
        self.exp_batch_pca = self.add_tool("denovo_rna_v2.exp_pca")
        if self.option('has_batch') == True:
            self.exp_batch_pca.set_options({
                'exp': self.batch.option('count_batch').prop['path']
            })
        else:
            self.exp_batch_pca.set_options({
                'exp': self.SVA.option('count_batch').prop['path']
            })
        if self.option('ellipse') == 'yes':
            self.exp_batch_pca.on('end', self.run_ellipse)
        if self.option('ellipse') == 'no':
            self.exp_batch_pca.on('end', self.run_exp_corr)
        self.exp_batch_pca.run()

    def run_ellipse(self):
        self.ellipse = self.add_tool('graph.ellipse')
        self.ellipse.set_options({
            'analysis': 'pca',
            'group_table': self.option('group_table').prop['path'],
            'pc_table': os.path.join(self.exp_batch_pca.output_dir, 'PCA.xls'),
        })
        self.ellipse.on('end', self.run_exp_corr)
        self.ellipse.run()

    def run_exp_corr(self):
        self.exp_batch_corr = self.add_tool("small_rna.exp_corr")
        if self.option('has_batch') == True:
            self.exp_batch_corr.set_options({
                'express_matrix': self.batch.option('count_batch').prop['path']
            })
        else:
            self.exp_batch_corr.set_options({
                'express_matrix': self.SVA.option('count_batch').prop['path']
            })
        self.exp_batch_corr.on('end', self.set_output)
        self.exp_batch_corr.run()

    def set_output(self):
        if self.option('has_batch') == True:
            p = self.batch.option('count_batch').path
        else:
            p = self.SVA.option('count_batch').path
        link_names = os.path.join(self.output_dir, os.path.basename(p))
        if os.path.exists(link_names):
            os.remove(link_names)
        os.link(p, link_names)
        self.option('count_batch').set_path(link_names)
        self.set_db()

    def set_db(self):
         """
         保存结果标准化数据到mongo数据库中
         """
         record_id = self.option('exp_id')
         exp_batch = self.api.api('small_rna_v2.batch')
         exp_batch.add_exp_batch(self.option('count_batch').path, self.option('exp_other').path, record_id,
                                           self.option('sg_exp_batch_id'),self.option('task_id'), self.option('params'))




         # all_exp = self.api.api("ref_rna_v2.all_exp")
         # add result info
         exp_batch.add_exp_batch_pca(self.exp_batch_pca.output_dir, main_id=self.option('sg_exp_batch_id'), record_id=record_id )
         if self.option('ellipse') == 'yes':
            exp_batch.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', main_id=self.option('sg_exp_batch_id'), record_id=record_id)
         if self.option('ellipse') == 'no':
             pass
         exp_batch.add_exp_batch_corr(self.exp_batch_corr.work_dir, main_id=self.option('sg_exp_batch_id'), record_id=record_id)
         self.end()


    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["01 Express", "", "表达量结果目录",0],
            ["01 Express/04 Exp_Batch", "", "批次效应评估", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "批次效应评估文件", 0],
            ["count_batch.txt", "txt", "批次效应处理表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(BatchWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_kit_standard import TableKitStandardWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'standard_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_kit_standard',
            'options': {
                'sep' : 'tab',
                'observe' : 'sum',
                'feature' : 'standard_scale',
                'table' : '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/table.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf =TableKitStandardWorkflow(wsheet)
        wf.sheet.id = 'table_kit_standard'
        wf.sheet.project_sn = 'table_kit_standard'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
