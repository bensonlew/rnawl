# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""otu表的样本距离计算"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class BatchWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BatchWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # {'name': 'count', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'batch_method', 'type': 'string'},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'count_batch', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'ellipse', 'type': 'string'},
            {'name': 'pca', 'type':'string', 'default': 'no'},
            {'name': 'main_id', 'type': 'string'}

            # {'name': 'group_dict', 'type': 'string'}

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        if self.option('has_batch') == True:
            self.run_batch()
        else:
            self.run_SVA()
        super(BatchWorkflow, self).run()

    def run_batch(self):
        self.batch = self.add_tool("tool_lab.batch_effect.batch")
        self.batch.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'batch_matrix': self.option('batch_matrix'),
            'batch_method': self.option('batch_method')
        })
        if self.option('pca') == 'yes':
            self.batch.on('end', self.run_exp_pca)
        if self.option('pca') == 'no':
            self.batch.on('end', self.set_output)
        self.batch.run()



    def run_SVA(self):
        self.SVA = self.add_tool("tool_lab.batch_effect.sva")
        self.SVA.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
        })
        if self.option('pca') == 'yes':
            self.SVA.on('end', self.run_exp_pca)
        if self.option('pca') == 'no':
            self.SVA.on('end', self.set_output)
        self.SVA.run()

    def run_exp_pca(self):
        self.exp_batch_pca = self.add_tool("ref_rna_v2.exp_pca")
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
            self.exp_batch_pca.on('end', self.set_output)
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
        self.exp_batch_corr = self.add_tool("ref_rna_v2.exp_corr")
        if self.option('has_batch') == True:
            self.exp_batch_corr.set_options({
                'exp': self.batch.option('count_batch').prop['path']
            })
        else:
            self.exp_batch_corr.set_options({
                'exp': self.SVA.option('count_batch').prop['path']
            })
        self.exp_batch_corr.on('end', self.set_output)
        self.exp_batch_corr.run()

    def set_output(self):
        if self.option('has_batch') == True:
            p = self.batch.option('count_batch').path
        else:
            p = self.SVA.option('count_batch').path
        link_names = os.path.join(self.output_dir, os.path.basename(p))
        os.link(p, link_names)
        self.option('count_batch').set_path(link_names)
        self.end()
        # self.set_db()

    def set_db(self):
         """
         保存结果标准化数据到mongo数据库中
         """
         record_id = self.option('exp_id')
         exp_batch = self.api.api('ref_rna_v3.batch')
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
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "批次效应处理结果目录", 0],
        ])
        super(BatchWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.batch import BatchWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'batch_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.batch',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/unigene.tpm.matrix.annot.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/group_table.txt',
                'batch_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/batch_new.txt',
                'batch_method': 'combat',
                'has_batch': True
            }
        }
        wsheet = Sheet(data=data)
        wf =BatchWorkflow(wsheet)
        wf.sheet.id = 'batch_effect'
        wf.sheet.project_sn = 'batch_effect'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
