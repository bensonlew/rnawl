# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""batch"""

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
from mbio.packages.prok_rna.chart import Chart


class BatchWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BatchWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'prok_rna.common'},
            {'name': 'count', 'type': 'infile', 'format': 'prok_rna.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'prok_rna.common'},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'prok_rna.common'},
            {'name': 'batch_method', 'type': 'string'},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'exp_id', 'type': 'string'},
            {'name': 'sg_exp_batch_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'count_batch', 'type': 'outfile', 'format': 'prok_rna.common'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'exp_other', 'type': 'infile', 'format': 'prok_rna.common'},
            {'name': 'params', 'type': 'string'},
            {'name': 'ellipse', 'type': 'string'}
            # {'name': 'group_dict', 'type': 'string'}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Express/04 Exp_Batch')
        # self.inter_dirs = []

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
        # self.get_run_log()
        if self.option('has_batch') == True:
            self.run_batch()
        else:
            self.run_SVA()
        super(BatchWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("prok_rna", table="sg_exp_batch", main_id=self.option('sg_exp_batch_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_batch(self):
        self.batch = self.add_tool('prok_rna.batch.batch')
        self.batch.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'batch_matrix': self.option('batch_matrix'),
            'batch_method': self.option('batch_method')
        })
        self.batch.on('end',self.run_exp_pca)
        self.batch.run()



    def run_SVA(self):
        self.SVA = self.add_tool('prok_rna.batch.sva')
        self.SVA.set_options({
            'count_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
        })

        self.SVA.on('end', self.run_exp_pca)
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
        self.set_db()

    def set_db(self):
         """
         保存结果标准化数据到mongo数据库中
         """
         record_id = self.option('exp_id')
         exp_batch = self.api.api('prok_rna.batch')
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

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'

        # PCA chart
        group_dict = dict()
        with open(self.option('group_table').prop['path'], 'r') as g:
            g.readline()
            for line in g:
                sample, group = line.strip().split('\t')
                if group not in group_dict.keys():
                    group_dict[group] = list()
                group_dict[group].append(sample)

        exp_pca_file = os.path.join(self.exp_batch_pca.output_dir, 'PCA.xls')
        exp_pca_var_file = os.path.join(self.exp_batch_pca.output_dir, 'Explained_variance_ratio.xls')
        if self.option('ellipse') == 'yes':
            exp_pca_ellipse = os.path.join(self.ellipse.output_dir, 'ellipse_out.xls')
        else:
            exp_pca_ellipse = None
        chart.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"])

        # Correlation Chart
        exp_corr_file = os.path.join(self.exp_batch_corr.work_dir, 'sample_correlation.xls')
        exp_corr_tree_file = os.path.join(self.exp_batch_corr.work_dir, 'sample.cluster_tree.txt')
        group_dict_json = json.dumps(group_dict)
        chart.chart_exp_corr(exp_corr_file, exp_corr_tree_file, group_dict_json)

        chart.to_pdf()

        # move pdf
        pca_pdf = os.path.join(self.work_dir, 'all.exp_relation_pca.scatter.pdf')
        new_path = os.path.join(self.output_dir, 'sample_batch_pca.pdf')
        self.move_pdf(pca_pdf, new_path)
        if os.path.exists(os.path.join(self.work_dir, 'all.exp_relation_pca_ell.scatter.pdf')):
            new_path = os.path.join(self.output_dir, 'sample_batch_pca_ell.pdf')
            self.move_pdf(os.path.join(self.work_dir, 'all.exp_relation_pca_ell.scatter.pdf'), new_path)

        self.move_pdf(os.path.join(self.work_dir, 'exp.heatmap.heat_corr.pdf'),
                      os.path.join(self.output_dir, 'sample_batch_correlation.pdf'))

    def end(self):
        self.chart()
        # if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
        #     os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        # os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        # self.inter_dirs = [
        #     ["01 Express", "", "表达量分析结果目录",0],
        #     ["01 Express/04 Exp_Batch", "", "批次效应评估", 0]
        # ]
        result_dir.add_relpath_rules([
            [".", "", "批次效应评估结果目录", 0],
            ["count_batch.txt", "txt", "批次效应处理表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['sample_batch_pca.pdf', 'pdf', '批次效应评估PCA图', 0],
            ['sample_batch_pca_ell.pdf', 'pdf', '批次效应评估PCA图(置信圈)', 0],
            ['sample_batch_correlation.pdf', 'pdf', '批次效应评估相关性热图图', 0],
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
