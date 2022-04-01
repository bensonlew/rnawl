# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 15:33
@file    : diff_gene_anno.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import glob
import os
import json
import time
import unittest

import pandas as pd

from biocluster.workflow import Workflow
from biocluster.file import getsize, exists
from biocluster.file import download

# from src.biocluster.workflow import Workflow
from bson import ObjectId


class DiffGeneAnnoWorkflow(Workflow):

    def __init__(self, wsheet_object):
        # 初始化网端参数
        self._sheet = wsheet_object
        super(DiffGeneAnnoWorkflow, self).__init__(wsheet_object)

        options = [
            # 基因列表文件和注释文件
            dict(name="gene_list", type="infile", format="prok_rna.diff_gene_list"),
            dict(name="anno_matrix", type="infile", format="prok_rna.diff_anno_matrix"),
            dict(name='task_id', type='string', default='tsg_32038'),
            dict(name='diff_main_id', type='string'),
            dict(name='update_info', type='string'),
        ]
        # 获取参数
        self.add_option(options)
        self.set_options(self._sheet.options())

        # 输出设置
        self.filepath = os.path.join(self.output_dir, 'diff_gene_annotation.xls')
        # self.option('outdir', outdir)
        # self.option('outfile_name', os.path.basename(self.filepath))

        self.module = self.add_module("prok_rna.diff_gene_anno")
        self.db_tool = self.api.api("prok_rna.all_exp")

    def run(self):
        self.module.on("end", self.set_db)
        self.run_module()
        super(DiffGeneAnnoWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        def diff_gene_anno_to_db(
            self, outpath, task_id=None, main_id=None, query_dict: dict = None,
            project_sn='prok_rna', main_table_name='diff_geen_anno'):
        """
        # add result info
        self.db_tool.diff_gene_anno_to_db(
            self.filepath, task_id=self.option('task_id'), main_id=self.option('diff_main_id'),
            main_table_name='diff_gene_anno_extr'
        )
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.tool.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "差异分析结果目录"],
        # ])
        super(DiffGeneAnnoWorkflow, self).end()

    def run_module(self):
        options = {
            'gene_list': self.option('gene_list'),
            'anno_matrix': self.option('anno_matrix'),
            'outdir': self.output_dir,
            'outfile_name': os.path.basename(self.filepath),
            'pool': 1
        }
        self.module.set_options(options)
        self.module.run()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random

    main_id = str(ObjectId("5bd907a8a4e1af0a8255a566"))

    data = {
        "id": "diff_gene_extract_" + str(random.randint(1, 10000)),
        "type": "workflow",
        "name": "prok_rna.diff_gene_anno",
        "instant": False,
        "options": {
            'gene_list': r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/gene.list',
            'anno_matrix': r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/all_anno_detail.xls',
            'task_id': 'tsg_32038',
            'diff_main_id': str(main_id),
            'update_info': json.dumps({'main_id': str(main_id)})
        }
    }

    wsheet = Sheet(data=data)
    wf = DiffGeneAnnoWorkflow(wsheet)
    wf.run()
