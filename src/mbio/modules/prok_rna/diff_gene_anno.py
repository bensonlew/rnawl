# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 13:52
@file    : diff_gene_anno.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import glob
import shutil
import unittest

# 临时使用
# from src.biocluster.module import Module



class DiffGeneAnnoModule(Module):
    def __init__(self, work_id):
        super(DiffGeneAnnoModule, self).__init__(work_id)
        options = [
            # 基因列表文件和注释文件
            dict(name="gene_list", type="infile", format="prok_rna.diff_gene_list"),
            dict(name="anno_matrix", type="infile", format="prok_rna.diff_anno_matrix"),

            # 输出设置
            dict(name="outdir", type="string", default=None),
            dict(name="outfile_name", type="string", default='diff_gene_anno'),

            dict(name='pool', type='int', default=1)
        ]
        self.tool_list = []
        self.add_option(options)

    def check_options(self):
        if not self.option('gene_list'):
            raise OptionError('必须输入gene_list文件')
        if not self.option('anno_matrix'):
            raise OptionError('必须输入anno_matrix文件')
        if not self.option('outdir'):
            raise OptionError('设置tool输出文件夹')
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()
        self.end()

    def anno_run(self):
        self.anno_tool = self.add_tool('prok_rna.diff_gene_anno')
        options = {
            "gene_list": self.option('gene_list'),
            "anno_matrix": self.option('anno_matrix'),
            "outdir": self.option('outdir'),
            "outfile_name": self.option('outfile_name'),
        }
        self.step.add_steps('diff_gene_anno')
        step = getattr(self.step, 'diff_gene_anno')
        step.start()
        self.anno_tool.set_options(options)
        self.anno_tool.on("end", self.finish_update, "diff_gene_anno")
        self.anno_tool.run()
        self.tool_list.append(self.anno_tool)

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        """
        运行
        :return:
        """
        super(DiffGeneAnnoModule, self).run()
        self.anno_run()


    def end(self):
        # result_dir = self.add_upload_dir(self.anno_tool.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "差异分析结果目录"],
        # ])
        super(DiffGeneAnnoModule, self).end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/diff_gene'
            data = {
                "id": "diff_gene_extract_" + str(random.randint(1, 10000)),
                "type": "module",
                "name": "prok_rna.diff_gene_anno",
                "instant": False,
                "options": {
                    'gene_list': r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/gene.list',
                    'anno_matrix': r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/all_anno_detail.xls',
                    # 'task_id': 'tsg_32038',
                    'pool': 1,
                    'outdir': r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng',
                    'outfile_name': 'diff_gene_annotation.xls'
                }
            }

            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
