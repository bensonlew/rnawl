# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 10:53
@file    : diff_gene_anno.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import glob
import os

import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
# # 临时添加，写好后删除
# from src.biocluster.agent import Agent
# from src.biocluster.tool import Tool


class DiffGeneAnnoAgent(Agent):

    def __init__(self, parent):
        super(DiffGeneAnnoAgent, self).__init__(parent)
        options = [
            # 基因列表文件和注释文件
            dict(name="gene_list", type="infile", format="prok_rna.diff_gene_list"),
            dict(name="anno_matrix", type="infile", format="prok_rna.diff_anno_matrix"),

            # 输出设置
            dict(name="outdir", type="string", default=None),
            dict(name="outfile_name", type="string", default='diff_gene_anno'),

            dict(name='pool', type='int', default=1),

        ]
        self.add_option(options)

    def set_resource(self):
        self._cpu = self.option("pool") + 2
        self._memory = "{}G".format(self.option("pool")*20)

    def end(self):
        out_dir = self.add_upload_dir(self.output_dir)
        out_dir.add_relpath_rules([
            [".", "", "差异基因注释文件输出目录"]
        ])
        super(DiffGeneAnnoAgent, self).end()


class DiffGeneAnnoTool(Tool):
    """
        Differential analysis based on edgeR, DEGseq, DEseq2.
    """

    def __init__(self, config):
        super(DiffGeneAnnoTool, self).__init__(config)
        self.python_path = os.path.join(
            self.config.SOFTWARE_DIR, 'program', 'Python', 'bin', 'python')
        self.package_path = os.path.join(
            self.config.PACKAGE_DIR, 'prok_rna', 'gene_anno_extr.py')

    def run(self):
        super(DiffGeneAnnoTool, self).run()

        # 计算
        self.run_package()
        self.set_output()
        self.end()

    def run_package(self):
        cmd = ' '.join([
            self.python_path, self.package_path,
            '-g', self.option('gene_list').prop['path'],
            '-a', self.option('anno_matrix').prop['path'],
            '-f', self.option('outfile_name'),
            '-o', self.option('outdir')
        ] )
        cmd_name = 'diff_gene_annotation_package'
        # 运行command并检测错误
        for i in range(2):
            command = self.add_command(cmd_name, cmd)
            command.software_dir = ''
            if i == 0:
                command.run()
            else:
                command.rerun()

            self.wait()

            if command.return_code == 0:
                self.logger.info(cmd_name + ' 运行成功, 命令: ' + cmd)
                break
            elif command.return_code != 0:
                if i == 2:
                    self.set_error(cmd_name + ' command 运行错误, 命令: ' + cmd)
                else:
                    self.logger.info(cmd_name + ' 运行出错, 开始重新尝试' + cmd)

    def set_output(self):
        file = os.path.join(self.option("outdir"), self.option('outfile_name'))
        link = os.path.join(self.output_dir, self.option('outfile_name'))
        if os.path.exists(link):
            os.remove(link)
        os.link(file, link)


if __name__ == '__main__':

    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """
        test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "DiffExpRxtr" + str(random.randint(1, 10000)),
                "type": "tool",
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
