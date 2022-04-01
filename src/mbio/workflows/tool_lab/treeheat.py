# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210831'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class TreeheatWorkflow(Workflow):
    """
    Treeheat base R package Treeheat
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TreeheatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "target_lab", "type": "string", 'default': ''},
            {"name": "show_apart", "type": "string", 'default': "heat-tree"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.Treeheat = self.add_tool("tool_lab.treeheat")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_Treeheat()
        super(TreeheatWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set:
            raise OptionError("必须设置输入表达量文件")
        return True
        if not self.option("target_lab"):
            raise OptionError("必须设置目标列名")
        return True

    def run_Treeheat(self):
        opts = {
            "infile":self.option('infile').prop['path'],
            "target_lab":self.option('target_lab'),
            "show_apart":self.option('show_apart')
        }
        self.Treeheat.set_options(opts)
        self.Treeheat.on('end', self.set_db)
        self.Treeheat.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        if os.path.exists(os.path.join(self.Treeheat.output_dir, 'treeheat.pdf')):
            pdf_s3_path = os.path.join(self._sheet.output,"treeheat.pdf")
        else:
            pdf_s3_path = None
        mongo_api = self.api.api('tool_lab.treeheat')
        mongo_api.add_treeheat(
            s3_output=pdf_s3_path,
            main_id = self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.Treeheat.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "热图决策树",0],
            ['./treeheat.pdf', '', '热图决策树', 0],
        ])
        super(TreeheatWorkflow, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitoollabtest.py '
        cmd += 'post toollabpipeline '
        cmd += '-c {} '.format("client03")
        cmd += "-b http://bcl.tsg.com "
        cmd += "-n \"params;basis\" -d \"{"
        args = dict(
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/bsa/pop.final.vcf',
            # wp='',
            # mp='HQS1',
            mb='XS11_1',
            wb='F44_mix',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.bsa",
            main_table_name="sg_bsa",
            task_id="bsa",
            project_sn="bsa",
            submit_location="bsa"
        )
        for arg in args:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += args[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "};{"
        for arg in config:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += config[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "}\""

        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()