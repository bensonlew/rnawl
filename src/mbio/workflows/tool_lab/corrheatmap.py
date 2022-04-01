# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210729'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class CorrheatmapWorkflow(Workflow):
    """
    Corrheatmap base R package Corrheatmap
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CorrheatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "type", "type": "string", 'default': 'full'},
            {"name": "method", "type": "string", 'default': "circle"},
            {"name": "color", "type": "string", 'default': "Spectral"},
            {"name": "corr_method", "type": "string", 'default': "pearson"},
            {"name": "corr_adjust", "type": "string", 'default': 'holm'},
            {"name": "show_significance", "type": "string", 'default': 'yes'},
            {"name": "in_significance", "type": "string", 'default': 'pch'},
            {"name": "show_coefficient", "type": "string", 'default': 'no'},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.Corrheatmap = self.add_tool("tool_lab.corrheatmap")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_Corrheatmap()
        super(CorrheatmapWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set:
            raise OptionError("必须设置输入表达量文件")
        return True

    def run_Corrheatmap(self):
        opts = {
            "infile":self.option('infile').prop['path'],
            "type":self.option('type'),
            "method":self.option('method'),
            "color":self.option('color'),
            "corr_method":self.option('corr_method'),
            "corr_adjust":self.option('corr_adjust'),
            "show_significance":self.option('show_significance'),
            "in_significance":self.option('in_significance'),
            "show_coefficient":self.option('show_coefficient')
        }
        self.Corrheatmap.set_options(opts)
        self.Corrheatmap.on('end', self.set_db)
        self.Corrheatmap.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        if os.path.exists(os.path.join(self.Corrheatmap.output_dir, 'corr_heatmap.pdf')):
            pdf_s3_path = os.path.join(self._sheet.output,"corr_heatmap.pdf")
        else:
            pdf_s3_path = None
        Circularbar_api = self.api.api('tool_lab.corrheatmap')
        Circularbar_api.add_corrheatmap(
            s3_output=pdf_s3_path,
            main_id = self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.Corrheatmap.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性热图2.0",0],
            ['./corr_heatmap.pdf', '', '相关性热图2', 0],
        ])
        super(CorrheatmapWorkflow, self).end()


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