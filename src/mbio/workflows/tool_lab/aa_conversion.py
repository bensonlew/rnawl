# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class AaConversionWorkflow(Workflow):
    """
    3-letter to 1-letter/ 1-letter to 3-letter amino acid conversion
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AaConversionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'seq', 'type': 'string'},
            {'name': 'type', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.conversion = self.add_tool("tool_lab.aa_conversion")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_conversion()
        super(AaConversionWorkflow, self).run()

    def check_options(self):
        if not self.option("type"):
            raise OptionError("必须设置输入类型选择")
        if not self.option('seq'):
            raise OptionError("必须设置输入原始序列")
        return True

    def run_conversion(self):
        opts = {
            'aa_str': self.option('seq'),
            'aa_type': self.option('type')
        }
        self.conversion.set_options(opts)
        self.conversion.on('end', self.set_db)
        self.conversion.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        aa_api = self.api.api('tool_lab.aa_conversion')
        result = os.path.join(self.conversion.output_dir, 'result.fa')
        aa_api.add_aa_main(result=result, main_id=self.option('main_id'))
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.conversion.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基于VCF进行BSA分析",0],
        #     [r'./pop\.result\..*', '', '关联区域变异位点统计表', 0],
        #     [r'.*\.png|pdf', '', 'Manhattan图', 0],
        # ])
        super(AaConversionWorkflow, self).end()


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
            seq='AlaCysAspGluPheGlyHisIleLysAsxXaaGlx',
            type='three',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.aa_conversion",
            main_table_name="sg_aa_conversion",
            task_id="aa_conversion",
            project_sn="aa_conversion",
            submit_location="aa_conversion"
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