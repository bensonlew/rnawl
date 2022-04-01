# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from biocluster.core.exceptions import OptionError


class SeqConversionWorkflow(Workflow):
    """
    (reverse-/)transcription between DNA and RNA
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqConversionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'type', 'type': 'string'},
            {"name": "infile", "type": "infile", 'format': 'ref_rna_v2.fasta'},     # file with miRNA in the first column
            {"name": "instr", "type": "string", 'default': ''},     # String of miRNA names/accessions. separate by comma
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_conversion")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(SeqConversionWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set and not self.option("instr"):
            raise OptionError("必须设置输入DNA/RNA序列")
        return True

    def run_tool(self):
        opts = {}
        if self.option('infile').is_set:
            opts.update({'infile': self.option('infile')})
        elif self.option('instr'):
            opts.update({'instr': self.option('instr')})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        # """
        # 保存结果表到mongo数据库中
        # """
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "DNA/RNA转换结果文件",0],
            [r'*.fa', 'FASTA', 'DNA/RNA转换结果文件', 0],
        ])
        super(SeqConversionWorkflow, self).end()


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
            infile='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_29032021/test_mix_seq.fa',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.sg_seq_conversion",
            main_table_name="seq_conversion",
            task_id="seq_conversion",
            project_sn="seq_conversion",
            submit_location="seq_conversion"
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

