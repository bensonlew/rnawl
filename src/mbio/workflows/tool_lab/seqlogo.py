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


class SeqlogoWorkflow(Workflow):
    """
    (reverse-/)transcription between DNA and RNA
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqlogoWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'file_type', 'type': 'string', 'default': 'fasta'},
            {"name": "infile", "type": "infile", 'format': 'ref_rna_v2.common'},
            {'name': 'seq_type', 'type': 'string', 'default': 'nt'},
            {'name': 'matrix_type', 'type': 'string', 'default': 'none'},
            {'name': 'unit', 'type': 'string', 'default': 'both'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seqlogo")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(SeqlogoWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set:
            raise OptionError("必须设置输入数据文件")
        return True

    def run_tool(self):
        opts = {
            'infile': self.option('infile'),
            'seq_type': self.option('seq_type'),
        }
        if self.option('file_type') == 'fasta':
            m_type = 'none'
        elif self.option('matrix_type') in ['pfm-4-columns', 'pfm-4-rows']:
            m_type = self.option('matrix_type').replace('4', 'four')
        else:
            m_type = self.option('matrix_type')
        unit_n = self.option('unit').split(',')
        if len(unit_n) > 1:
            unit = 'both'
        else:
            unit = self.option('unit')
        opts.update({'matrix_type': m_type, 'unit': unit})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "SeqLogo图结果文件",0],
            [r'*.pdf', 'pdf', 'SeqLogo图结果文件', 0],
        ])
        super(SeqlogoWorkflow, self).end()


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
            file_type='matrix',
            infile='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_29032021/test_4rows.txt',
            seq_type='nt',
            matrix_type='pfm-4-rows',
            unit='bits,probability',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.seqlogo",
            main_table_name="seqlogo",
            task_id="seqlogo",
            project_sn="seqlogo",
            submit_location="seqlogo"
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