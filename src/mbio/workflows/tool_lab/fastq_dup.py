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


class FastqDupWorkflow(Workflow):
    """
    (reverse-/)transcription between DNA and RNA
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FastqDupWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.fastq_dup = self.add_module("tool_lab.fastq_dup")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_module()
        super(FastqDupWorkflow, self).run()

    def check_options(self):
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须设置输入fq序列文件夹")
        return True

    def run_module(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
        }
        self.fastq_dup.set_options(opts)
        self.fastq_dup.on('end', self.set_db)
        self.fastq_dup.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dup_api = self.api.api('tool_lab.fastq_dup')
        dup_data = os.path.join(self.fastq_dup.output_dir, 'dup.xls')
        dup_api.add_fastq_dup(data=dup_data, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.fastq_dup.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "序列冗余度统计结果文件",0],
            ['./dup.xls', 'xls', '序列冗余度统计表', 0],
        ])
        super(FastqDupWorkflow, self).end()


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
            fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_042021/fastq_dup_2",
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.fastq_dup",
            main_table_name="sg_fastq_dup",
            task_id="fastq_dup",
            project_sn="fastq_dup",
            submit_location="fastq_dup"
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