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


class ReadsRepairWorkflow(Workflow):
    """
    reads repair using bbmap
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ReadsRepairWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "r1", "type": "infile", "format": "ref_rna_v2.fastq"},
            {'name': 'r2', "type": 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.tool = self.add_tool("tool_lab.reads_repair")
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ReadsRepairWorkflow, self).run()

    def check_options(self):
        if not self.option("r1").is_set:
            raise OptionError("请上传R1序列")
        if not self.option("r2").is_set:
            raise OptionError("请上传R2序列")
        return True

    def run_tool(self):
        opts = {
            'r1': self.option('r1'),
            'r2': self.option('r2'),
        }
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
            [".", "", "Reads数据修复结果目录",0],
            ['*', '', 'Reads数据修复结果文件', 0],
        ])
        super(ReadsRepairWorkflow, self).end()


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
            r1='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/reads_repair/example.R1.fq.gz',
            r2='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/reads_repair/example.R2.fq.gz',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.reads_repair",
            main_table_name="sg_reads_repair",
            task_id="reads_repair",
            project_sn="reads_repair",
            submit_location="reads_repair"
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