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


class PositionDepthWorkflow(Workflow):
    """
    read depth at each position or region
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PositionDepthWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "in_bed", "type": "infile", "format": "ref_rna_v2.bed"},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.position_depth")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(PositionDepthWorkflow, self).run()

    def check_options(self):
        if not self.option("in_bam").is_set:
            raise OptionError("必须设置输入BAM文件")
        if not self.option("in_bed").is_set:
            raise OptionError("必须设置输入BED文件")
        return True

    def run_tool(self):
        opts = {
            "in_bed": self.option('in_bed'),
            "in_bam": self.option('in_bam'),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        depth_api = self.api.api('tool_lab.position_depth')
        out_table = glob.glob(os.path.join(self.tool.output_dir, '*_depth.xls'))
        print(out_table)
        depth_api.add_position_depth(result=out_table, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "位点深度统计结果文件",0],
            ['./*_depth.xls', 'xls', '位点深度统计结果文件', 0],
        ])
        super(PositionDepthWorkflow, self).end()


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
            in_bed='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/position_depth/hs_cdk1.bed',
            in_bam='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/position_depth/E_3.bam'
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.position_depth",
            main_table_name="sg_position_depth",
            task_id="position_depth",
            project_sn="position_depth",
            submit_location="position_depth"
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