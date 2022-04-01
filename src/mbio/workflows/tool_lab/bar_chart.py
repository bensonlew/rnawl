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


class BarChartWorkflow(Workflow):
    """
    (reverse-/)transcription between DNA and RNA
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BarChartWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "group_file", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'set_group', 'type': 'string', 'default': ''},
            {'name': 'ishape', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.bar_chart")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(BarChartWorkflow, self).run()

    def check_options(self):
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if not self.option("set_group"):
            raise OptionError("必须设置是否分组")
        return True

    def run_tool(self):
        opts = {
            'data_file': self.option('data_file'),
        }
        if self.option('group_file').is_set:
            opts.update({'group_file': self.option('group_file')})
        if self.option('ishape'):
            opts.update({'ishape': self.option('ishape')})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        bar_api = self.api.api('tool_lab.bar_chart')
        data4bar = self.tool.option('out_file').prop['path']
        if self.option('group_file').is_set:
            bar_api.add_bar_detail(data=data4bar, main_id=self.option('main_id'), group_file=self.option('group_file').prop['path'])
        else:
            bar_api.add_bar_detail(data=data4bar, main_id=self.option('main_id'), data_file=self.option('data_file').prop['path'])
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.tool.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "柱状图结果文件",0],
        #     [r'*.txt', 'txt', '柱状图结果文件', 0],
        # ])
        super(BarChartWorkflow, self).end()


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
            data_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_29032021/data_test.txt',
            # group_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_29032021/group_test.txt',
            ishape='none',
            set_group='no',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.bar_chart",
            main_table_name="bar_chart",
            task_id="bar_chart",
            project_sn="bar_chart",
            submit_location="bar_chart"
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