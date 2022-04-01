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


class MirbaseConversionWorkflow(Workflow):
    """
    hypergeometric test function for gene set enrichment analysis that are designed to accept user defined annotation
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MirbaseConversionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'type', 'type': 'string'},
            {"name": "infile", "type": "infile", 'format': 'medical_transcriptome.common'},     # file with miRNA in the first column
            {"name": "instr", "type": "string", 'default': ''},     # String of miRNA names/accessions. separate by comma
            {"name": "target_ver", "type": "string", "default": "v22"},     # target miRBase version
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.mirbase_conversion")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(MirbaseConversionWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set and not self.option("instr"):
            raise OptionError("必须设置输入含有miRNA name/Accession的文件或字符串")
        if self.option("infile").is_set:
            file_form = os.path.basename(self.option("infile").prop['path']).split('.')[-1]
            if file_form not in ['txt', 'TXT']:
                raise OptionError("必须设置输入txt格式文件")
        if not self.option("target_ver"):
            raise OptionError("必须设置输入目标miRBase版本")
        return True

    def run_tool(self):
        opts = {
            'target_ver': self.option("target_ver"),
        }
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
        # mirbase_conversion = self.api.api("tool_lab.mirbase_conversion")
        # # add result info
        # mirbase_conversion.add_conversion(self.option('main_id'), result=self._sheet.output)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "miRBase ID转换结果文件",0],
            [r'*.xls', 'XLS', 'miRBase ID转换结果文件', 0],
        ])
        super(MirbaseConversionWorkflow, self).end()


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
            infile="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_190221/mirid_test.xls",
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.mirbase_conversion",
            main_table_name="mirbase_conversion",
            task_id="mirbase_conversion",
            project_sn="mirbase_conversion",
            submit_location="mirbase_conversion"
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
