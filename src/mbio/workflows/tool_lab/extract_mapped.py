# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import glob
import unittest
import pandas as pd
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class ExtractMappedWorkflow(Workflow):
    """
    This script is used to extracts mapped/unmapped reads from input bam/sam file.
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExtractMappedWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "extract_type", "type": "string", "default": 'mapped'},
            {"name": "fq_type", "type": "string", "default": 'PE'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.extract_mapped")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(ExtractMappedWorkflow, self).run()

    def check_options(self):
        if not self.option('input_file').is_set:
            raise OptionError('SAM/BAM文件必须输入')
        if self.option('fq_type').upper() not in ["PE", "SE"]:
            raise OptionError('测序类型参数输入有误')
        if self.option('extract_type').lower() not in ["mapped", "unmapped"]:
            raise OptionError('暂不支持提取该类型: {}'.format(self.option('extract_type')))
        return True

    def run_tool(self):
        opts = {
            'input_file': self.option('input_file'),
            'extract_type': self.option('extract_type').lower(),
            'fq_type': self.option('fq_type').upper(),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Mapped/Unmapped reads 提取结果文件",0],
        ])
        super(ExtractMappedWorkflow, self).end()


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
            input_file="/mnt/ilustre/users/sanger-dev/workspace/20200109/Refrna_tsg_36796/RnaseqMapping/output/bam/A1_1.bam",
            extract_type="mapped",
            fq_type="PE"
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.extract_mapped",
            main_table_name="extract_mapped",
            task_id="extract_mapped",
            project_sn="extract_mapped",
            submit_location="extract_mapped"
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
