# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210901'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class MulticorrWorkflow(Workflow):
    """
    Multicorr base R package Multicorr
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MulticorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile_a", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "infile_b", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "infile_group", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "r_threshold", "type": "string", 'default': "0.8"},
            {"name": "p_threshold", "type": "string", 'default': "0.05"},
            {"name": "corr_method", "type": "string", 'default': "pearson"}, # "pearson","spearman","kendall"
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.Multicorr = self.add_tool("tool_lab.multicorr")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_Multicorr()
        super(MulticorrWorkflow, self).run()

    def check_options(self):
        if not self.option("infile_a").is_set:
            raise OptionError("必须设置输入表达量文件")
        if not self.option("infile_b").is_set:
            raise OptionError("必须设置输入表达量文件")
        return True

    def run_Multicorr(self):
        opts = {
            "infile_a":self.option('infile_a').prop['path'],
            "infile_b":self.option('infile_b').prop['path'],
            "infile_group":self.option('infile_group').prop['path'],
            "r_threshold":self.option('r_threshold'),
            "p_threshold":self.option('p_threshold'),
            "corr_method":self.option('corr_method'),
        }
        self.Multicorr.set_options(opts)
        self.Multicorr.on('end', self.set_db)
        self.Multicorr.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        mongo_api = self.api.api('tool_lab.multicorr')
        mongo_api.add_multicorr(
            outfile_dir=self.Multicorr.output_dir,
            s3_dir = self._sheet.output,
            main_id = self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.Multicorr.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "组间相关性分析",0],
            ['./cor*', '', '显著相关组', 0],
            ['./result_cor.txt', '', '相关性总览', 0],
        ])
        super(MulticorrWorkflow, self).end()


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