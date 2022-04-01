# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210908'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
import pandas as pd


class MirnasearchWorkflow(Workflow):
    """
    Mirnasearch base python package Mirnasearch
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MirnasearchWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "direction", "type": "string", 'default': 'gene_to_mirna'}, # gene_to_mirna/mirna_to_gene
            {"name": "speice", "type": "string", 'default': 'Human'}, # Human/Mouse/Rat
            {"name": "database", "type": "string", 'default': 'tarbase'}, # tarbase/mirdb/mirtarbase/all
            {"name": "search_string", "type": "string", 'default': ""},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.mirnasearch_tools = []
        self.set_options(self._sheet.options())

    def run(self):
        self.run_Mirnasearch()
        super(MirnasearchWorkflow, self).run()

    def check_options(self):
        if not self.option("search_string"):
            raise OptionError("必须设置搜索数据")
        return True

    def run_Mirnasearch(self):
        if self.option("direction") == "gene_to_mirna":
            search_item = "gene"
        if self.option("direction") == "mirna_to_gene":
            search_item = "mirna"
        for search_thing in str(self.option("search_string")).split(","):
            one_mirnasearch_tool = self.add_tool("tool_lab.mirnasearch")
            opts = {
                "speice":self.option('speice'),
                search_item:search_thing,
                "database":self.option('database')
            }
            one_mirnasearch_tool.set_options(opts)
            self.mirnasearch_tools.append(one_mirnasearch_tool)
        self.on_rely(self.mirnasearch_tools, self.set_db)
        self.mirnasearch_tools_all_output_files = []
        for tool in self.mirnasearch_tools:
            self.mirnasearch_tools_all_output_files.append(os.path.join(tool.output_dir, "result.txt"))
            tool.run()
        # self.Mirnasearch.on('end', self.set_db)
        # self.Mirnasearch.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        mirnasearch_tools_all_output_files_df = pd.concat([pd.read_table(i,sep="\t") for i in self.mirnasearch_tools_all_output_files if os.path.exists(i)])
        mirnasearch_tools_all_output_files_df.to_csv(os.path.join(self.output_dir,"all_results.txt"), sep='\t',index=False)
        mongo_api = self.api.api('tool_lab.mirnasearch')
        mongo_api.add_mirnasearch(
            result_file=os.path.join(self.output_dir,"all_results.txt"),
            main_id = self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "预测结果",0],
            ['./all_results.txt', '', '预测结果', 0],
        ])
        super(MirnasearchWorkflow, self).end()


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