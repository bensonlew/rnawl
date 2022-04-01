# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import glob
import unittest
import pandas as pd
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class PadjustWorkflow(Workflow):
    """
    Given a set of p-values, returns p-values adjusted using one of several methods.
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PadjustWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pvalue", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "method", "type": "string", "default": 'BH'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = list()
        self.methods = self.option("method").split(",")

    def run(self):
        self.run_tool()
        super(PadjustWorkflow, self).run()

    def check_options(self):
        if not self.option('pvalue').is_set:
            raise OptionError('pvalue文件必须输入')
        for method in self.option("method").split(","):
            if method not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", 'BH/fdr']:
                raise OptionError('暂不支持该校正方法')
        return True

    def run_tool(self):
        for method in self.methods:
            if method == "BH/fdr":
                method = "BH"
            tool = self.add_tool('tool_lab.padjust')
            opts = {
                'pvalue': self.option('pvalue'),
                'method': method,
            }
            tool.set_options(opts)
            self.tools.append(tool)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            self.tools[0].run()
        else:
            self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()

    def set_output(self):
        results = list()
        results_pd = list()
        for each_tool in self.tools:
            results.append(glob.glob(each_tool.output_dir + "/*.txt")[0])
        for result in results:
            result_pd = pd.read_table(result, index_col=0, header=0)
            results_pd.append(result_pd)
        result_table = pd.concat(results_pd, axis=1)
        result_table.to_csv(self.output_dir + '/padjust.txt', sep='\t')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "pvalue校正结果文件",0],
            [r'padjust.txt', 'txt', 'pvalue校正结果文件', 0],
        ])
        super(PadjustWorkflow, self).end()


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
            pvalue="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/pvalue",
            method="BH/fdr,BY,hommel,holm",
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.padjust",
            main_table_name="padjust",
            task_id="padjust",
            project_sn="padjust",
            submit_location="padjust"
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
