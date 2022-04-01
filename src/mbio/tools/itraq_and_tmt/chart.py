# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

import os
import unittest
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import json
from mbio.packages.itraq_and_tmt.chart import Chart

from mbio.packages.ref_rna_v2.xml2table2 import xml2table
from mbio.packages.rna.annot_config import AnnotConfig
# from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
# from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
from mbio.packages.itraq_and_tmt.chart_report import ChartReport


class ChartAgent(Agent):
    """
    Blast比对到Rfam库
    """

    def __init__(self, parent):
        super(ChartAgent, self).__init__(parent)
        options = [
            {"name": "file_json", "type": "infile", "format": "ref_rna_v2.common"},  # 输入序列文件
            {"name": "chart_list", "type": "string", "default": None},  #
            {'default': 10, 'type': 'int', 'name': 'cpu'},
            {"name": "chart_report", "type": "bool", "default": True}

        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("file_json").is_set:
            raise OptionError("必须提供file_json输入文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '180G'

    def end(self):
        super(ChartAgent, self).end()


class ChartTool(Tool):
    def __init__(self, config):
        super(ChartTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.parafly = software_dir + '/bioinfo/denovo_rna_v2/trinityrnaseq-2.8.5/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self._LD_LIBRARY_PATH = "$LD_LIBRARY_PATH:" + software_dir + "/bioinfo/sg_chart/miniconda2/lib"
        self._NODE_PATH = software_dir + "/bioinfo/sg_chart/node-v14.16.0-linux-x64/lib/node_modules"
        self.set_environ(LD_LIBRARY_PATH=self._LD_LIBRARY_PATH, NODE_PATH=self._NODE_PATH)


    def run(self):
        """
        运行
        :return:
        """ 
        super(ChartTool, self).run()
        self.chart(self.option("file_json").prop['path'])
        self.parafly_pdf()
        if self.option('chart_report'):
            self.chart_report()
        self.end()

    def chart(self, json_file):
        """
        ;json 配置文件
        """
        with open(json_file, 'r') as f:
            a = json.loads(f.read())
            if a["type"] == "workflow":
                chart_obj = Chart()
            elif a["type"] == "geneset":
                chart_obj = Chart()
            elif a["type"] == "advance":
                chart_obj = Chart()
            chart_obj.chart_json_batch(a)

    def parafly_pdf(self):
        cmd_list = open('para_run.sh','r').read().strip().split('\n')
        # cmd_str = ""
        i = 1
        for cmd in cmd_list:
            if "puppeteer" in cmd:
                cmd_pup = cmd.strip()
                cmd_name = 'cmd_pup' + str(i)
                try:
                    return_code = os.system(cmd_pup)
                    if return_code == 0:
                        self.logger.info("{} Finished successfully".format(cmd_name))
                    else:
                        self.logger.info("{} unfinished command".format(cmd_name))
                except:
                    try:
                        return_code = os.system(cmd_pup)
                        if return_code == 0:
                            self.logger.info("{} Finished successfully".format(cmd_name))
                        else:
                            self.logger.info("{} unfinished command".format(cmd_name))
                    except:
                        self.set_error("{} failed".format(cmd_name))
                # command = self.add_command(cmd_name, cmd_pup, ignore_error=True, shell=True)
                # command.run()
                # self.wait()
                # if command.return_code == 0:
                #     self.logger.info("{} Finished successfully".format(cmd_name))
                # else:
                #     self.logger.info("{} unfinished command".format(cmd_name))
                i += 1
                # cmd_str += " && sleep 8 && " + cmd.strip()
        os.system("sed -i '/puppeteer/d' para_run.sh")
        cmd = '{} '.format(self.parafly)
        cmd += '-{} {} '.format("c", "para_run.sh")
        cmd += '-{} {} '.format("CPU", self.option("cpu"))
        cmd += '-v -shuffle'
        # cmd += cmd_str
        cmd_name = 'phatomjstopdf'
        command = self.add_command(cmd_name, cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.logger.info("{} unfinished command".format(cmd_name))

    def chart_report(self):
        report = ChartReport()
        report.run()




class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "chart" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.chart",
            "instant": False,
            "options": dict(
                file_json="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2/chart_workflow.json"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
