# -*- coding: utf-8 -*-
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.itraq_and_tmt.go_graph import draw_GO
import glob
import unittest

class DoGeneGmxAgent(Agent):
    def __init__(self, parent):
        super(DoGeneGmxAgent, self).__init__(parent)
        options = [
            {"name": "do_list", "type": "string", "default": None},
            {"name": "min_num", "type": "string", "default": "15"},
            {"name": "max_num", "type": "string", "default": "500"},
            {"name": "do_gmx", "type": "string", "default": "do.gmt"},
            {"name": "do_sets", "type": "string", "default": "all"},

        ]
        self.add_option(options)
        self.step.add_steps("do_dag")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.do_dag.start()
        self.step.update()

    def stepfinish(self):
        self.step.do_dag.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["do.gmx", "gmx", "do对应基因列表"]
        ])
        super(DoGeneGmxAgent,self).end()

class DoGeneGmxTool(Tool):
    def __init__(self,config):
        super(DoGeneGmxTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.do2gene = self.config.PACKAGE_DIR + "/medical_transcriptome/do2gene.py"

    def run_dag(self):
        if self.option("do_sets") == "" or self.option("do_sets") == "all":
            do_sets = "all"
        else:
            do_sets = self.option("do_sets")
        cmd = self.python_path + " {} {} {} {} {} {}".format(self.do2gene,
                                                            self.option("do_list"),
        self.option("min_num"),
        self.option("max_num"),
        do_sets,
        self.option("do_gmx"))
        cmd_obj = self.add_command("cmd", cmd, ignore_error=True).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("cmd list执行成功")
        self.logger.info("do2gene生成完毕")


    def set_output(self):
        do_dag = glob.glob(self.work_dir + '/do.gmt')
        for each in do_dag:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(DoGeneGmxTool, self).run()
        self.run_dag()
        self.set_output()
        self.end()
