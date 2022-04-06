# -*- coding: utf-8 -*-
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.itraq_and_tmt.go_graph import draw_GO
import glob
import unittest

class GoGeneGmxAgent(Agent):
    def __init__(self, parent):
        super(GoGeneGmxAgent, self).__init__(parent)
        options = [
            {"name": "go_list", "type": "string", "default": None},
            {"name": "min_num", "type": "string", "default": "15"},
            {"name": "max_num", "type": "string", "default": "500"},
            {"name": "go_gmx", "type": "string", "default": "go.gmt"},
            {"name": "go_type", "type": "string", "default": "all"},
            {"name": "go_sets", "type": "string", "default": "all"},

        ]
        self.add_option(options)
        self.step.add_steps("go_dag")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.go_dag.start()
        self.step.update()

    def stepfinish(self):
        self.step.go_dag.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option("go_list")):
            raise OptionError("%s not exist", variables = (self.option("go_list")), code="33711102")

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["go.gmx", "gmx", "go对应基因列表"]
        ])
        super(GoGeneGmxAgent,self).end()

class GoGeneGmxTool(Tool):
    def __init__(self,config):
        super(GoGeneGmxTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.go2gene = self.config.PACKAGE_DIR + "/medical_transcriptome/go2gene.py"
        self.obo = self.config.SOFTWARE_DIR + '/database/GO/go-basic.obo'

    def run_dag(self):
        if self.option("go_sets") == "" or self.option("go_sets") == "all":
            go_sets = "all"
        else:
            go_sets = self.option("go_sets")
        cmd = self.python_path + " {} {} {} {} {} {} {}".format(self.go2gene,
                                                            self.option("go_list"),
        self.option("go_type"),
        self.option("min_num"),
        self.option("max_num"),
        go_sets,
        self.option("go_gmx"))
        cmd_obj = self.add_command("cmd", cmd, ignore_error=True).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("cmd list执行成功")
        self.logger.info("go2gene生成完毕")


    def set_output(self):
        go_dag = glob.glob(self.work_dir + '/go.gmt')
        for each in go_dag:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(GoGeneGmxTool, self).run()
        self.run_dag()
        self.set_output()
        self.end()
