#!/usr/bin/python
# -*- coding:utf-8 -*-
# __author__ = "konghualei, 20170426"
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.ref_rna.express.single_sample import *
from mbio.packages.ref_rna.express.set_strand import set_strand
from mbio.packages.denovo_rna.express.express_distribution import *
from mbio.packages.ref_rna.express.genesetVenn import *
import shutil
import os
import re

class GenesetVennAgent(Agent):
    def __init__(self,parent):
        super(GenesetVennAgent, self).__init__(parent)
        options = [
            {"name":"fpkm","type":"infile","format":"rna.express_matrix"},
            {"name":"graph","type":"outfile","format":"rna.express_matrix"},
            {"name":"table","type":"outfile","format":"rna.express_matrix"}
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps("genesetvenn")
        self.on("start",self.stepstart)
        self.on("end",self.stepfinish)

    def stepstart(self):
        self.step.genesetvenn.start()
        self.step.update()

    def stepfinish(self):
        self.step.genesetvenn.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fpkm").is_set:
            raise Exception("请输入fpkm文件！")

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["./venn_graph.xls", "xls", "venn图输出目录1"],
            ["./venn_table.xls", "xls", "venn图输出目录2"]
        ])
        super(GenesetVennAgent, self).end()

class GenesetVennTool(Tool):
    def __init__(self, config):
        super(GenesetVennTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR+"/miniconda2/bin:$PATH"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path1 = "/program/R-3.3.1/bin/Rscript "


    def genesetvenn_run(self):
        rfile=self.work_dir+"/genesetvenn.r"
        ExpressVenn(rfile,self.option("fpkm").prop['path'],self.output_dir)
        cmd1=self.r_path1+rfile
        cmd = self.add_command("genesetvenn", cmd1, ignore_error=True).run()
        self.wait(cmd)
        if cmd.return_code==0:
            self.logger.info("%s运行完成" % cmd1)
        elif cmd.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("%s运行出错" % cmd1)

    def set_output(self):
        tmp = self.output_dir + "/tmp"
        with open(self.option("fpkm").prop['path'],'r+') as f1,open(tmp,'w+') as f2:
            f2.write("#group_name\tspecies_name\n")
            for lines in f1:
                f2.write(lines)
        shutil.copy2(tmp,self.output_dir+"/venn_graph.xls")
        os.remove(tmp)
        for files in os.listdir(self.output_dir):
            print files
            if re.search(r'table',files):
                self.option("table").set_path(self.output_dir+"/"+files)

    def run(self):
        super(GenesetVennTool,self).run()
        self.genesetvenn_run()
        self.set_output()
        self.end()
