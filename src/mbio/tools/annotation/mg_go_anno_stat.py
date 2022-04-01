# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class MgGoAnnoStatAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: gaohao
    last_modified: 20180929
    """

    def __init__(self, parent):
        super(MgGoAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "go_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "level", "type": "string","default":"all"},  # 各个分析计算单个水平丰度用
            {"name": "level_out", "type": "outfile", "format": "meta.otu.otu_table"} # 单水平输出丰度文件
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40 # 每次重运行增加内存40G by qingchen.zhang @ 20190527

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if not self.option('go_anno').is_set:
            raise OptionError("Must provide go_anno result file !", code="31204101")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("Must provide reads_profile_table result file !", code="31204102")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(MgGoAnnoStatAgent, self).end()


class MgGoAnnoStatTool(Tool):

    def __init__(self, config):
        super(MgGoAnnoStatTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.script3 = self.config.PACKAGE_DIR + '/annotation/all_go_func_anno.pl'
        self.result_name = 'GO'
        self.anno = self.option('go_anno').prop['path']
        self.abund = self.option('reads_profile_table').prop['path']
        self.out = self.work_dir + '/sample'


    def run(self):
        super(MgGoAnnoStatTool, self).run()
        self.run_go_anno_all()
        self.set_output()
        self.end()

    def run_go_anno_all(self):
        cmd = '{} {} {} {}'.format(self.perl_path, self.script3, self.anno,self.abund)
        if self.option("level"):
            cmd += " '{}' ".format(self.option("level"))
        command = self.add_command('run_go_anno_all', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            if self.option("level") == "all":
                self.get_level4(self.work_dir + "/all.go1234.function.xls",self.work_dir + "/go.level4.xls")
            self.logger.info("run_go_anno_all succeed")
        elif command.return_code in [-9, 1]:  # change return_code by qingchen.zhang @ 20190604
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20190604
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("run_go_anno_all failed", code="31204101")

    def set_output(self):
        if self.option("level") == "all":
            for file in ["all.go1.function.xls","all.go12.function.xls","all.go123.function.xls","all.go1234.function.xls","go.level4.xls"]:
                if os.path.exists(self.output_dir + '/' + file):
                    os.remove(self.output_dir + '/' + file)
                os.link(self.work_dir + '/' + file,self.output_dir + '/' + file)
        else:
            self.option("level_out", self.work_dir + "/level.xls")

    def get_level4(self,file,output):
        with open (file,'r') as f,open (output,'w') as p:
            files =f.readlines()
            for line in files:
                temp =line.split("\t")
                des="\t".join(temp[7:])
                p.write(temp[5] + "\t" + des)