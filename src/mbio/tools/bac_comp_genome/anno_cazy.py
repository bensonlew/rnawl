# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.align.blast.blastout_statistics import *
import pandas as pd


class AnnoCazyAgent(Agent):
    """
    利用脚本对比对结果进行统计
    """
    def __init__(self, parent):
        super(AnnoCazyAgent, self).__init__(parent)
        options = [
            {"name": "hmmscan_result", "type": "infile", "format": "meta_genomic.hmmscan_table"},  # 比对结果文件
            {"name": "best", "type": "bool", "default": True},  # 是否一条序列只挑选一个结果
            {"name": "add_score", "type": "string", "default": 'False'},  # 是否添加score和identity
        ]
        self.add_option(options)
        self.step.add_steps("cazy_anno")
        self.on("start", self.step_start)
        self.on("end", self.step_end)
        self._cpu = 1
        self._memory = ''
        self._memory_increase_step = 30

    def step_start(self):
        self.step.cazy_anno.start()
        self.step.update()

    def step_end(self):
        self.step.cazy_anno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("hmmscan_result").is_set:
            raise OptionError("必须提供hmmscan比对结果作为输入文件", code="31200801")
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = "15G"

    def end(self):
        super(AnnoCazyAgent, self).end()


class AnnoCazyTool(Tool):
    def __init__(self, config):
        super(AnnoCazyTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/cazy_anno_comp.py"
        self.class_def = self.config.SOFTWARE_DIR + "/database/CAZyDB/class_definition.txt"
        self.FamInfo = self.config.SOFTWARE_DIR + "/database/CAZyDB/FamInfo.txt"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin')
        self.set_environ(PERLBREW_ROOT=self.config.SOFTWARE_DIR + '/program/perl')

    def run(self):
        super(AnnoCazyTool, self).run()
        self.run_annot()
        self.end()

    def run_annot(self):
        cmd1 = '{} {} --out.dm {} --output_dir {} --class_def {} --FamInfo {}'.\
            format(self.python_path, self.script_path, self.option('hmmscan_result').prop['path'],
                   self.output_dir + "/gene_", self.class_def, self.FamInfo)
        if self.option("best"):
            cmd1 += ' -best True '
        if "add_score" in self.get_option_object().keys():
            cmd1 += ' -add_score {}'.format(self.option("add_score"))
        self.logger.info("start cazy_anno")
        command = self.add_command("cazy_anno", cmd1, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cazy_anno done")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("cazy_anno error", code="31200801")
            raise Exception("cazy_anno error")

        if "add_score" in self.get_option_object().keys():
            if self.option("add_score") == 'True':
                if os.path.getsize(self.output_dir+"/gene_dbCAN.hmmscan.out.dm.ds") >0:
                    data = pd.read_table(self.output_dir + "/gene_dbCAN.hmmscan.out.dm.ds")
                    data['Coverd_fraction'] = (data['Query_end'] - data['Query_start']) / data['Query_length'] * 100
                    data.to_csv(self.output_dir + "/gene_dbCAN.hmmscan.out.dm.ds", sep='\t', index=False)
