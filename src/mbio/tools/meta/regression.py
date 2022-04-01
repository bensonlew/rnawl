# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError
import types



class RegressionAgent(Agent):
    """
    version v1.0
    author: gaohao
    last_modified:2017.10.30
    """

    def __init__(self, parent):
        super(RegressionAgent, self).__init__(parent)
        options = [
            {"name": "taxon_table", "type": "infile","format":"meta.otu.otu_table"},
            {"name": "func_table", "type": "infile","format":"meta.otu.otu_table"},
            {"name": "pcnum", "type": "int","default": 3},  # by houshuang 20191009 alpha多样性时只有一个指数
        ]
        self.add_option(options)
        self.step.add_steps('regression')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.regression.start()
        self.step.update()

    def step_end(self):
        self.step.regression.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('taxon_table').is_set:
            raise OptionError('必须提供物种层次的表', code="32706401")
        if not self.option('func_table').is_set:
            raise OptionError('必须提供功能层次的表', code="32706402")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
        ])
        super(RegressionAgent, self).end()


class RegressionTool(Tool):
    def __init__(self, config):
        super(RegressionTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        #self.script1_path = self.config.PACKAGE_DIR + "/statistical/regression_analysis.pl"    
        self.script1_path = self.config.PACKAGE_DIR + "/statistical/regression_analysis_meta.pl"    # guanqing 20180518 替换底层脚本
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run_regression(self):
        cmd = self.perl_path +' %s -i1 %s -i2 %s -c T -o %s -pcnum %d' % (                  # guanqing 20180518 增加pcnum
        self.script1_path, self.option('taxon_table').prop['path'], self.option('func_table').prop['path'], self.work_dir + '/Regression', self.option("pcnum"))
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Regression的r文件生成成功")
        else:
            self.set_error("Regression的r文件生成失败", code="32706403")
            self.set_error("Regression的r文件生成失败", code="32706405")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/RegressionAnalysis.cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算Regression成功")
        else:
            self.set_error("R程序计算Regression失败", code="32706404")
            self.set_error("R程序计算Regression失败", code="32706406")

    def run(self):
        super(RegressionTool, self).run()
        self.run_regression()
        self.set_output()

    def set_output(self):
        linksites_data = os.path.join(self.output_dir, 'Regression.data.xls')
        linksites_r = os.path.join(self.output_dir, 'Regression.message.xls')
        for i in linksites_data,linksites_r:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + '/Regression.data.xls', linksites_data)
        os.link(self.work_dir + '/Regression.message.xls', linksites_r)
        self.end()



