# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess

class PlotAccumAgent(Agent):
    """
    to perform bubble plot
    author: yiru
    modified at date 2017/3/16
    """
    def __init__(self,parent):
        super(PlotAccumAgent,self).__init__(parent)
        options = [
            {"name": "otu_table","type": "infile","format": "meta.otu.otu_table","default":None}
        ]
        self.add_option(options)
        self.step.add_steps('plot_accum') 
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.plot_accum.start()
        self.step.update()

    def step_end(self):
        self.step.plot_accum.finish()
        self.step.update()

    def check_options(self):
        if not self.option('otu_table').is_set:
            raise OptionError("必须提供otu表")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".","","结果目录"],
            ["./bubble.pdf","pdf","Accum boxplot"]
        ])
        super(PlotAccumAgent,self).end()



class PlotAccumTool(Tool):

    def __init__(self,config):
        super(PlotAccumTool,self).__init__(config)
        self._version = "1.0" 
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        super(PlotAccumTool,self).run()
        self.run_pdf()

    def run_pdf(self):
        cmd = '{}/program/perl/perls/perl-5.24.0/bin/perl {}/bioinfo/meta/scripts/shared_otu_curves.pl'.format(self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR) 
        cmd += ' -i %s ' %(self.option('otu_table').prop['path'])
        cmd += ' -o %s ' %('test')##prefix for o
        print cmd
        if not os.path.exists(self.work_dir + '/accum/'):
            os.mkdir(self.work_dir + '/accum/')
        self.logger.info('运行累积曲线pdf生成程序')
        self.logger.info(cmd) 
        
        try:
            subprocess.check_output(cmd,shell=True)
            self.logger.info('生成 cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 失败')
            self.set_error('无法生成 cmd.r 文件')
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/R --restore --no-save < %s/test.cmd.r' % (self.work_dir + '/accum/'), shell=True) ##add /accum/
            self.logger.info('累积曲线生成成功')
        except subprocess.CalledProcessError:
            self.logger.info('累积曲线生成失败')
            self.set_error('R运行生成累积曲线失败')
            raise "运行R脚本生成累积曲线失败"
        #self.logger.info('运行plot-bubble.pl程序进行npca计算完成')
        print self.output_dir
        print self.work_dir
        if os.path.exists(self.output_dir + '/accum.pdf'):
            os.remove(self.output_dir + '/accum.pdf')
            os.link(self.work_dir + '/accum/accum.pdf', self.output_dir + '/accum.pdf')
        self.logger.info('运行shared_otu_curves.pl程序进行计算完成')
        self.end()


