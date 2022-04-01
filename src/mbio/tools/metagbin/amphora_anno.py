# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.01.14

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class AmphoraAnnoAgent(Agent):
    """
    amphora进行bin的物种注释
    """
    def __init__(self, parent):
        super(AmphoraAnnoAgent, self).__init__(parent)
        options = [
            {"name": "amphora_dir", "type": "infile", "format": "metagbin.amphora_dir"},  # 每个bin注释的*.anno.xls文件目录
            {"name":"out","type":"outfile","format":"metagbin.file_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("amphora_dir").is_set:
            raise OptionError("必须设置参数binning的注释文件目录!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AmphoraAnnoAgent, self).end()

class AmphoraAnnoTool(Tool):
    def __init__(self, config):
        super(AmphoraAnnoTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.result = self.work_dir + '/'
        self.summary =self.config.PACKAGE_DIR + "/metagbin/amphora_summary.pl"

    def run_sum_anno(self):
        cmd = '{} {} {} {}'.format(self.perl_path, self.summary,self.option("amphora_dir").prop['path'],self.result)
        command = self.add_command("run_sum_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("汇总bin注释结果文件运行完成!")
            self.set_output()
        else:
            self.set_error("汇总bin注释结果文件运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/summary.anno.xls'):
            os.remove(self.output_dir + '/summary.anno.xls')
        os.link(self.work_dir + '/summary.anno.xls',self.output_dir + '/summary.anno.xls')
        self.option('out', self.output_dir + '/summary.anno.xls')

    def run(self):
        super(AmphoraAnnoTool, self).run()
        self.run_sum_anno()
        self.end()