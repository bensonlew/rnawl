# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re

class PhiAnnoStatAgent(Agent):
    """
    宏基因phi注释结果丰度统计表
    author: guhaidong
    last_modify:
    """

    def __init__(self, parent):
        super(PhiAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "phi_anno_table", "type": "infile", "format": "sequence.profile_table"},
            # 基因注释具体信息结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("phi_anno_table").is_set:
            raise OptionError("必须设置注释文件", code="31204601")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置基因丰度文件", code="31204602")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'  # 改回5G by GHD @20180428

    def end(self):
        super(PhiAnnoStatAgent, self).end()

class PhiAnnoStatTool(Tool):
    def __init__(self, config):
        super(PhiAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        self.script = self.config.PACKAGE_DIR + '/annotation/scripts/phi_anno_abudance.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(PhiAnnoStatTool, self).run()
        self.run_phi_stat()
        self.set_output()
        self.end()

    def run_phi_stat(self):
        self.logger.info("start phi_stat")
        cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.script, self.option('phi_anno_table').prop['path'],self.option('reads_profile_table').prop['path'], self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('phi_profile', cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("phi_stat succeed")
        elif command.return_code in [-9, 1]:  # change return_code by ghd @ 20190110
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.set_error("phi_stat failed", code="31204601")
            self.set_error("phi_stat failed", code="31204602")

    def set_output(self):
        self.logger.info("start set_output")
