# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.06.04

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class HeterStatAgent(Agent):
    """
    真菌基因组杂合率评估统计
    """
    def __init__(self, parent):
        super(HeterStatAgent, self).__init__(parent)
        options = [
            {"name": "frequency", "type": "infile", "format": "sequence.profile_table"},  ##计算时kmer=17时的频率
            {"name": "heter1", "type": "infile", "format": "sequence.profile_table"},
            {"name": "heter2", "type": "infile", "format": "sequence.profile_table"},
            {"name": "heter3", "type": "infile", "format": "sequence.profile_table"},
            {"name": "heter4", "type": "infile", "format": "sequence.profile_table"},
            {"name": "heter5", "type": "infile", "format": "sequence.profile_table"},
            {"name": "heter6", "type": "infile", "format": "sequence.profile_table"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("frequency").is_set:
            raise OptionError("必须添加frequency文件！", code="32101301")
        if not self.option("heter1").is_set:
            raise OptionError("必须添加heter1的文件！", code="32101302")
        if not self.option("heter2").is_set:
            raise OptionError("必须添加heter2的文件！", code="32101303")
        if not self.option("heter3").is_set:
            raise OptionError("必须添加heter3的文件！", code="32101304")
        if not self.option("heter4").is_set:
            raise OptionError("必须添加heter4的文件！", code="32101305")
        if not self.option("heter5").is_set:
            raise OptionError("必须添加heter5的文件！", code="32101306")
        if not self.option("heter6").is_set:
            raise OptionError("必须添加heter6的文件！", code="32101307")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(HeterStatAgent, self).end()

class HeterStatTool(Tool):
    def __init__(self, config):
        super(HeterStatTool, self).__init__(config)
        self.fre = self.option("frequency").prop['path']
        self.heter1 = self.option("heter1").prop['path']
        self.heter2 = self.option("heter2").prop['path']
        self.heter3 = self.option("heter3").prop['path']
        self.heter4 = self.option("heter4").prop['path']
        self.heter5 = self.option("heter5").prop['path']
        self.heter6 = self.option("heter6").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/fungi_genome/"
        self.new_file = self.work_dir + '/' + 'frequency.new.xls'
        self.result = self.work_dir + '/' + 'all_heter_kmer_freq.xls'

    def run_frequency(self):
        cmd = '{} {}frequency_combine.pl {} {}'.format(self.perl_path,self.perl_script,self.fre,self.new_file)
        command = self.add_command("run_frequency", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_frequency运行完成")
        else:
            self.set_error("run_frequency运行出错!", code="32101301")

    def run_combine(self):
        cmd = '{} {}file_combine.pl {} {} {} {} {} {} {} -l "0:0;1:3;2:3;3:3;4:3;5:3;6:3" -s depth,Target,H0.001,H0.002,H0.003,H0.005,H0.01,H0.02 -o {}'.format(self.perl_path,self.perl_script,self.new_file,self.heter1,self.heter2,self.heter3,self.heter4,self.heter5,self.heter6,self.result)
        command = self.add_command("run_combine", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_combine运行完成")
        else:
            self.set_error("run_combine运行出错!", code="32101302")

    def set_output(self):
        if os.path.exists(self.output_dir + '/all_heter_kmer_freq.xls'):
            os.remove(self.output_dir + '/all_heter_kmer_freq.xls')
        os.link(self.result,self.output_dir + '/all_heter_kmer_freq.xls')

    def run(self):
        super(HeterStatTool, self).run()
        self.run_frequency()
        self.run_combine()
        self.set_output()
        self.end()