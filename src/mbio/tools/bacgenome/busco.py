# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.05.28

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class BuscoAgent(Agent):
    """
    真菌基因组bus基因组评估
    """
    def __init__(self, parent):
        super(BuscoAgent, self).__init__(parent)
        options = [
            {"name": "scaf_fa", "type": "infile", "format": "sequence.fasta"},  ##计算时使用base数量
            {"name": "database", "type": "string"},
            {"name": "sample_name", "type": "string"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("scaf_fa").is_set:
            raise OptionError("必须添加scaf_fa的fa文件！", code="32100301")
        if not self.option("database"):
            raise OptionError("必须添加database参数！", code="32100302")

    def set_resource(self):
        self._cpu = 12
        self._memory = '20G'

    def end(self):
        super(BuscoAgent, self).end()

class BuscoTool(Tool):
    def __init__(self, config):
        super(BuscoTool, self).__init__(config)
        self.scaf = self.option("scaf_fa").prop['path']
        self.sample = self.option("sample_name")
        if self.option('database') == "bacteria":
            self.database = self.config.SOFTWARE_DIR + '/database/BUSCO/bacteria'
        else:
            self.database = self.option("database")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/augustus/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/augustus/scripts')
        self.set_environ(AUGUSTUS_CONFIG_PATH=self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/augustus/config')
        self.set_environ(LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + '/program/miniconda2/pkgs/boost-cpp-1.64.0-1/lib')  #add for ningbo zouguanqing 20181029
        self.python_path = "program/Python/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/fungi_genome/"
        self.busco = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/busco-master/scripts/run_BUSCO.py"
        self.file = self.work_dir + '/run_' + self.sample + '/' + 'short_summary_' + self.sample + '.txt'

    def run_busco(self):
        self.logger.info('正在进行busco完整性预测')
        cmd = '{} {} --in {} -m genome -l {} -o {} -f'.format(self.python_path,self.busco,self.scaf,self.database,self.sample)
        command = self.add_command("run_busco", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_busco运行完成")
        else:
            self.set_error("run_busco运行出错!", code="32100301")

    def run_stat(self):
        cmd = '{} {}tiqu_busco.pl {} {}'.format(self.perl_path,self.perl_script,self.file,self.work_dir + '/' + self.sample + '_busco.xls')
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_stat运行完成")
        else:
            self.set_error("run_stat运行出错!", code="32100302")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if os.path.exists(self.output_dir + '/' + self.sample + '_busco.xls'):
            os.remove(self.output_dir + '/' + self.sample + '_busco.xls')
        os.link(self.work_dir + '/' + self.sample + '_busco.xls', self.output_dir + '/' + self.sample + '_busco.xls')

    def run(self):
        super(BuscoTool, self).run()
        self.run_busco()
        self.run_stat()
        self.set_output()
        self.end()