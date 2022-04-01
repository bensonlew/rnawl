# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.08

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class IsPredictAgent(Agent):
    """
    根据输入文件的不同进行分析，进行IS的预测
    """

    def __init__(self, parent):
        super(IsPredictAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "file_dir", "type": "string"},  # 该文件夹下有*.gff,*.faa,*.fna(基因的gff文件，基因蛋白文件，基因的核酸文件)
            {"name": "gene_gff", "type": "string"},  # 预测的基因gff文件
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        if not self.option("input_fa").is_set:
            raise OptionError("必须设置参数input_fa文件!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(IsPredictAgent, self).end()


class IsPredictTool(Tool):
    def __init__(self, config):
        super(IsPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin:" + self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/ISEScan-1.7.2.1:" + self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.python = "/program/Python35/bin/python3"
        self.python_script = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/ISEScan-1.7.2.1/isescan.py"
        self.fasta = self.option("input_fa").prop['path']
        self.fa_file = os.path.basename(self.fasta)
        self.sample = self.fa_file.split(".fna")[0]
        self.python_path = '/program/Python/bin/python'
        self.script_path = self.config.PACKAGE_DIR + '/mobile_genetic_elements/'

    def run_is(self):
        if os.path.exists(self.work_dir + "/" + self.fa_file):
            os.remove(self.work_dir + "/" + self.fa_file)
        os.link(self.fasta, self.work_dir + "/" + self.fa_file)
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        if self.option("file_dir"):
            shutil.copytree(self.option("file_dir"), self.work_dir + "/temp")
        cmd = "{} {} {} {} {} --nthread 4".format(self.python, self.python_script, self.fa_file, "temp", "hmm")
        command = self.add_command("run_is", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_is运行完成！")
        elif command.return_code == 1:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("run_is运行完成运行出错!")

    def get_result(self):
        """
        整理预测出的IS结果
        :return:
        """
        cmd = "{} {}is_chuli.py --g {} --o {} ".format(self.python_path, self.script_path, self.work_dir + "/prediction/" + self.fa_file + ".gff", self.work_dir + '/is_result.xls')
        command = self.add_command("get_result", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_is运行完成！")
        else:
            self.set_error("run_is运行完成运行出错!")

    def chuli_result(self):
        """
        根据整理的结果与gff编码取比较，去除部分错误IS
        :return:
        """
        cmd = "{} {}gff_mobile.py --g {} --m {} --o {} --s 2 --e 3 --i 1".format(self.python_path, self.script_path, self.option("gene_gff"), self.work_dir + '/is_result.xls', self.output_dir + '/' + self.sample + ".is.xls")
        command = self.add_command("chuli_result", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("chuli_result运行完成！")
        else:
            self.set_error("chuli_result运行完成运行出错!")

    def run(self):
        super(IsPredictTool, self).run()
        self.run_is()
        if os.path.exists(self.work_dir + "/prediction/" + self.fa_file + ".gff"):
            self.get_result()
            self.chuli_result()
        self.end()