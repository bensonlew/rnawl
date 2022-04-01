# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir,link_file


class IsPredictAgent(Agent):
    """
    根据输入文件的不同进行分析，进行IS的预测
    """

    def __init__(self, parent):
        super(IsPredictAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "file_dir", "type": "string"},  # 该文件夹下有*.gff,*.faa,*.fna(基因的gff文件，基因蛋白文件，基因的核酸文件)
            {"name": "sample", "type": "string"}, ## 传入的样本名称
        ]
        self.add_option(options)

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
        self.sample = self.option("sample")
        self.python_path = '/program/Python/bin/python'
        self.script_path = self.config.PACKAGE_DIR + '/bacgenome/'

    def run_is(self):
        """
        运行ISEScan软件进行计算，注意这里修改了软件计算方式，要具体看软件内容
        """
        if os.path.exists(self.work_dir + "/" + self.sample):
            os.remove(self.work_dir + "/" + self.sample)
        os.link(self.fasta, self.work_dir + "/" + self.sample)
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        if self.option("file_dir"):
            shutil.copytree(self.option("file_dir"), self.work_dir + "/temp")
        cmd = "{} {} {} {} {} --nthread 4".format(self.python, self.python_script, self.sample, "temp", "hmm")
        command = self.add_command("run_is", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_is运行完成！")
        else:
            self.set_error("run_is运行完成运行出错!")

    def get_result(self):
        """
        整理预测出的IS结果
        :return:
        """
        cmd = "{} {}is_chuli.py --g {} --o {} ".format(self.python_path, self.script_path, self.work_dir + "/prediction/" + self.sample + ".gff", self.work_dir + '/is_result.xls')
        self.logger.info(cmd)
        command = self.add_command("get_result", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("format_gff运行完成！")
        else:
            self.set_error("format_gff运行完成运行出错!")

    def pick_sequence(self):
        """
        提取序列
        :return:
        """
        stat_file = self.work_dir + '/' + self.sample + ".is.xls"
        # gene_ffn = os.path.join(self.option("file_dir"), self.sample + ".ffn")
        cmd = "{} {}tiqu_sequence2.py {} {} {}".format(self.python_path, self.script_path, stat_file, self.fasta, self.sample)
        self.logger.info(cmd)
        command = self.add_command("pick_sequence", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pick_sequence运行完成！")
        else:
            self.set_error("pick_sequence运行出错!")

    def chuli_result(self):
        """
        根据整理的结果与gff编码取比较，去除部分错误IS
        :return:
        """
        cmd = "{} {}gff_mobile.py --g {} --m {} --o {} --s 2 --e 3 --i 1".format(self.python_path, self.script_path, os.path.join(self.option("file_dir"), self.sample+".gff"), self.work_dir + '/is_result.xls', self.work_dir + '/' + self.sample + ".is.xls")
        command = self.add_command("chuli_result", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("chuli_result运行完成！")
        else:
            self.set_error("chuli_result运行完成运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录!")
        link_dir(self.work_dir + "/prediction/", self.output_dir)
        # link_file(self.work_dir +"/"+ self.sample + ".gene.fna", self.work_dir +"/"+ self.sample + ".gene.fna")
        # link_file(self.work_dir +"/"+ self.sample + "_sequence.fna", self.work_dir +"/"+ self.sample + "_sequence.fna")
        if os.path.exists(self.work_dir + '/' + self.sample + ".is.xls"):
            link_file(self.work_dir + '/' + self.sample + ".is.xls", self.output_dir + '/' + self.sample + ".is.xls")
        if os.path.exists(self.work_dir +"/"+ self.sample + ".gene.fna"):
            link_file(self.work_dir +"/"+ self.sample + ".gene.fna", self.output_dir +"/"+ self.sample + ".gene.fna")
        if os.path.exists(self.work_dir +"/"+ self.sample + "_sequence.fna"):
            link_file(self.work_dir +"/"+ self.sample + "_sequence.fna", self.output_dir +"/"+ self.sample + "_sequence.fna")
        if os.path.exists(self.work_dir +"/"+ self.sample + ".stat.xls"):
            link_file(self.work_dir +"/"+ self.sample + ".stat.xls", self.output_dir +"/"+ self.sample + ".stat.xls")

    def run(self):
        super(IsPredictTool, self).run()
        self.run_is()
        if os.path.exists(self.work_dir + "/prediction/" + self.sample + ".gff"):
            self.get_result()
            self.chuli_result()
            if os.path.exists(self.work_dir + '/' + self.sample + ".is.xls"):
                self.pick_sequence()
            self.set_output()
        self.end()