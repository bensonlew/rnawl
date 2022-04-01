# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2021.01.18

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import os,shutil


class RrnaTreeAgent(Agent):
    """
    16s序列进行比对对齐和掩盖低置信度的位置处理
    """

    def __init__(self, parent):
        super(RrnaTreeAgent, self).__init__(parent)
        options = [
            {"name": "s16_fa", "type": "infile", "format": "sequence.fasta"},  # 所有16s的序列
            {"name": "sample", "type": "string"},  ##样品名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("s16_fa").is_set:
            raise OptionError("必须输入16s_fa文件夹！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(RrnaTreeAgent, self).end()


class RrnaTreeTool(Tool):
    def __init__(self, config):
        super(RrnaTreeTool, self).__init__(config)
        self._version = "1.0"
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/ssu-align-0.1.1/bin:"
        self.manpath = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/ssu-align-0.1.1/share/man"
        self.ssu = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/ssu-align-0.1.1/share/ssu-align-0.1.1"
        self.set_environ(PATH=self.path, MANPATH=self.manpath, SSUALIGNDIR=self.ssu)
        self.sof_path ="/bioinfo/compare_genome/software/ssu-align-0.1.1/bin/"
        self.s16_fa = self.option("s16_fa").prop["path"]
        self.iqtree = "/bioinfo/compare_genome/software/iqtree-1.6.12-Linux/bin/iqtree"

    def run(self):
        """
        运行
        :return:
        """
        super(RrnaTreeTool, self).run()
        self.run_align()
        self.run_mask()
        self.run_tree()
        self.set_output()
        self.end()

    def run_align(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        cmd = '{}ssu-align -n bacteria {} {}'.format(self.sof_path, self.s16_fa, "temp")
        command = self.add_command("run_align", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_align运行完成！")
        else:
            self.set_error("run_align运行完成运行出错!")

    def run_mask(self):
        cmd = "{}ssu-mask --afa {}".format(self.sof_path, "temp")
        command = self.add_command("run_mask", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_mask运行完成！")
        else:
            self.set_error("run_mask运行完成运行出错!")

    def run_tree(self):
        if os.path.exists(self.work_dir + "/temp2"):
            shutil.rmtree(self.work_dir + "/temp2")
        os.mkdir(self.work_dir + "/temp2")
        cmd = '{} -s {} -pre {} -m MFP -nt AUTO -ntmax 8 -bb {}'.format(self.iqtree, self.work_dir + "/temp/temp.bacteria.mask.afa", self.work_dir + "/temp2/all", 1000)
        command = self.add_command("run_tree", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_tree运行完成！")
        else:
            self.set_error("run_tree运行完成运行出错!")

    def set_output(self):
        link_file(self.work_dir + "/temp/temp.bacteria.mask.afa", self.output_dir + "/all.16s.align.fasta")
        link_file(self.work_dir + "/temp2/all.treefile", self.output_dir + "/"+self.option("sample") +".16s.nwk")