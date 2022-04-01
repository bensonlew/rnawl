# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2021.01.19

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import re,os
import shutil

class HouseTreeAgent(Agent):
    """
    看家基因比对、去低值区、构建进化树
    """

    def __init__(self, parent):
        super(HouseTreeAgent, self).__init__(parent)
        options = [
            {"name": "seq", "type": "infile", "format": "sequence.fasta"}, #输入文件
            {"name": "sample", "type": "string"},  ##样品名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq").is_set:
            raise OptionError("必须输入seq文件！")
        return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(HouseTreeAgent, self).end()


class HouseTreeTool(Tool):
    def __init__(self, config):
        super(HouseTreeTool, self).__init__(config)
        self._version = "1.0"
        self.seq = self.option("seq").prop["path"]
        self.sh = "../../../../.." + self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.mafft = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/mafft-7.429-without-extensions/binaries/mafft"
        self.trimal = "/bioinfo/compare_genome/software/trimal-trimAl/source/trimal"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.num = self.get_num(self.seq)
        self.align = self.work_dir + "/all.align.fasta"
        self.trimal_fa = self.work_dir + "/all.align_last.fasta"
        self.iqtree = "/bioinfo/compare_genome/software/iqtree-1.6.12-Linux/bin/iqtree"

    def run(self):
        """
        运行
        :return:
        """
        super(HouseTreeTool, self).run()
        self.run_align()
        self.run_trimal()
        self.run_tree()
        self.set_output()
        self.end()

    def run_align(self):
        if self.num <=200:
            cmd = '{}mafft_accuracy.sh {} {} {} '.format(self.sh, self.mafft, self.seq, self.align)
        else:
            cmd = '{}mafft_fast.sh {} {} {} '.format(self.sh, self.mafft, self.seq, self.align)
        command = self.add_command("run_align", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_align运行完成！")
        else:
            self.set_error("run_align运行完成运行出错!")

    def run_trimal(self):
        cmd = "{} -in {} -out {} -automated1".format(self.trimal, self.align, self.trimal_fa)
        self.logger.info(cmd)
        self.logger.info("开始运行run_trimal")
        command = self.add_command("run_trimal", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_trimal完成")
        else:
            self.set_error("运行run_trimal运行出错!")

    def run_tree(self):
        if os.path.exists(self.work_dir + "/temp2"):
            shutil.rmtree(self.work_dir + "/temp2")
        os.mkdir(self.work_dir + "/temp2")
        cmd = '{} -s {} -pre {} -m MFP -nt AUTO -ntmax 8 -bb {}'.format(self.iqtree, self.trimal_fa, self.work_dir + "/temp2/all", 1000)
        command = self.add_command("run_tree", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_tree运行完成！")
        else:
            self.set_error("run_tree运行完成运行出错!")


    def set_output(self):
        link_file(self.trimal_fa, self.output_dir + "/all.align_last.fna")
        link_file(self.work_dir + "/temp2/all.contree", self.output_dir + "/"+ self.option("sample") + ".house_keeping.nwk")


    def get_num(self, file):
        list= []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if re.search("^>",line):
                    list.append(line)
        return len(list)