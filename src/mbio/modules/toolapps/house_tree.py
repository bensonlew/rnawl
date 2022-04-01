# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.19

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class HouseTreeModule(Module):
    """
    从NCBI下载构建看家基因
    """
    def __init__(self, work_id):
        super(HouseTreeModule, self).__init__(work_id)
        options = [
            {"name": "input", "type": "infile", "format": "sequence.fasta"},  # 样品的基因组的看家基因序列
            {"name": "genomes", "type": "infile", "format": "sequence.fasta_dir"},#数据库的参考基因组文件夹
            {"name": "sample", "type": "string"},  ##样品名称
        ]
        self.add_option(options)
        self.genomes_house = self.add_module('toolapps.genomes_house')
        self.house_tree = self.add_tool('toolapps.house_tree')


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("input").is_set:
            raise OptionError("必须设置参数input!")
        if not self.option("genomes").is_set:
            raise OptionError("必须设置参数genomes!")

    def run_housekeeping(self):

        opts = {
            'genomes': self.option('genomes'),
        }
        self.genomes_house.set_options(opts)
        self.genomes_house.on('end', self.run_house_tree)
        self.genomes_house.run()

    def run_house_tree(self):
        """
        house_keeping进化树
       :return:
        """
        os.system("cat {} {} >{}".format(self.option("input").prop['path'], self.genomes_house.option("out").prop['path'], self.work_dir + "/all.house_keeping.fasta"))
        opts = {
            'seq': self.work_dir + "/all.house_keeping.fasta",
            "sample": self.option("sample"),
        }
        self.house_tree.set_options(opts)
        self.house_tree.on("end", self.set_output)
        self.house_tree.run()

    def run(self):
        """
        运行
        :return:
        """
        super(HouseTreeModule, self).run()
        self.run_housekeeping()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.output_dir+"/"+ self.option("sample") +".house_keeping.nwk"):
            os.remove(self.output_dir+"/"+ self.option("sample") +".house_keeping.nwk")
        self.logger.info(self.house_tree.output_dir+"/"+ self.option("sample") +".house_keeping.nwk")
        os.link(self.house_tree.output_dir+"/"+ self.option("sample") +".house_keeping.nwk", self.output_dir+"/"+ self.option("sample") +".house_keeping.nwk")
        self.end()

    def end(self):
        super(HouseTreeModule, self).end()