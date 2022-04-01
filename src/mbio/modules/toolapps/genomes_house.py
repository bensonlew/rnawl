# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.14

from biocluster.module import Module
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class GenomesHouseModule(Module):
    """
    基于参考基因组和查询基因组，得到看家基因的连起来的序列
    """
    def __init__(self, work_id):
        super(GenomesHouseModule, self).__init__(work_id)
        options = [
            {"name": "genomes", "type": "outfile", "format": "sequence.fasta_dir"},
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)
        self.modules =[]


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genomes").is_set:
            raise OptionError("必须设置参数genomes的目录！")

    def run_housekeeping(self):
        for file in os.listdir(self.option("genomes").prop['path']):
            sample = file.split(".fasta")[0]
            self.housekeeping = self.add_tool('toolapps.genome_housekeeping')
            opts = {
                'input_genome': self.option("genomes").prop['path'] + "/" + file,
                'sample': sample,
            }
            self.housekeeping.set_options(opts)
            self.modules.append(self.housekeeping)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(GenomesHouseModule, self).run()
        self.run_housekeeping()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.work_dir+"/house_keeping"):
            shutil.rmtree(self.work_dir+"/house_keeping")
        os.mkdir(self.work_dir+"/house_keeping")
        for module in self.modules:
            link_dir(module.output_dir, self.work_dir+"/house_keeping")
        files = os.listdir(self.work_dir+"/house_keeping")
        de_list =[]
        for file in files:
            if re.search(".core_gene.fa", file):
                de_list.append(self.work_dir+"/house_keeping/"+file)
        os.system("cat {} >{}".format(" ".join(de_list), self.work_dir+"/all.house_keeping.fa"))
        if os.path.exists(self.output_dir+"/all.house_keeping.fa"):
            os.remove(self.output_dir+"/all.house_keeping.fa")
        os.link(self.work_dir+"/all.house_keeping.fa", self.output_dir+"/all.house_keeping.fa")
        self.option("out", self.output_dir+"/all.house_keeping.fa")
        self.end()

    def end(self):
        super(GenomesHouseModule, self).end()