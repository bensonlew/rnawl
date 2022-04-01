# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# __modify__ = '2020.06.30'

import os
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file, link_dir


class AssembleSpadesModule(Module):
    """
    单个样品进行三代数据库reads提取，spades的的组装
    """
    def __init__(self, work_id):
        super(AssembleSpadesModule, self).__init__(work_id)
        option = [
            {"name": "data_type", "type": "string"},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            {"name": "base", "type": "string"},
            {"name": "depth_pacbio_num", "type": "int", "default": 200}, # 抽取数据的乘数
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table", "required": True},  # 二代PE reads
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 二代MP reads
        ]
        self.add_option(option)
        self.seqtk = self.add_tool("bacgenome.seqtk")
        self.spades = self.add_tool("assemble.assemble_spades")
        self.step.add_steps('spades')

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(AssembleSpadesModule, self).run()
        self.run_check()
        if self.depth > self.option("depth_pacbio_num"):
            self.spades.on("end", self.end)
            self.seqtk.on("end", self.run_spades)
            self.run_seqtk()
        else:
            self.spades.on("end", self.end)
            self.run_spades()

    def run_check(self):
        self.gsize_in_bp = float(self.option("genomeSize")) * 1024 * 1024
        self.depth = float(self.option("base")) / self.gsize_in_bp


    def run_seqtk(self):
        scale = self.option("depth_pacbio_num") * self.gsize_in_bp / float(self.option("base"))
        self.seqtk.set_options({
            "fastq": self.option("subread_fq").prop["path"],
            "outfastq": self.option("sample_name") + ".pacbio.fq",
            "scale": scale
        })
        self.seqtk.run()

    def run_spades(self):
        if self.depth >= self.option("depth_pacbio_num"):
            subread_fq = self.seqtk.output_dir + "/" + self.option("sample_name") + ".pacbio.fq"
        else:
            subread_fq = self.option("subread_fq").prop["path"]
        opts = {
            "sample_name": self.option("sample_name"),
        }
        if self.option("PE_list").is_set:
            opts['PE_list'] = self.option("PE_list")
        elif self.option("MP_list").is_set:
            opts['MP_list'] = self.option("MP_list")
        if self.option("data_type") in ['pacbio', 'Pacbio']:
            opts['pacbio'] = subread_fq
        elif self.option("data_type") in ['nanopore', 'Nanopore']:
            opts['nanopore'] = subread_fq
        self.spades.set_options(opts)
        self.spades.on("end", self.set_output, "spades")
        self.spades.run()


    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event['data'] == 'spades':
            link_dir(self.spades.output_dir, self.output_dir + "/spades")
        self.logger.info("设置结果成功")

    def end(self):
        super(AssembleSpadesModule, self).end()
