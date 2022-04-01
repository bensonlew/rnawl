# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# __modify__ = '2021.05.19'

import os
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file, link_dir


class Unicycler2Module(Module):
    """
    单个样品进行三代数据库reads提取，unicycler的组装
    """

    def __init__(self, work_id):
        super(Unicycler2Module, self).__init__(work_id)
        option = [
            {"name": "data_type", "type": "string"},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "sample_name", "type": "string"},
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
        ]
        self.add_option(option)
#        self.bam_to_fq = self.add_tool("bacgenome.pacbio_convert")
#        self.seqtk = self.add_tool("bacgenome.seqtk")
        self.unicycler = self.add_tool("bacgenome.unicycler")

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("data_type"):
            raise OptionError("请提供三代数据类型！nanopore or pacbio")
        return True

    def run(self):
        super(Unicycler2Module, self).run()
#        self.num = float(os.path.getsize(self.option("subread_fq").prop["path"])) / (1024 * 1024 * 1024)
#        self.logger.info(self.num)
#        if self.num > 3:
#            self.seqtk.on("end", self.run_unicycler)
#            self.run_seqtk()
#        else:
        self.run_unicycler()



    def run_seqtk(self):
        fq = self.option("subread_fq").prop["path"]
        scal = 3 /self.num
        if self.option("data_type") in['nanopore']:
            self.seqtk.set_options({
                "fastq": fq,
                "outfastq": self.option("sample_name") + ".nanopore.fq",
                "scale": scal
            })
        elif self.option("data_type") in['pacbio']:
            self.seqtk.set_options({
                "fastq": fq,
                "outfastq": self.option("sample_name") + ".pacbio.fq",
                "scale": scal
            })
        self.seqtk.run()


    def run_unicycler(self):
#        subread_fq = ''
#        if self.option("data_type") in ['nanopore']:
#            if self.num > 3:
#                subread_fq = self.seqtk.output_dir + "/" + self.option("sample_name") + ".nanopore.fq"
#            else:
#            subread_fq = self.option("subread_fq").prop["path"]
#        elif self.option("data_type") in ['pacbio']:
#            if self.num > 3:
#                subread_fq = self.seqtk.output_dir + "/" + self.option("sample_name") + ".pacbio.fq"
#            else:
        subread_fq = self.option("subread_fq").prop["path"]
        self.unicycler.set_options({
            "read1": self.option("read1"),
            "read2": self.option("read2"),
            "long": subread_fq,
            "sample_name": self.option("sample_name"),
        })
        self.unicycler.on("end", self.set_output)
        self.unicycler.run()

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        link_dir(self.unicycler.output_dir, self.output_dir + "/Unicycler")
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        super(Unicycler2Module, self).end()