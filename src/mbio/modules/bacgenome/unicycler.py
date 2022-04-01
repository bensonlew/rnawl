# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# __modify__ = '2020.06.30'

import os
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file, link_dir


class UnicyclerModule(Module):
    """
    单个样品进行三代数据库reads提取，unicycler的组装
    """

    def __init__(self, work_id):
        super(UnicyclerModule, self).__init__(work_id)
        option = [
            {"name": "data_type", "type": "string"},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            {"name": "base", "type": "string"},
            {"name": "depth_pacbio_num", "type": "int", "default": 200}, # 抽取数据的乘数
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
        ]
        self.add_option(option)
#        self.seqtk = self.add_tool("bacgenome.seqtk")
        self.unicycler = self.add_tool("assemble.unicycler")
        self.step.add_steps('unicycler')

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(UnicyclerModule, self).run()
#        self.run_check()
#        if self.depth > self.option("depth_pacbio_num"):
#            self.unicycler.on("end", self.end)
#            self.seqtk.on("end", self.run_unicycler)
#            self.run_seqtk()
#        else:
        self.unicycler.on("end", self.end)
        self.run_unicycler()


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

    def run_unicycler(self):
#        if self.depth >= self.option("depth_pacbio_num"):
#            subread_fq = self.seqtk.output_dir + "/" + self.option("sample_name") + ".pacbio.fq"
#        else:
        subread_fq = self.option("subread_fq").prop["path"]
        self.unicycler.set_options({
            "read1": self.option("read1"),
            "read2": self.option("read2"),
            "long": subread_fq,
            "sample_name": self.option("sample_name"),
        })
        self.unicycler.on("end", self.set_output, "unicycler")
        self.unicycler.run()

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event['data'] == 'unicycler':
            link_dir(self.unicycler.output_dir, self.output_dir + "/Unicycler")
        self.logger.info("设置结果成功")

    def end(self):
        super(UnicyclerModule, self).end()
