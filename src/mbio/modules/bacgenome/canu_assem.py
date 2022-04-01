# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# __modify__ = '2020.03.30'

import os
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file, link_dir


class CanuAssemModule(Module):
    """
    单个样品进行三代数据库reads提取，canu的reads质控，canu的组装
    """

    def __init__(self, work_id):
        super(CanuAssemModule, self).__init__(work_id)
        option = [
            {"name": "data_type", "type": "string"},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "sample_name", "type": "string"},
            {"name": "base", "type": "string"},
            {"name": "depth_pacbio_num", "type": "int", "default": 200}, # 抽取数据的乘数
            {"name": "error_rate", "type": "string", "default": "0.025"},  # canu拼接的错误率参数
            {"name": "cor_min_coverage", "type": "string", "default": "default",
             'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "cor_mhap_sensitivity", "type": "string", "default": "default",
             "choose": ["default", "low", "normal", "high"]}  # canu param
        ]
        self.add_option(option)
        self.seqtk = self.add_tool("bacgenome.seqtk")
        self.canu = self.add_tool("assemble.canu_assem")
        self.gunzip_fastq = self.add_tool('bacgenome.fastq_ungz')
        self.qc_third = self.add_tool("bacgenome.pacbio_num")  # 对三代数据进行长度统计
        self.step.add_steps('pacbio_stat', 'canu')

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(CanuAssemModule, self).run()
        self.run_check()
        if self.depth > self.option("depth_pacbio_num"):
            self.qc_third.on("end", self.end)
            self.gunzip_fastq.on("end", self.run_pacbio_assess)
            self.canu.on("end", self.run_gizp_pacbio)
            self.seqtk.on("end", self.run_canu)
            self.run_seqtk()
        else:
            self.qc_third.on("end", self.end)
            self.gunzip_fastq.on("end", self.run_pacbio_assess)
            self.canu.on("end", self.run_gizp_pacbio)
            self.run_canu()

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

    def run_canu(self):
        if self.depth >= self.option("depth_pacbio_num"):
            subread_fq = self.seqtk.output_dir + "/" + self.option("sample_name") + ".pacbio.fq"
        else:
            subread_fq = self.option("subread_fq").prop["path"]
        self.canu.set_options({
            "data_type": self.option("data_type"),
            "subread_fq": subread_fq,
            "genomeSize": str(self.option("genomeSize")) + 'm',
            "sample_name": self.option("sample_name"),
            "error_rate": self.option("error_rate"),
            "cor_min_coverage": self.option("cor_min_coverage"),
            "cor_mhap_sensitivity": self.option("cor_mhap_sensitivity"),
        })
        self.canu.on("end", self.set_output, "canu")
        self.canu.run()

    def run_gizp_pacbio(self):
        self.gunzip_fastq.set_options({
            "fastq": self.canu.work_dir + "/assemble/" + self.option("sample_name") + ".trimmedReads.fasta.gz",
            "sample_name": self.option("sample_name"),
            "direction": '1',
            "lib_type": 'pacbio',
            "result_path": self.gunzip_fastq.work_dir,
            "nozip": "false"
        })
        self.gunzip_fastq.run()

    def run_pacbio_assess(self):
        if os.path.exists(self.output_dir + "/" + self.option("sample_name") + "_pacbio.clean.fa"):
            os.remove(self.output_dir + "/" + self.option("sample_name") + "_pacbio.clean.fa")
        os.link(self.gunzip_fastq.work_dir + "/" + self.option("sample_name") + "_pacbio.1.fq", self.output_dir + "/" + self.option("sample_name") + "_pacbio.clean.fa")
        self.qc_third.set_options({
            "input_fa": self.output_dir + "/" + self.option("sample_name") + "_pacbio.clean.fa",
            "sample_name": self.option("sample_name")
        })
        self.qc_third.on("end", self.set_output, "pacbio_stat")
        self.qc_third.run()

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event['data'] == 'canu':
            link_dir(self.canu.output_dir, self.output_dir + "/Canu")
        if event['data'] == 'pacbio_stat':
            link_dir(self.qc_third.output_dir, self.output_dir + "/pacbio_assess")
        self.logger.info("设置结果成功")

    def end(self):
        super(CanuAssemModule, self).end()
