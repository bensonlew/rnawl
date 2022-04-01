# -*- coding: utf-8 -*-
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class InputProcessModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(InputProcessModule, self).__init__(work_id)
        options = [
            {"name": "input_dir", "type": "infile", "format": "sequence.fastq_dir2"},
            {"name": "sample_info", "type": "infile", "format": "meta_genomic.specimen_info"},
            {"name": "qc_quality", "type": "int", "default": 20},
            {"name": "qc_length", "type": "int", "default": 50},
            {"name": "qc", "type": "bool", "default": False},
        ]
        self.modules = []
        self.add_option(options)
        self.cat_gz = self.add_tool("sequence.cat_file")
        self.fastp = self.add_module("meta.qc.fastp_qc")
        self.fp_stat = self.add_tool("metagenomic.fastp_stat")
        self.sample_path = defaultdict(list)
        self.fq_path = ""
        self.stat_path = ""

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('input_dir').is_set:
            raise OptionError('必须输入input_dir文件')

    def run_cat_gz(self):
        self.cat_gz.set_options({'input_dir': self.option("input_dir")})
        self.cat_gz.run()

    def run_qc(self):
        opts = {
            'fastq_dir': self.cat_gz.output_dir,
            'qualified_quality_phred': str(self.option('qc_quality')),
            'length_required': str(self.option('qc_length')),
            'length_required': str(self.option('qc_length')),
            'compression': '5',
        }
        self.fastp.set_options(opts)
        self.fastp.run()

    def get_stat(self):
        self.fp_stat.set_options(
            {"json_dir": self.fastp.output_dir + "/qc_stat",
             "sample_info": self.option("sample_info")})
        self.fp_stat.run()

    def run(self):
        """
        运行
        :return:
        """
        super(InputProcessModule, self).run()
        self.cat_gz.on("end", self.run_qc)
        self.fastp.on("end", self.get_stat)
        self.fp_stat.on("end", self.set_output)
        self.run_cat_gz()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """

        if self.option("qc"):
            self.fq_path = self.fastp.output_dir + "/after_qc_dir"
        else:
            self.fq_path = self.cat_gz.output_dir
        if not os.path.exists(self.output_dir + "/rawdata"):
            os.mkdir(self.output_dir + "/rawdata")
        self.qc_stat = self.fp_stat.output_dir + "/reads.cleanData.stat"
        self.link(self.fp_stat.output_dir + "/reads.rawData.stat", "output/rawdata/")
        self.link(self.fp_stat.output_dir + "/base_info.txt", "output/rawdata/")
        self.raw_stat = self.fp_stat.work_dir + "/reads.rawData2.stat"
        self.rawdata = self.output_dir + "/rawdata"
        self.base_info = self.fp_stat.output_dir + 'base_info.txt'
        self.end()

    def end(self):
        super(InputProcessModule, self).end()
