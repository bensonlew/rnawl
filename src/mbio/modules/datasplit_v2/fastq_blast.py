# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171215

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class FastqBlastModule(Module):
    """
    对单端或者双端的fastq进行随机抽取10000条序列，
    根据参数database选择比对到Nt/Rfam数据库，对比对的xml结果进行Nt/Rfam分类统计
    """
    def __init__(self, work_id):
        super(FastqBlastModule, self).__init__(work_id)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # fastq文件
            {"name": "fastq_l", "type": "infile", "format": "datasplit.fastq"},  # 样本左端fastq文件
            {"name": "fastq_r", "type": "infile", "format": "datasplit.fastq"},  # 样本右端fastq文件
            {"name": "num", "type": "string", "default": "5000"},  # 总共抽取出来的序列数
            {"name": "database", "type": "string"},  # 数据库，nt或者rfam
            {"name": "evalue", "type": "float", "default": 0.01},  # evalue值
            {"name": "num_threads", "type": "int", "default": 2},  # cpu数
            {"name": "memory", "type": "int", "default": 25},  # 内存G
            {"name": "num_alignment", "type": "int", "default": 5},
        ]
        self.add_option(options)
        self.fastq_extract = self.add_tool("datasplit_v2.sequence_extract")

    def check_options(self):
        if self.option("fastq").is_set:
            pass
        elif self.option("fastq_l").is_set and self.option("fastq_r").is_set:
            pass
        else:
            raise OptionError("请设置fastq文件或者fastq_l和fastq_r文件")
        if not self.option("database"):
            raise OptionError("请设置比对的数据库nt/rfam")
        if self.option("database").lower() not in ["nt", "rfam", "nt_v20200604", "nt_v202107","nt_v20211009", "rfam_v14.6"]:
            raise OptionError("数据库{}不在nt、rfam内".format(self.option("database")))

    def run_fastq_extract(self):
        if self.option("fastq").is_set:
            options = {
                "fastq": self.option("fastq"),
                "num": self.option("num")
            }
        else:
            options = {
                "fastq_l": self.option("fastq_l"),
                "fastq_r": self.option("fastq_r"),
                "num": self.option("num")
            }
        self.fastq_extract.set_options(options)
        self.fastq_extract.on("end", self.run_blast)
        self.fastq_extract.run()

    def run_blast(self):
        options = {
            "query": self.fastq_extract.option("out_fasta"),
            "database": self.option("database"),
            "query_type": "nucl",
            "reference_type": "nucl",
            "blast": "blastn",
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "memory": self.option("memory"),
            "num_alignment": self.option("num_alignment"),
        }
        if self.option("database") in ["nt_v20211009", "rfam_v14.6"] :
            self.blast = self.add_tool("align.blast_nt")
        else:
            self.blast = self.add_tool("align.blast") 
        self.blast.set_options(options)
        if self.option("database").startswith("nt"):
            self.blast.on("end", self.run_nt_stat)
        else:
            self.blast.on("end", self.run_rfam_stat)
        self.blast.run()

    def run_nt_stat(self):
        options = {
            "nt_xml": self.blast.option("outxml"),
            "query": self.fastq_extract.option("out_fasta")
        }
        self.nt_stat = self.add_tool("datasplit_v2.blast_nt_stat")
        self.nt_stat.set_options(options)
        self.nt_stat.on('end', self.set_output, "nt_stat")
        self.nt_stat.run()

    def run_rfam_stat(self):
        options = {
            "rfam_xml": self.blast.option("outxml")
        }
        self.rfam_stat = self.add_tool("datasplit_v2.blast_rfam_stat")
        self.rfam_stat.set_options(options)
        self.rfam_stat.on('end', self.set_output, "rfam_stat")
        self.rfam_stat.run()

    def set_output(self, event):
        self.logger.info("set output nt/rfam")
        if event["data"] == "nt_stat":
            f1 = self.nt_stat.output_dir + "/nt_species_stat.xls"
            f2 = self.output_dir + "/nt_species_stat.xls"
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
            self.end()
        if event["data"] == "rfam_stat":
            f1 = self.rfam_stat.output_dir + "/rfam_summary.xls"
            f2 = self.output_dir + "/rfam_summary.xls"
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
            self.end()
        # self.end()

    def run(self):
        super(FastqBlastModule, self).run()
        self.run_fastq_extract()

    def end(self):
        super(FastqBlastModule, self).end()

