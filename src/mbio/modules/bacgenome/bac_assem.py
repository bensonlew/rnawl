# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '20200314'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir,link_file
import pandas as pd


class BacAssemModule(Module):
    """
    细菌基因组流程组装总模块，根据总的样品sample_info分别进行三种情况的组装，PE，PE+Pacbio,Pacbio
    """

    def __init__(self, work_id):
        super(BacAssemModule, self).__init__(work_id)
        option = [
            {"name": "data_type", "type": "string"},
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "sample", "type": "string"},
            {"name": "kmers", "type": "string"},
            {'name': 'depth_pacbio_num', 'type': 'int', 'default': 200},
            {'name': 'pe_assem_tool', 'type': 'string', 'default': 'soapdenovo','choose': ['soapdenovo', 'velvet']},  # 二代数据拼接工具
            {"name": "library_type", "type": "string"},
            {"name": "sample_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},
            {"name": "subread_fq", "type": "infile", "format": "sequence.fastq"},  ##三代数据subreads.fq文件
            {"name": "genomeSize", "type": "string"},  ##基因组大小
            {"name": "base", "type": "string"},
            {"name": "error_rate", "type": "string", "default": "0.025"},  # canu拼接的错误率参数
            {"name": "cor_min_coverage", "type": "string", "default": "default",
             'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "cor_mhap_sensitivity", "type": "string", "default": "default",
             "choose": ["default", "low", "normal", "high"]}  ,# canu param
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table", "required": True},  # 二代PE reads
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 二代MP reads
        ]
        soapdenovo_opt = [
            {"name": "soapdenovo_D", "type": "int", "default": 1, "min": 1, "max": 10},
            {"name": "soapdenovo_d", "type": "string", "default": "3,5,10"},
            {"name": "soapdenovo_M", "type": "int", "default": 1, "min":0, "max": 3},
            {"name": "soapdenovo_R", "type": "bool", "default": True},
            {"name": "soapdenovo_F", "type": "bool", "default": True},
            {"name": "soapdenovo_u", "type": "string", "default": "unmask", "choose": ["mask", "unmask"]},
            {"name": "soapdenovo_G", "type": "int", "default": 50, "min": 0}
        ]
        velvet_opt = [
            {"name": "velvet_min_contig_lgth", "type": "int", "default": 200, "min": 100},
            {"name": "velvet_min_pair_count", "type": "int", "default": 15, "min": 5}
        ]
        option += soapdenovo_opt + velvet_opt
        self.add_option(option)
        self.assem_tools = []
        self.step.add_steps('soapdenovo', 'canu', 'spades', 'unicycler')

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(BacAssemModule, self).run()
        if self.option("library_type") in ['PE']:
            self.pe_assem = self.add_module("bacgenome.assemble_hiseq")
            self.run_pe_assemble()
        elif self.option("library_type") in ['Pacbio',"Nanopore"]:
            self.canu_assem = self.add_module("bacgenome.canu_assem")
            self.run_canu_assemble()
        elif self.option("library_type") in ['PE+Pacbio', "PE+Nanopore", "Nanopore+PE", "Pacbio+PE"]:
            self.pe_assem = self.add_module("bacgenome.assemble_hiseq")
            self.canu_assem = self.add_module("bacgenome.canu_assem")
            self.spades = self.add_module("bacgenome.assemble_spades")
            self.unicycler = self.add_module("bacgenome.unicycler")
            self.on_rely([self.pe_assem,self.canu_assem,self.unicycler,self.spades], self.end)
            self.run_pe_assemble()
            self.run_canu_assemble()
            self.run_spades_assemble()
            self.run_unicycler_assemble()

    def run_pe_assemble(self):
        self.pe_assem.set_options({
            "fq_dir": self.option("fq_dir"),
            "sample_name": self.option('sample'),
            "assem_tool": self.option("pe_assem_tool"),
            "sample_info": self.option("sample_info"),
            "kmers": self.option("kmers"),
            "soapdenovo_D": self.option("soapdenovo_D"),
            "soapdenovo_d": self.option("soapdenovo_d"),
            "soapdenovo_M": self.option("soapdenovo_M"),
            "soapdenovo_R": self.option("soapdenovo_R"),
            "soapdenovo_F": self.option("soapdenovo_F"),
            "soapdenovo_u": self.option("soapdenovo_u"),
            "soapdenovo_G": self.option("soapdenovo_G"),
            "velvet_min_contig_lgth": self.option("velvet_min_contig_lgth"),
            "velvet_min_pair_count": self.option("velvet_min_pair_count")
        })
        self.pe_assem.on("end", self.set_output, "soapdenovo")
        self.pe_assem.run()

    def run_canu_assemble(self):
        self.canu_assem.set_options({
            "data_type":self.option("data_type"),
            "base": self.option("base"),
            "depth_pacbio_num":self.option("depth_pacbio_num"),
            "subread_fq": self.option("subread_fq"),
            "genomeSize": self.option("genomeSize"),
            "sample_name": self.option("sample"),
            "error_rate": self.option("error_rate"),
            "cor_min_coverage": self.option("cor_min_coverage"),
            "cor_mhap_sensitivity": self.option("cor_mhap_sensitivity"),
        })
        self.canu_assem.on("end", self.set_output, "canu")
        self.canu_assem.run()

    def run_spades_assemble(self):
        opts = {
            "data_type": self.option("data_type"),
            "base": self.option("base"),
            "depth_pacbio_num": self.option("depth_pacbio_num"),
            "sample_name": self.option("sample"),
            "genomeSize": self.option("genomeSize"),
            "subread_fq": self.option("subread_fq"),
        }
        if self.option("PE_list").is_set:
            opts['PE_list'] = self.option("PE_list")
        elif self.option("MP_list").is_set:
            opts['MP_list'] = self.option("MP_list")
        self.spades.set_options(opts)
        self.spades.on("end", self.set_output, "spades")
        self.spades.run()

    def run_unicycler_assemble(self):
        self.unicycler.set_options({
            "data_type": self.option("data_type"),
            "base": self.option("base"),
            "depth_pacbio_num": self.option("depth_pacbio_num"),
            "genomeSize": self.option("genomeSize"),
            "read1": self.option("read1"),
            "read2": self.option("read2"),
            "subread_fq": self.option("subread_fq"),
            "sample_name": self.option("sample"),
        })
        self.unicycler.on("end", self.set_output, "unicycler")
        self.unicycler.run()

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        if self.option("library_type") in ['PE']:
            if event['data'] == 'soapdenovo':
                link_dir(self.pe_assem.output_dir, self.output_dir + "/" + self.option("sample") + "/soapdenovo")
            self.end()
        elif self.option("library_type") in ['Pacbio',"Nanopore"]:
            if event['data'] == 'canu':
                link_dir(self.canu_assem.output_dir, self.output_dir + "/" + self.option("sample") +  "/canu" )
            self.end()
        elif self.option("library_type") in ['PE+Pacbio', "PE+Nanopore", "Nanopore+PE", "Pacbio+PE"]:
            if event['data'] == 'soapdenovo':
                link_dir(self.pe_assem.output_dir, self.output_dir + "/" + self.option("sample") +  "/soapdenovo")
            if event['data'] == 'canu':
                link_dir(self.canu_assem.output_dir + "/Canu", self.output_dir + "/" + self.option("sample") + "/canu")
                link_dir(self.canu_assem.output_dir + "/pacbio_assess", self.output_dir + "/pacbio_assess")
                link_file(self.canu_assem.output_dir + "/" + self.option("sample") + "_pacbio.clean.fa", self.output_dir + "/pacbio_assess/" + self.option("sample") + ".pacbio_clean.fa")
            if event['data'] == 'spades':
                link_dir(self.spades.output_dir + "/spades", self.output_dir + "/" + self.option("sample") + "/spades")
            if event['data'] == 'unicycler':
                link_dir(self.unicycler.output_dir + "/Unicycler", self.output_dir + "/" + self.option("sample") + "/unicycler")

    def end(self):
        super(BacAssemModule, self).end()