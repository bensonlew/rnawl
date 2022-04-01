# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '20200314'
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir
import pandas as pd


class BacAssem2Module(Module):
    """
    细菌基因组流程组装总模块，根据总的样品信息提取有用信息进行分析
    """

    def __init__(self, work_id):
        super(BacAssem2Module, self).__init__(work_id)
        option = [
            {"name": "sample_list", "type": "string"},
            {"name": "kmers", "type": "string"},
            {"name": "error_rate", "type": "string", "default": "0.025"},  # canu拼接的错误率参数
            {"name": "cor_min_coverage", "type": "string", "default": "default", 'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "cor_mhap_sensitivity", "type": "string", "default": "default", "choose": ["default", "low", "normal", "high"]}, # canu param
            # 其他的拼接参数一大堆......
            {'name': 'depth_pacbio_num', 'type': 'int', 'default': 200},  # 抽取序列的覆盖度
            {"name": "pe_list", "type": "string"},  # 哪些样本是做二代数据
            {"name": "third_list", "type": "string"},  # 哪些样本是做二代数据
            {"name": "extract_reads_dir", "type": "infile", "format": "sequence.fastq_dir"}, #二代质控后提取的reads
            {"name": "third_dir", "type": "string"}, #三代质控fq格式的reads
            {'name': 'pe_assem_tool', 'type': 'string', 'default': 'soapdenovo', 'choose': ['soapdenovo', 'velvet']},
            {"name": "third_stat", "type": "string"},
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
        self.modules = []
        self.samples = []
        self.dict ={}

    def run(self):
        super(BacAssem2Module, self).run()
        self.logger.info(self.option("pe_list"))
        if self.option("pe_list"):
            self.pe_dict, self.pe_fq_dict = self.get_pe_info(self.option("pe_list"))
        else:
            pass
        if self.option("third_list"):
            self.third_dict, self.third_fq_dict, self.third_size_dict = self.get_third_info(self.option("third_list"))
            self.base_dict = self.get_third_base(self.option("third_stat"))
        else:
            pass
        self.run_assemble()


    def run_assemble(self):
        self.logger.info(self.option("sample_list"))
        if re.search(",", self.option("sample_list")):
            self.samples = self.option("sample_list").split(",")
        else:
            self.samples.append(self.option("sample_list"))
        self.logger.info(self.samples)
        self.logger.info(self.pe_dict)
        for sample in self.samples:
            assem = self.add_module("bacgenome.bac_assem")
            self.dict[assem] = sample
            if self.option("pe_list") and not self.option("third_list"):
                library_type = self.pe_dict[sample]
                assem.set_options({
                    "library_type": library_type,
                    "fq_dir": self.option("extract_reads_dir"),
                    "sample": sample,
                    "pe_assem_tool": self.option("pe_assem_tool"),
                    "sample_info": self.option("pe_list"),
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
            elif not self.option("pe_list") and self.option("third_list"):
                library_type = self.third_dict[sample]
                assem.set_options({
                    "data_type": self.third_dict[sample],
                    "library_type": library_type,
                    "subread_fq": self.third_fq_dict[sample],
                    "genomeSize": self.third_size_dict[sample],
                    "sample": sample,
                    "base": self.base_dict[sample],
                    "depth_pacbio_num": self.option("depth_pacbio_num"),
                    "error_rate": self.option("error_rate"),
                    "cor_min_coverage": self.option("cor_min_coverage"),
                    "cor_mhap_sensitivity": self.option("cor_mhap_sensitivity"),
                })
            else:
                if sample in self.pe_dict and sample not in self.third_dict:
                    library_type = self.pe_dict[sample]
                    assem.set_options({
                        "library_type": library_type,
                        "fq_dir": self.option("extract_reads_dir"),
                        "sample": sample,
                        "pe_assem_tool": self.option("pe_assem_tool"),
                        "sample_info": self.option("pe_list"),
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
                elif not sample in self.pe_dict and sample in self.third_dict:
                    library_type = self.third_dict[sample]
                    assem.set_options({
                        "data_type": self.third_dict[sample],
                        "library_type": library_type,
                        "subread_fq": self.third_fq_dict[sample],
                        "genomeSize": self.third_size_dict[sample],
                        "sample": sample,
                        "base": self.base_dict[sample],
                        "depth_pacbio_num": self.option("depth_pacbio_num"),
                        "error_rate": self.option("error_rate"),
                        "cor_min_coverage": self.option("cor_min_coverage"),
                        "cor_mhap_sensitivity": self.option("cor_mhap_sensitivity"),
                    })

                elif sample in self.pe_dict and sample in self.third_dict:
                    library_type = self.pe_dict[sample] + "+" + self.third_dict[sample]
                    self.get_pe_list(sample, self.pe_fq_dict, self.work_dir + "/" + sample + ".pe_list.txt")
                    assem.set_options({
                        "data_type": self.third_dict[sample],
                        "library_type": library_type,
                        "fq_dir": self.option("extract_reads_dir"),
                        "sample": sample,
                        "pe_assem_tool": self.option("pe_assem_tool"),
                        "sample_info": self.option("pe_list"),
                        "kmers": self.option("kmers"),
                        "soapdenovo_D": self.option("soapdenovo_D"),
                        "soapdenovo_d": self.option("soapdenovo_d"),
                        "soapdenovo_M": self.option("soapdenovo_M"),
                        "soapdenovo_R": self.option("soapdenovo_R"),
                        "soapdenovo_F": self.option("soapdenovo_F"),
                        "soapdenovo_u": self.option("soapdenovo_u"),
                        "soapdenovo_G": self.option("soapdenovo_G"),
                        "velvet_min_contig_lgth": self.option("velvet_min_contig_lgth"),
                        "velvet_min_pair_count": self.option("velvet_min_pair_count"),
                        "subread_fq": self.third_fq_dict[sample],
                        "genomeSize": self.third_size_dict[sample],
                        "base": self.base_dict[sample],
                        "depth_pacbio_num": self.option("depth_pacbio_num"),
                        "error_rate": self.option("error_rate"),
                        "cor_min_coverage": self.option("cor_min_coverage"),
                        "cor_mhap_sensitivity": self.option("cor_mhap_sensitivity"),
                        "read1": self.pe_fq_dict[sample].split(",")[0],
                        "read2": self.pe_fq_dict[sample].split(",")[1],
                        "PE_list": self.work_dir + "/" + sample + ".pe_list.txt"
                    })
            self.modules.append(assem)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on('end', self.set_output)
        for module in self.modules:
            module.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        if len(self.modules) > 1:
            for module in self.modules:
                if os.path.exists(module.output_dir + "/pacbio_assess"):
                    link_dir(module.output_dir + "/pacbio_assess", self.output_dir + "/pacbio_assess")
                link_dir(module.output_dir + "/" + self.dict[module], self.output_dir + "/assemble/" + self.dict[module])
        elif len(self.modules) == 1:
            if os.path.exists(self.modules[0].output_dir + "/pacbio_assess"):
                link_dir(self.modules[0].output_dir + "/pacbio_assess", self.output_dir + "/pacbio_assess")
            link_dir(self.modules[0].output_dir + "/" + self.dict[self.modules[0]], self.output_dir + "/assemble/" + self.dict[self.modules[0]])
        self.end()

    def end(self):
        super(BacAssem2Module, self).end()

    def get_pe_info(self, input):
        dict ={}
        dict_fq = {}
        with open(input, "r") as f:
            lines =f.readlines()
            for line in lines[1:]:
                line = line.rstrip().split("\t")
                dict[line[5]] = line[6]
                fq = line[-2].split(",")
                dict_fq[line[5]] = self.option("extract_reads_dir").prop['path'] + "/" + fq[0] + "," + self.option("extract_reads_dir").prop['path'] + "/" + fq[1]
        return dict, dict_fq

    def get_third_info(self, input):
        dict ={}
        dict_fq ={}
        dict_size = {}
        with open(input, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip().split("\t")
                dict[line[0]] = line[1]
                dict_fq[line[0]] = self.option("third_dir") + "/" + line[2]
                dict_size[line[0]] = line[3]
        return dict, dict_fq, dict_size

    def get_third_base(self, input):
        dict ={}
        with open(input, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.rstrip().split("\t")
                dict[line[0]] = line[3]
        return dict

    def get_pe_list(self,sample, dict1, output):
        with open(output, "w") as f:
            if sample in dict1:
                fq = dict1[sample].split(",")
                f.write("{}\t{}\t{}\n".format(sample, fq[0], fq[1]))