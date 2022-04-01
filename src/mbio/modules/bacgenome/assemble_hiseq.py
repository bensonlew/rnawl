# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/15'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir
import pandas as pd


class AssembleHiseqModule(Module):
    """
    细菌基因组流程做二代数据拼接的模块，进行单样本的组装
    """

    def __init__(self, work_id):
        super(AssembleHiseqModule, self).__init__(work_id)
        option = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "sample_name", "type": "string", "required": True},  # 样本名
            {"name": "assem_tool", "type": "string", "choose": ['soapdenovo', 'velvet']},  # 用于本样本拼接的软件
            {"name": "sample_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},  # 样品相关信息,并非1.0的sample_info
            {"name": "kmers", "type": "string"}
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
        self.mp_cut = self.add_module("bacgenome.mp_cut_fq")
        self.config_file = self.add_tool("assemble.create_config")
        self.assemble = ""  # 根据具体参数来添加
        self.scaf_select = self.add_tool("assemble.scaf_select")
        self.gapcloser = self.add_tool("assemble.gapcloser_scaf")
        self.pe_list = os.path.join(self.work_dir, "pe_list")
        self.sample_info = os.path.join(self.work_dir, "sample_info")
        self.has_mp = True

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(AssembleHiseqModule, self).run()
        self.run_sample_info()
        if self.has_mp:
            self.run_mp_cut()
        elif self.option("assem_tool") == "spades":
            self.run_assemble()
        else:
            self.run_config()

    def run_sample_info(self):
        # 生成老模块需要的sample_info
        # self.pe_list, self.sample_info, self.has_mp
        data = pd.read_table(self.option("sample_info").prop["path"])
        data["Sample Name"] = data["Sample Name"].apply(str)
        data = data[data["Sample Name"]==self.option("sample_name")]
        if "MP" in data["Library"].tolist():
            self.has_mp = True
        else:
            self.has_mp = False
        pe_data = data[data["Library"]=="PE"]
        pe_list = pe_data[["Sample_lib", "Insert Size(bp)", "Read Len(bp)"]]
        pe_list.loc[:, "fq1"] = self.option("fq_dir").prop["path"] + "/" + pe_list["Sample_lib"] + ".clean.1.fq"
        pe_list.loc[:, "fq2"] = self.option("fq_dir").prop["path"] + "/" + pe_list["Sample_lib"] + ".clean.2.fq"
        pe_list.loc[:, "fqs"] = ""  # 为了兼容老的config程序
        pe_list.reindex(columns=["Sample_lib", "fq1", "fq2", "fqs", "Insert Size(bp)", "Read Len(bp)"]).\
            to_csv(self.pe_list, sep="\t", header=False, index=False)
        data[["Sample_lib", "Insert Size(bp)", "Read Len(bp)"]].to_csv(self.sample_info, sep="\t", header=False, index=False)


    def run_mp_cut(self):
        # 截短MP文库成50bp的reads
        self.mp_cut.set_options({
            "fq_dir": self.option("fq_dir"),
            "sample_info": self.sample_info,
            "sample_name": self.option("sample_name")
        })
        if self.option("assem_tool") == "soapdenovo":
            self.mp_cut.on('end', self.run_config)
        else:
            self.mp_cut.on("end", self.run_assemble)
        self.mp_cut.run()

    def run_config(self):
        self.config_file.set_options({
            "PE_list": self.pe_list,
            "MP_list": self.mp_cut.option("mp_list"),
            "sample_name": self.option("sample_name")
        })
        self.config_file.on('end', self.run_assemble)
        self.config_file.run()

    def run_assemble(self):
        if self.option("assem_tool") == "soapdenovo":
            self.assemble = self.add_module("bacgenome.assemble_soap_denovo")
            self.assemble.set_options({
                "config": self.config_file.option("config_file"),
                "sample_name": self.option("sample_name"),
                "kmers": self.option("kmers"),
                "D": self.option("soapdenovo_D"),
                "d": self.option("soapdenovo_d"),
                "M": self.option("soapdenovo_M"),
                "R": self.option("soapdenovo_R"),
                "F": self.option("soapdenovo_F"),
                "u": self.option("soapdenovo_u"),
                "G": self.option("soapdenovo_G")
            })
            self.assemble.on("end", self.run_select)
        elif self.option("assem_tool") == "velvet":
            self.assemble = self.add_module("bacgenome.assemble_velvet")
            self.assemble.set_options({
                "PE_list": self.pe_list,
                "MP_list": self.mp_cut.option("mp_list"),
                "sample_name": self.option("sample_name"),
                "kmers": self.option("kmers"),
                "min_contig_lgth": self.option("velvet_min_contig_lgth"),
                "min_pair_count": self.option("velvet_min_pair_count")
            })
            self.assemble.on("end", self.run_select)
        elif self.option("assem_tool") == "spades":
            self.assemble = self.add_tool("assemble.assemble_spades")
            self.assemble.set_options({
                "PE_list": self.pe_list,
                "MP_list": self.mp_cut.option("mp_list"),
                "sample_name": self.option("sample_name"),
                "kmer": self.option("kmers")
            })
            self.assemble.on("end", self.run_gaploser)
        self.assemble.run()

    def run_select(self):
        self.scaf_select.set_options({
            "seq_dir": self.assemble.option("scafSeq")
        })
        self.scaf_select.on("end", self.run_gaploser)
        self.scaf_select.run()

    def run_gaploser(self):
        opts = {
            "PE_list": self.pe_list,
            "sample_name": self.option("sample_name")
        }
        if self.option("assem_tool") == "spades":
            opts["seq_scaf"] = self.assemble.option("scf_seq")
        else:
            opts["seq_scaf"] = self.scaf_select.option('scf_seq')
        self.gapcloser.set_options(opts)
        self.gapcloser.on("end", self.set_output)
        self.gapcloser.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        link_dir(self.gapcloser.output_dir, self.output_dir)  # 会不会有别的文件？不会
        # 注意spades的前后结果会不会不一样
        self.logger.info("设置注释结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleHiseqModule, self).end()
