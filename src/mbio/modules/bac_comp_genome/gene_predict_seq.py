# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"
# 20190908

import os
import json
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class GenePredictSeqModule(Module):
    """
    细菌比较基因组seq序列预测模块
    分析内容：编码基因预测、tRNA预测、rRNA预测，及其相关信息统计
    将上传seq序列分为两种情形：chromosome和plasmid，这样方便进行统计和预测
    """

    def __init__(self, work_id):
        super(GenePredictSeqModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 染色体序列文件
            {"name": "genome_plas", "type": "infile", "format": "sequence.fasta"},  # 完成图时可能有的质粒序列
            {"name": "gene_prefix", "type": "string"},  # 基因前缀，自定义前缀,默认gene
            {"name": "plasmid_prefix", "type": "string"},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa", ‘plassmidB’: ‘plasb"}"
            {"name": "genome_name", "type": "string"}, # 基因组名称
            {"name": "genome_type", "type": "string"}, #基因组类型是属于完成图还是扫描图
            {"name": "sample", "type": "string"}, #样品名称
            {"name": "software_list", "type":"string", "default":"glimmer"},  #使用的软件，逗号分割
            {"name": "p_software_list", "type":"string", "default":"genemark"},  #使用的软件，逗号分割
            {"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            {"name": "seq_plas", "type": "outfile", "format": "sequence.fasta"},  # 完成图质粒预测的序列
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(option)
        self.genome_tools = []  # 一个基因组运行编码基因预测、tRNA预测、rRNA预测
        self.plas_tools = []  # 质粒运行编码基因预测、tRNA预测、rRNA预测
        #self.predict_tidy = [self.dna_predict_tidy]
        self.step.add_steps("trna", "trnap", "cds", "predict_tidy")
        self.tools = []
        self.database = ""

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        #if not (self.option("genome").is_set) and (self.option("genome_plas").is_set):
            #raise OptionError("必须提供基因组文件")
        if not self.option("genome_name"):
            raise OptionError("必须提供样品名")
        return True

    def run_cds_predict(self):
        """
        进行cds基因预测
        :return:
        """
        self.cds_predict = self.add_module("bac_comp_genome.gene_predict_cds")
        opts = {
            "input_genome": self.option("genome"),
            "software_list": self.option("software_list"),
            "trna": self.trna_predict.option("gff"),
            "rrna": self.rrna_predict.option("gff"),
            "orf_prefix": self.option("gene_prefix"),
            "genome_name": self.option("genome_name"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.run_dna_predict_tidy)
        self.cds_predict.run()

    def run_trna_predict(self):
        """
        进行tRNA功能预测
        :return:
        """
        self.trna_predict = self.add_tool("bac_comp_genome.trnascanse")
        opts = {
            "input_genome": self.option("genome"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
            "gene_tag": self.option("gene_prefix"),
        }
        self.trna_predict.set_options(opts)
        self.genome_tools.append(self.trna_predict)
        self.trna_predict.on("start", self.set_step, {"start": self.step.trna})
        self.trna_predict.on("end", self.set_step, {"end": self.step.trna})
        self.trna_predict.on("end", self.set_output, "trna")

    def run_rrna_predict(self):
        """
        进行rRNA功能预测
        :return:
        """
        self.rrna_predict = self.add_tool("bac_comp_genome.barrnap")
        opts = {
            "input_genome": self.option("genome"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
            "gene_tag": self.option("gene_prefix"),
        }
        self.rrna_predict.set_options(opts)
        self.genome_tools.append(self.rrna_predict)

    def run_dna_predict_tidy(self):
        """
        对染色体的预测结果重新矫正，统计基础信息
        :return:
        """
        self.dna_predict_tidy = self.add_tool("bac_comp_genome.dna_predict_tidy")
        opts = {
            "input_genome": self.option("genome"),
            "gene_prefix": self.option("gene_prefix"),
            "gene_gff": self.cds_predict.option("gff"),
            "trna_gff": self.trna_predict.option("gff"),
            "rrna_gff": self.rrna_predict.option("gff"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
        }
        self.dna_predict_tidy.set_options(opts)
        self.dna_predict_tidy.on("start", self.set_step, {"start": self.step.predict_tidy})
        self.dna_predict_tidy.on("end", self.set_step, {"end": self.step.predict_tidy})
        self.dna_predict_tidy.on("end", self.set_output, "predict_tidy")
        self.dna_predict_tidy.run()

    def run_cds_predict_plas(self):
        """
        进行cds基因预测genemark
        :return:
        """
        opts = {
            "input_genome": self.option("genome_plas"),
            "genome_type": "plas",
            "orf_prefix": self.option("plasmid_prefix"),
            "genome_name": self.option("genome_name"),
            'software' : self.option("p_software_list"),
        }
        self.cds_predict_p = self.add_tool("bac_comp_genome.glimmer_and_genemark")
        self.cds_predict_p.set_options(opts)
        self.cds_predict_p.on("end", self.run_dna_predict_tidy_pla)
        self.cds_predict_p.run()

    def run_trna_predict_plas(self):
        """
        进行tRNA预测
        :return:
        """
        opts = {
            "input_genome": self.option("genome_plas"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
            "gene_tag": self.option("plasmid_prefix"),
        }
        self.trna_predictp = self.add_tool("bac_comp_genome.trnascanse")
        self.trna_predictp.set_options(opts)
        self.plas_tools.append(self.trna_predictp)
        self.trna_predictp.on("start", self.set_step, {"start": self.step.trnap})
        self.trna_predictp.on("end", self.set_step, {"end": self.step.trnap})
        self.trna_predictp.on("end", self.set_output, "trnap")

    def run_rrna_predict_plas(self):
        """
        进行rRNA预测
        :return:
        """
        opts = {
            "input_genome": self.option("genome_plas"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
            "gene_tag": self.option("plasmid_prefix"),
        }
        self.rrna_predictp = self.add_tool("bac_comp_genome.barrnap")
        self.rrna_predictp.set_options(opts)
        self.plas_tools.append(self.rrna_predictp)

    def run_dna_predict_tidy_pla(self):
        """
        根据质粒预测结果进行重新矫正，统计基础信息
        :return:
        """
        self.dna_predict_tidyp = self.add_tool("bac_comp_genome.dna_predict_tidy")
        opts = {
            "input_genome": self.option("genome_plas"),
            "plasmid_prefix": json.dumps({self.option("genome_name"): self.option("plasmid_prefix")}),
            "gene_gff": self.cds_predict_p.option("gff"),
            "trna_gff": self.trna_predictp.option("gff"),
            "rrna_gff": self.rrna_predictp.option("gff"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
        }
        self.dna_predict_tidyp.set_options(opts)
        self.dna_predict_tidyp.on("start", self.set_step, {"start": self.step.predict_tidy})
        self.dna_predict_tidyp.on("end", self.set_step, {"end": self.step.predict_tidy})
        self.dna_predict_tidyp.on("end", self.set_output, "predict_tidy")
        self.dna_predict_tidyp.run()

    def run(self):
        super(GenePredictSeqModule, self).run()
        if self.option("genome").is_set:
            self.logger.info("正在对染色体进行基因预测")
            self.run_trna_predict()
            self.run_rrna_predict()
            self.on_rely(self.genome_tools, self.run_cds_predict)
            for tool in self.genome_tools:
                tool.run()

        elif self.option("genome_plas").is_set:
            self.logger.info("正在对质粒进行基因预测")
            self.run_trna_predict_plas()
            self.run_rrna_predict_plas()
            self.on_rely(self.plas_tools, self.run_cds_predict_plas)
            for tool_p in self.plas_tools:
                tool_p.run()

    def set_output(self, event):
        """
        设置结果目录
        :param event:
        :return:
        """
        self.logger.info("正在设置结果目录")
        if event["data"] == "predict_tidy":
            if self.option("genome").is_set:
                dir_path = self.dna_predict_tidy.output_dir
                for file in os.listdir(dir_path):
                    file_path = os.path.join(self.dna_predict_tidy.output_dir, file)
                    new_file = os.path.join(self.output_dir, file)
                    if os.path.exists(new_file):
                        os.remove(new_file)
                    os.link(file_path, new_file)
            elif self.option("genome_plas").is_set:
                dir_path = self.dna_predict_tidyp.output_dir
                for file in os.listdir(dir_path):
                    file_path = os.path.join(self.dna_predict_tidyp.output_dir, file)
                    new_file = os.path.join(self.output_dir, file)
                    if os.path.exists(new_file):
                        os.remove(new_file)
                    os.link(file_path, new_file)
            self.end()

    def end(self):
        super(GenePredictSeqModule, self).end()
