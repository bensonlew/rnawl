# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"
# 20190908

import os
import json
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class GenePredictGffModule(Module):
    """
    细菌比较基因组seq+gff序列预测模块
    分析内容：编码基因预测、tRNA预测、rRNA预测，及其相关信息统计
    这部分的基因预测不再区分染色体和质粒，因为主要是基因预测的软件不同，如果基因预测的结果cds数目为0，前端直接报错，不会传给后台；
    """

    def __init__(self, work_id):
        super(GenePredictGffModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 染色体序列文件
            {"name": "genome_plas", "type": "infile", "format": "sequence.fasta"},  # 质粒序列文件
            {"name": "gene_prefix", "type": "string"},  # 基因前缀，自定义前缀,默认gene
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa", ‘plassmidB’: ‘plasb"}"
            {"name": "genome_name", "type": "string"}, # 基因组名称
            {"name": "genome_type", "type": "string"}, #基因组类型是属于完成图还是扫描图
            {"name": "sample", "type": "string"}, #样品名称
            {"name": "plasmid_prefix", "type": "string"},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa", ‘plassmidB’: ‘plasb"}"
            {"name": "input_gff", "type": "infile", "format": "gene_structure.gff3"},  # 输入的gff文件
            {"name": "trna", "type": "string", "default": "-"},  # tRNA输入数目，如果没有填写就默认为-，如果为0或其他数字不执行预测
            {"name": "rrna", "type": "string", "default": "-"},  # rRNA输入数目，如果没有填写就默认为-，如果为0或其他数字不执行预测
            {"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            {"name": "seq_plas", "type": "outfile", "format": "sequence.fasta"},  # 完成图质粒预测的序列
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(option)
        self.genome_tools = []  # 一个基因组运行编码基因预测、tRNA预测、rRNA预测
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
        if (not self.option("genome").is_set) and (not self.option("genome_plas").is_set):
            raise OptionError("必须提供基因组文件")
        if not self.option("genome_name"):
            raise OptionError("必须提供样品名")
        if not self.option("input_gff").is_set:
            raise OptionError("必须提供所需的gff文件")
        return True

    def run_cds_predict(self):
        """
        进行cds基因预测,从gff文件中提取出所需要的信息
        :return:
        """
        self.cds_predict = self.add_tool("bac_comp_genome.get_seq_from_gff")
        opts = {
            "sample": self.option("sample"),
            "genome": self.option("genome_name"),
            "input_gff": self.option("input_gff"),
            "genome_type": self.option("genome_type"),
            "type": "cds"
        }
        if self.option("gene_prefix"):
            opts["gene_tag"] = self.option("gene_prefix")
        else:
            opts["gene_tag"] = self.option("plasmid_prefix")
        self.cds_predict.set_options(opts)
        self.genome_tools.append(self.cds_predict)

    def run_trna_predict(self):
        """
        进行tRNA功能预测，先对tRNA个数进行判断，然后选择方式方法
        :return:
        """
        if self.option("trna") != "-":
            self.trna_predict = self.add_tool("bac_comp_genome.get_seq_from_gff")
            opts = {
                "sample": self.option("sample"),
                "genome": self.option("genome_name"),
                # "gene_tag": self.option("gene_prefix"),
                "input_gff": self.option("input_gff"),
                "genome_type": self.option("genome_type"),
                "type": "tRNA"
            }
        else:
            self.trna_predict = self.add_tool("bac_comp_genome.trnascanse")
            opts = {
                "input_genome": self.option("genome"),
                "genome_name": self.option("genome_name"),
                "genome_type": self.option("genome_type"),
                # "gene_tag": self.option("gene_prefix")
            }
        if self.option("gene_prefix"):
            opts["gene_tag"] = self.option("gene_prefix")
        else:
            opts["gene_tag"] = self.option("plasmid_prefix")
        self.trna_predict.set_options(opts)
        self.genome_tools.append(self.trna_predict)

    def run_rrna_predict(self):
        """
        进行rRNA功能预测，先对rRNA个数进行判断，然后选择方法去做
        :return:
        """
        if self.option("rrna") != "-":
            self.rrna_predict = self.add_tool("bac_comp_genome.get_seq_from_gff")
            opts = {
                "sample": self.option("sample"),
                "genome": self.option("genome_name"),
                # "gene_tag": self.option("gene_prefix"),
                "input_gff": self.option("input_gff"),
                "genome_type": self.option("genome_type"),
                "type": "rRNA"
            }
        else:
            self.rrna_predict = self.add_tool("bac_comp_genome.barrnap")
            opts = {
                "input_genome": self.option("genome"),
                "genome_name": self.option("genome_name"),
                "genome_type": self.option("genome_type"),
                # "gene_tag": self.option("gene_prefix")
            }
        if self.option("gene_prefix"):
            opts["gene_tag"] = self.option("gene_prefix")
        else:
            opts["gene_tag"] = self.option("plasmid_prefix")
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
            "sample": self.option("sample"),
            "is_gff": "gff+seq"
        }
        self.dna_predict_tidy.set_options(opts)
        self.dna_predict_tidy.on("start", self.set_step, {"start": self.step.predict_tidy})
        self.dna_predict_tidy.on("end", self.set_step, {"end": self.step.predict_tidy})
        self.dna_predict_tidy.on("end", self.set_output, "predict_tidy")
        self.dna_predict_tidy.run()

    def run_dna_predict_tidy_pla(self):
        """
        根据质粒预测结果进行重新矫正，统计基础信息
        :return:
        """
        self.dna_predict_tidyp = self.add_tool("bac_comp_genome.dna_predict_tidy")
        opts = {
            "input_genome": self.option("genome_plas"),
            "plasmid_prefix": json.dumps({self.option("genome_name"): self.option("plasmid_prefix")}),
            "gene_gff": self.cds_predict.option("gff"),
            "trna_gff": self.trna_predict.option("gff"),
            "rrna_gff": self.rrna_predict.option("gff"),
            "genome_name": self.option("genome_name"),
            "genome_type": self.option("genome_type"),
            "is_gff": "gff+seq"
        }
        self.dna_predict_tidyp.set_options(opts)
        self.dna_predict_tidyp.on("start", self.set_step, {"start": self.step.predict_tidy})
        self.dna_predict_tidyp.on("end", self.set_step, {"end": self.step.predict_tidy})
        self.dna_predict_tidyp.on("end", self.set_output, "predict_tidy")
        self.dna_predict_tidyp.run()


    def run(self):
        self.logger.info("开始运行啦")
        super(GenePredictGffModule, self).run()
        if self.option("genome").is_set or self.option("genome_plas").is_set:
            self.logger.info("正在对染色体进行基因预测")
            self.run_cds_predict()
            self.run_trna_predict()
            self.run_rrna_predict()
            if self.option("plasmid_prefix"):
                self.on_rely(self.genome_tools, self.run_dna_predict_tidy_pla)
            else:
                self.on_rely(self.genome_tools, self.run_dna_predict_tidy)
            for tool in self.genome_tools:
                tool.run()

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
            else:
                self.logger.info("没有正确的生成结果文件")
            self.logger.info("设置结果目录文件完成")
            self.end()

    def end(self):
        super(GenePredictGffModule, self).end()
