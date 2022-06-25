# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.09.19

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GeneSampleSeqAgent(Agent):
    """
    基因详情页-样本序列
    """
    def __init__(self, parent):
        super(GeneSampleSeqAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},  # ref.fa
            {"name": "ref_gff", "type": "string"},  # gff文件
            {"name": "specimen_ids", "type": "string"},  # 样本
            {"name": "pop_final_vcf", "type": "infile", "format": "sequence.vcf"},  # pop_final_vcf
            {"name": "gene_id", "type": "string"},  # gene_id
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"},
            {"name": "fa_path", "type": "string"},
            {"name": "ssr_path", "type": "string"},   # 基因组路径
            {"name": "gene_name", "type": "string"},
            {"name": "is_wgs_result", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref_gff文件")
        if not self.option("specimen_ids"):
            raise OptionError("请设置specimen_ids")
        if not self.option("pop_final_vcf").prop["path"]:
            raise OptionError("请设置pop_final_vcf文件")
        if not self.option("gene_id"):
            raise OptionError("请设置gene_id")
        if not self.option("main_id"):
            raise OptionError("请设置main_id")
        if not self.option("project_type"):
            raise OptionError("请设置project_type")
        if not self.option("ssr_path"):
            raise OptionError("请设置ssr_path")
        if not self.option("gene_name"):
            raise OptionError("请设置gene_name")
        if not self.option("is_wgs_result"):
            raise OptionError("请设置is_wgs_result")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(GeneSampleSeqAgent, self).end()


class GeneSampleSeqTool(Tool):
    def __init__(self, config):
        super(GeneSampleSeqTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.gene_samples_seq = self.config.PACKAGE_DIR + "/dna_evolution/final.version.pl"
        self.get_exon_seq = self.config.PACKAGE_DIR + "/dna_evolution/01.get-exonseq.pl"
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_geneome"  # 参考组配置文件
        # self.ref_fa = os.path.join(self.json_path, self.option("ref_fa"))
        # self.ref_gff = os.path.join(self.json_path, self.option("ref_gff"))
        self.ref_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome", self.option("ref_fa"))  # 测试
        self.ref_gff = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome", self.option("ref_gff"))  # 测试
        self.ref_protein_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome",
                                           os.path.join(self.option("ssr_path"), "ref.protein.fa"))
        self.out_protein_fa = os.path.join(self.output_dir, "protein.fa")
        self.s3_protein_fa = os.path.join(self.option("fa_path"), "gene_sample_seq/protein.fa")
        self.ref_mRNA_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome",
                                        os.path.join(self.option("ssr_path"), "ref.mRNA.fa"))
        self.out_mRNA_fa = os.path.join(self.output_dir, "mRNA.fa")
        self.s3_mRNA_fa = os.path.join(self.option("fa_path"), "gene_sample_seq/mRNA.fa")
        self.ref_cds_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_geneome",
                                       os.path.join(self.option("ssr_path"), "ref.cds.fa"))
        if self.option("is_wgs_result") == "yes":
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome"
            self.ref_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome", self.option("ref_fa"))
            self.ref_gff = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome", self.option("ref_gff"))
            self.ref_protein_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome",
                                               os.path.join(self.option("ssr_path"), "ref.protein.fa"))
            self.ref_mRNA_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome",
                                            os.path.join(self.option("ssr_path"), "ref.mRNA.fa"))
            self.ref_cds_fa = os.path.join(self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome",
                                           os.path.join(self.option("ssr_path"), "ref.cds.fa"))
        self.out_cds_fa = os.path.join(self.output_dir, "cds.fa")
        self.s3_cds_fa = os.path.join(self.option("fa_path"), "gene_sample_seq/cds.fa")
        self.s3_gene_path = os.path.join(self.option("fa_path"), ("gene_sample_seq/" + self.option("gene_id") + ".fa"))

    def run_gene_samples_seq(self):
        """
        运行final.version.pl脚本，得到样本列表对应的gene_id序列信息。
        """
        with open(self.work_dir + "/sample.list", "w") as w:
            for s in self.option("specimen_ids").split(";"):
                w.write(s + "\n")
        cmd = "{} {} -ref {}".format(self.perl_path, self.gene_samples_seq, self.ref_fa)
        cmd += " -gff {} -vcf {}".format(self.ref_gff, self.option("pop_final_vcf").prop["path"])
        cmd += " -list {} -gene {} -out {}".format(self.work_dir + "/sample.list", self.option("gene_id"),
                                                   self.output_dir + "/" + self.option("gene_id") + ".fa")
        command = self.add_command("get_samples_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gene_samples_seq.pl 运行完成")
        else:
            self.set_error("gene_samples_seq.pl 运行失败")
        seq_api = self.api.api("dna_evolution.gene_sample_seq")
        if self.option("project_type"):
            seq_api._project_type = self.option("project_type")
        gene_path = self.output_dir + "/" + self.option("gene_id") + ".fa"
        input_desc = "运行结束"
        seq_api.add_sg_gene_sample_seq_detail(seq_id=self.option("main_id"), gene_id=self.option("gene_id"),
                                              gene_path=self.s3_gene_path, type="gene", status="end", desc=input_desc)

    def run_get_exon_seq(self):
        """
        运行01.get-exonseq.pl脚本，得到基因结构图输入参数信息。
        """
        cmd = "{} {} -ref {} -gff {} -gene {} -out {}".format(self.perl_path, self.get_exon_seq, self.ref_fa,
                                                              self.ref_gff, self.option("gene_id"), self.output_dir)
        command = self.add_command("get_exon_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gene_samples_seq.pl 运行完成")
        else:
            self.set_error("gene_samples_seq.pl 运行失败")
        seq_api = self.api.api("dna_evolution.gene_sample_seq")
        if self.option("project_type"):
            seq_api._project_type = self.option("project_type")
        exon_path = self.output_dir + "/" + self.option("gene_id") + ".exon.fa"
        with open(exon_path, "r")as fr:
            line = fr.readline()
            if line.startswith(">"):
                seq_api.add_sg_structure(self.option("gene_id"), self.option("main_id"), exon_path)

    def get_fa(self, input_file, out_file, type, s3_file):
        """
        从ref.cds.fa、ref.protein.fa、ref.mRNA.fa中提取序列。
        """
        write_file = 0  # 如果为1就将该行写入文件
        write_times = 0  # 如果为1就是第一行，open写为w
        with open(out_file, "w") as fw:
            fw.write(" ")
        with open(input_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(">"):
                    temp = line.strip().split("gene=")
                    gene_id_list = temp[0].strip().split(">")
                    gene_id = gene_id_list[1]
                    if gene_id == self.option("gene_id"):
                        write_file = 1
                        write_times += 1
                    else:
                        write_file = 0
                    if write_file == 1:
                        if write_times == 1:
                            with open(out_file, "w") as fw:
                                fw.write(line)
                        else:
                            with open(out_file, "a") as fw:
                                fw.write(line)
                else:
                    if write_file == 1:
                        with open(out_file, "a") as fw:
                            fw.write(line)
        seq_api = self.api.api("dna_evolution.gene_sample_seq")
        if self.option("project_type"):
            seq_api._project_type = self.option("project_type")
        if write_times == 0:
            input_desc = type + "文件中没有该gene id：" + self.option("gene_id")
            seq_api.add_sg_gene_sample_seq_detail(seq_id=self.option("main_id"), gene_id=self.option("gene_id"),
                                                  gene_path=s3_file, type=type, status="error", desc=input_desc)
        else:
            input_desc = ""
            seq_api.add_sg_gene_sample_seq_detail(seq_id=self.option("main_id"), gene_id=self.option("gene_id"),
                                                  gene_path=s3_file, type=type, status="end", desc=input_desc)

    def run(self):
        super(GeneSampleSeqTool, self).run()
        self.run_gene_samples_seq()
        self.run_get_exon_seq()
        self.get_fa(self.ref_protein_fa, self.out_protein_fa, "protein", self.s3_protein_fa)
        self.get_fa(self.ref_mRNA_fa, self.out_mRNA_fa, "mRNA", self.s3_mRNA_fa)
        self.get_fa(self.ref_cds_fa, self.out_cds_fa, "cds", self.s3_cds_fa)
        self.end()
