bam_addrg.py#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
from mbio.packages.ref_rna_v3.snp_anno import snp_anno
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re
import unittest


class AnnovarAgent(Agent):
    """
    Annovar:用对处理vcf格式文件/注释突变信息
    version 1.0
    author: qindanhua
    last_modify: 2016.12.30
    """

    def __init__(self, parent):
        super(AnnovarAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组类型
            {"name": "input_file", "type": "infile", "format": "gene_structure.vcf,gene_structure.vcf_dir"},  # 输入文件
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件
            {"name": "combine_vcf", "type": "bool", "default": False},  # 输入文件
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('annovar')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.annovar.start()
        self.step.update()

    def step_end(self):
        self.step.annovar.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("input_file").is_set:
            raise OptionError("请输入VCF格式文件", code="35600802")
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 11
        self._memory = '30G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./snp_anno.xls", "xls", "snp注释结果表"]
        ])
        super(AnnovarAgent, self).end()


class AnnovarTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(AnnovarTool, self).__init__(config)
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK.jar"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        self.perl_path = "miniconda2/bin/"
        self.perl_full_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.annovar_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/annovar/"
        self.gtfToGenePred_path = "/bioinfo/gene-structure/annovar/"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.ref_fasta = ''
        self.ref_gtf = ''
        self.vcf_path = ''
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def dict(self):
        """
        使用picard对参考基因组构建字典
        """
        if 'path' in self.option("ref_fasta").prop.keys():
            ref_fasta = self.option("ref_fasta").prop["path"]
        else:
            ref_fasta = self.option("ref_fasta")
        dict_name = os.path.dirname(ref_fasta) + "/" + ".".join(os.path.basename(ref_fasta).split(".")[:-1]) + ".dict"
        cmd = "program/sun_jdk1.8.0/bin/java -jar {}picard.jar CreateSequenceDictionary R={} O={}" \
            .format(self.picard_path, ref_fasta, dict_name)
        if os.path.exists(dict_name):
            os.remove(dict_name)
        self.logger.info("开始用picard对参考基因组构建字典")
        self.logger.info("cmd")
        command = self.add_command("dict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("参考基因组构建dict done!")
        else:
            self.set_error("构建dict过程error！", code="35600807")

    def combine_vcf(self):
        samples_option = ''
        # 因为samtools产生的就是一个vcf,在module的参数,不能作为glob的输入
        if os.path.isdir(self.option("input_file").prop["path"]):
            vcf_files = glob.glob('{}/*.vcf'.format(self.option("input_file").prop["path"]))
            for vf in sorted(vcf_files):
                samples_option += ' --variant {}'.format(vf)
        else:
            vcf_files = self.option("input_file").prop["path"]
            samples_option += ' --variant {}'.format(vcf_files)

        cmd = '{}java -jar {} -R {} -T CombineVariants {} -o Combine_Variants.vcf -genotypeMergeOptions UNIQUIFY' \
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop["path"], samples_option)
        command = self.add_command("combine_variants", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行combine_variants结束")
        else:
            self.set_error("运行combine_variants出错", code="35600808")

    def get_genome(self):
        if self.option("ref_genome") == "customer_mode":
            self.ref_fasta = self.option("ref_fasta").prop["path"]
            self.ref_gtf = self.option("ref_gtf").prop["path"]
        else:
            ref_genome_json = self.config.SOFTWARE_DIR + "/database/refGenome/scripts/ref_genome.json"
            with open(ref_genome_json, "r") as f:
                ref_dict = json.loads(f.read())
                self.ref_fasta = ref_dict[self.option("ref_genome")]["ref_genome"]
                self.ref_gtf = ref_dict[self.option("ref_genome")]["gtf"]

    def gtf_to_genepred(self):
        ret = self.check_gtf(self.ref_gtf)
        if ret:
            cmd = "{}gtfToGenePred -genePredExt {} {}.genes_refGene.tmp.txt".format(
                self.gtfToGenePred_path, self.ref_gtf, self.option("ref_genome")
            )
        else:
            cmd = "{}gtfToGenePred {} {}.genes_refGene.tmp.txt".format(
                self.gtfToGenePred_path, self.ref_gtf, self.option("ref_genome")
            )
            # 2019.04.12 bug 虽然可以生成genepred，但是转换成的refgene并不能被annover接受
            self.set_error('find format error in %s, abord', variables=(self.ref_gtf), code="35600809")
        command = self.add_command("gtftogenepred", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gtfToGenePred完成!")
            self.logger.info("awk \'{print 1\"\\t\"$0}\'" + " {}.genes_refGene.tmp.txt  > {}_refGene.tmp.txt"
                             .format(self.option("ref_genome"), self.option("ref_genome")))
            self.logger.info(
                "awk -F \"\t\" \'{if ($4 == \"+\" || $4 == \"-\") {print $0}}\'" + " {}_refGene.tmp.txt > {}_refGene.txt"
                .format(self.option("ref_genome"), self.option("ref_genome")))
            try:
                subprocess.check_output("awk \'{print 1\"\\t\"$0}\'" + " {}.genes_refGene.tmp.txt  > {}_refGene.tmp.txt"
                                        .format(self.option("ref_genome"), self.option("ref_genome")), shell=True)
                subprocess.check_output(
                    "awk -F \"\t\" \'{if ($4 == \"+\" || $4 == \"-\") {print $0}}\'" + " {}_refGene.tmp.txt > {}_refGene.txt"
                    .format(self.option("ref_genome"), self.option("ref_genome")), shell=True)
                self.logger.info("提取gtfToGenePred结果信息完成")
                return True
            except subprocess.CalledProcessError:
                self.logger.info("提取gtfToGenePred结果信息出错")
                return False
        else:
            self.set_error("运行gtfToGenePred出错！", code="35600810")

    def retrieve_seq_from_fasta(self):
        cmd = "{}perl {}retrieve_seq_from_fasta.pl -format refGene -seqfile {} --outfile {}_refGeneMrna.fa " \
              "{}_refGene.txt".format(self.perl_path, self.annovar_path, self.ref_fasta,
                                      self.option("ref_genome"), self.option("ref_genome"))
        command = self.add_command("retrieve_seq_from_fasta", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0 or None:
            self.logger.info("运行retrieve_seq_from_fasta完成!")
            if os.path.exists("./geneomedb"):
                shutil.rmtree("./geneomedb")
            os.system("mkdir geneomedb")
            os.system("mv {}.genes_refGene.tmp.txt {}_refGene.txt {}_refGeneMrna.fa ./geneomedb/"
                      .format(self.option("ref_genome"), self.option("ref_genome"), self.option("ref_genome")))
        else:
            self.set_error("运行retrieve_seq_from_fasta出错！", code="35600811")

    def convert2annovar(self):
        cmd = "{}perl {}convert2annovar.pl -format vcf4 {} > clean.snp.avinput" \
            .format(self.perl_full_path, self.annovar_path, self.vcf_path)
        cmd = "{}perl {}convert2annovar.pl -format vcf4old {} > clean.snp.avinput" \
            .format(self.perl_full_path, self.annovar_path, self.vcf_path)
        self.logger.info(self.vcf_path)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("提取convert2annovar结果信息完成")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("提取convert2annovar结果信息出错")
            return False

    def annotate_variation(self):
        cmd = "{}perl {}annotate_variation.pl -buildver {} -dbtype refGene clean.snp.avinput ./geneomedb" \
            .format(self.perl_path, self.annovar_path, self.option("ref_genome"), self.option("input_file"))
        command = self.add_command("annotate_variation", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行annotate_variation完成!")
        else:
            self.set_error("运行annotate_variation出错！", code="35600812")

    def set_output(self):
        self.logger.info("snp annotation")
        snp_anno("./clean.snp.avinput.variant_function", "./clean.snp.avinput.exonic_variant_function", self.vcf_path,
                 "./snp_anno.xls")
        self.logger.info("snp annotation done")
        self.logger.info("set ouptput")
        self.logger.info(self.output_dir + "/snp_anno.xls")
        if os.path.exists(self.output_dir + "/snp_anno.xls"):
            os.remove(self.output_dir + "/snp_anno.xls")
        os.link(self.work_dir + "/snp_anno.xls", self.output_dir + "/snp_anno.xls")
        # os.link(self.work_dir + "/snp_annotation.xls", self.output_dir + "/snp_annotation.xls")
        self.logger.info("set ouptput done")
        # if not self.option('combine_vcf'):
        #     vcf_f = open(self.vcf_path, "r")
        #     vcf_used = vcf_f.read().split("#CHROM")[1]
        #     genetype_type_ori = open(self.work_dir + "/genotype_ori.xls", "w")
        #     genetype_type_ori.write("#CHROM" + vcf_used)
        #     genetype_type = open(self.output_dir + "/genotype.xls", "w")
        #     genetype_type_ori.close()
        #     genetype_infos = open(self.work_dir + "/genotype_ori.xls", "r")
        #     genetype_info = genetype_infos.readline()
        #     genetype_type.write("#CHROM" + "\t" + "Start\t" + "Ref\t" + "Alt")
        #     for i in genetype_info.strip().split("\t")[9:]:
        #         i_info = "\t" + i + "genotype"
        #         genetype_type.write(i_info)
        #     genetype_type.write("\n")
        #     for line in genetype_infos.readlines():
        #         line_infos = line.strip().split("\t")
        #         if len(line_infos[4].split(",")) > 1:
        #             for altn in range(0,len(line_infos[4].split(","))):
        #                 genetype_type.write(line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[altn])
        #                 for detail in line_infos[9:]:
        #                     detail = detail.split(":")
        #                     if detail[0] == "./.":
        #                         genetype_type.write("\t" + "./.")
        #                     elif detail[0] == "0/0":
        #                         genetype_type.write("\t" + line_infos[3] + line_infos[3])
        #                     # detail_depth = "."
        #                     else:
        #                         gts = detail[0].split("/")
        #                         if gts[0] == "0":
        #                             if int(gts[0])-1> altn:
        #                                 genetype_type.write("\t" + "./.")
        #                             else:
        #                                 if int(gts[1]) - 1 >= altn:
        #                                     genetype_type.write("\t" + line_infos[3])
        #                                     genetype_type.write(line_infos[4].split(",")[altn])
        #                                 else:
        #                                     genetype_type.write("\t" + "./.")
        #                         else:
        #                             if int(gts[0]) - 1 > altn:
        #                                 genetype_type.write("\t" + "./.")
        #                             else:
        #                                 if int(gts[1]) - 1 >= altn:
        #                                     genetype_type.write("\t" + line_infos[4].split(",")[altn])
        #                                     genetype_type.write(line_infos[4].split(",")[altn])
        #                                 else:
        #                                     genetype_type.write("\t" +"./.")
        #                     # else:
        #                     # detail_depth = detail[1].split(",")[0] + "," + detail[1].split(",")[1]
        #                 genetype_type.write("\n")
        #             # genetype_type.write(
        #             #     line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[
        #             #         0])
        #             # for detail in line_infos[9:]:
        #             #     detail = detail.split(":")
        #             #     # detail_depth = detail[3].split(",")[0] + "," + detail[3].split(",")[1]
        #             #     genetype_type.write("\t" + detail[0])
        #             # genetype_type.write("\n")
        #             # genetype_type.write(
        #             #     line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[
        #             #         1])
        #             # for detail in line_infos[9:]:
        #             #     detail = detail.split(":")
        #             #     # detail_depth = detail[3].split(",")[0] + "," + detail[3].split(",")[1]
        #             #     genetype_type.write("\t" + detail[0])
        #             # genetype_type.write("\n")
        #         else:
        #             genetype_type.write(line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[0])
        #             for detail in line_infos[9:]:
        #                 detail = detail.split(":")
        #                 # detail_depth = detail[3]
        #                 if detail[0] == "./." :
        #                     genetype_type.write("\t" + "./.")
        #                 elif detail[0] == "0/0":
        #                     genetype_type.write("\t" + line_infos[3] + line_infos[3])
        #                 elif detail[0] =="0/1":
        #                     genetype_type.write("\t" + line_infos[3] + line_infos[4])
        #                 elif detail[0] == "1/1":
        #                     genetype_type.write("\t" + line_infos[4] + line_infos[4])
        #             genetype_type.write("\n")
        # else:
        #     vcf_f = open(self.vcf_path, "r")
        #     vcf_used = vcf_f.read().split("#CHROM")[1]
        #     genetype_type_ori = open(self.work_dir + "/genotype_ori.xls", "w")
        #     genetype_type_ori.write("#CHROM" + vcf_used)
        #     genetype_type = open(self.output_dir + "/genotype.xls", "w")
        #     genetype_type_ori.close()
        #     genetype_infos = open(self.work_dir + "/genotype_ori.xls", "r")
        #     genetype_info = genetype_infos.readline()
        #     genetype_type.write("#CHROM" + "\t" + "Start\t" + "Ref\t" + "Alt")
        #     for i in genetype_info.strip().split("\t")[9:]:
        #         i_info = "\t" + i.split(".")[0] + "_genotype"
        #         genetype_type.write(i_info)
        #     genetype_type.write("\n")
        #     for line in genetype_infos.readlines():
        #         line_infos = line.strip().split("\t")
        #         if len(line_infos[4].split(",")) > 1:
        #             for altn in range(0,len(line_infos[4].split(","))):
        #                 genetype_type.write(line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[altn])
        #                 for detail in line_infos[9:]:
        #                     detail = detail.split(":")
        #                     if detail[0] == "./.":
        #                         genetype_type.write("\t" +"./.")
        #                     elif detail[0] == "0/0":
        #                         genetype_type.write("\t" + line_infos[3] + line_infos[3])
        #                     # detail_depth = "."
        #                     else:
        #                         gts=detail[0].split("/")
        #                         if gts[0] == "0":
        #                             if int(gts[0])-1> altn:
        #                                 genetype_type.write("\t" + "./.")
        #                             else:
        #                                 if int(gts[1])-1 >= altn:
        #                                     genetype_type.write("\t" + line_infos[3])
        #                                     genetype_type.write(line_infos[4].split(",")[altn])
        #                                 else:
        #                                     genetype_type.write("\t" + "./.")
        #                         else:
        #                             if int(gts[0]) - 1 > altn:
        #                                 genetype_type.write("\t" + "./.")
        #                             else:
        #                                 if int(gts[1]) - 1 >= altn:
        #                                     genetype_type.write("\t" + line_infos[4].split(",")[altn])
        #                                     genetype_type.write(line_infos[4].split(",")[altn])
        #                                 else:
        #                                     genetype_type.write("\t" + "./.")
        #                     #genetype_type.write("\t" + detail[0])
        #                     # else:
        #                     # detail_depth = detail[1].split(",")[0] + "," + detail[1].split(",")[1]
        #                 genetype_type.write("\n")
        #                 # genetype_type.write(line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[
        #                 #         1])
        #                 # for detail in line_infos[9:]:
        #                 #     detail = detail.split(":")
        #                 #     # if detail[0] =="./.":
        #                 #     # detail_depth = "."
        #                 #     genetype_type.write("\t" + detail[0])
        #                 #     # else:
        #                 #     # detail_depth = detail[1].split(",")[0] + "," + detail[1].split(",")[2]
        #                 # genetype_type.write("\n")
        #         else:
        #             genetype_type.write(
        #                 line_infos[0] + "\t" + line_infos[1] + "\t" + line_infos[3] + "\t" + line_infos[4].split(",")[
        #                     0])
        #             for detail in line_infos[9:]:
        #                 detail = detail.split(":")
        #                 # if detail[0] =="./.":
        #                 # detail_depth = "."
        #                 if detail[0] == "./." :
        #                     genetype_type.write("\t" +"./.")
        #                 elif detail[0] == "0/0":
        #                     genetype_type.write("\t" + line_infos[3] + line_infos[3])
        #                 elif detail[0] =="0/1":
        #                     genetype_type.write("\t" + line_infos[3] + line_infos[4])
        #                 elif detail[0] == "1/1":
        #                     genetype_type.write("\t" + line_infos[4] + line_infos[4])
        #                 # else:
        #                 #  detail_depth = detail[1]
        #                 # genetype_type.write("\t" + detail_depth)
        #             genetype_type.write("\n")
        # genetype_type.close()
        # df1 = pd.read_table(self.output_dir + "/snp_anno.xls",error_bad_lines=False)
        # df2 = pd.read_table(self.output_dir + "/genotype.xls",error_bad_lines=False)
        # df22 = df2.iloc[:, 4:]
        # df3 = pd.concat([df1, df22], axis=1)
        # df3.to_csv(self.output_dir + "/snp_anno.xls", sep="\t", index=False)
        self.logger.info("set output done")

    def run(self):
        super(AnnovarTool, self).run()
        self.vcf_path = self.option("input_file").prop["path"]
        if self.option('combine_vcf'):
            self.dict()
            self.combine_vcf()
            self.vcf_path = "Combine_Variants.vcf"
        self.get_genome()
        self.gtf_to_genepred()
        self.retrieve_seq_from_fasta()
        self.convert2annovar()
        self.annotate_variation()
        self.set_output()
        self.desc(self.option("des"), self.option("des_type"))
        self.add_name_des(self.work_dir + "/snp_anno.xls", self.work_dir + "/id_name_desc.txt")
        self.end()

    def desc(self, des, des_type):
        """
        获取蛋白名或描述信息
        """
        with open(des, "rb") as f, open(self.work_dir + "/id_name_desc.txt", "w") as w:
            seq1 = ("gene_id", "symbol", "des")
            w.write("\t".join(seq1) + "\n")
            for line in f:
                line = line.strip().split('\t')
                gene_id = line[0]
                if des_type == "type3":
                    symbol = "-"
                    des = line[3] if line[3] else "-"
                    seq = (gene_id, symbol, des)
                    w.write("\t".join(seq) + "\n")
                elif des_type == "type2":
                    symbol = line[2] if line[2] else "-"
                    des = line[5] if line[5] else "-"
                    seq = (gene_id, symbol, des)
                    w.write("\t".join(seq) + "\n")
                elif des_type == "type1":
                    symbol = line[2] if line[2] else "-"
                    des = line[7] if line[7] else "-"
                    seq = (gene_id, symbol, des)
                    w.write("\t".join(seq) + "\n")
                else:
                    pass
        df1 = pd.read_table(self.work_dir + "/id_name_desc.txt", header=0, sep="\t")
        df1.drop_duplicates(keep='first', inplace=True)
        df1 = df1.fillna('-')
        os.remove(self.work_dir + "/id_name_desc.txt")
        df1.to_csv(self.work_dir + "/id_name_desc.txt", index=False, sep="\t", header=True)

    def add_name_des(self, snp_anno_new, id_name_desc):
        with open(snp_anno_new) as f1, open(id_name_desc) as f2, open(self.work_dir + "/snp_annotation.xls", "w") as w1:
            dict_id_symbol = collections.defaultdict(dict)
            dict_id_desc = collections.defaultdict(dict)
            _ = f2.readline()
            for line in f2:
                gene_id, symbol, des = line.strip().split("\t")  # 此处一定要strip, 要不然没个des后面都有一个换行符
                dict_id_symbol[gene_id] = symbol
                dict_id_desc[gene_id] = des
            f1_col_names = f1.readline().strip().split("\t")
            w1.write("\t".join(f1_col_names[0:8]) + "\t" + "Gene name" + "\t" + "Gene description" + "\t" + "\t".join(
                f1_col_names[8:]) + "\n")
            count = 0
            for line in f1:
                count += 1
                whole = line.strip().split("\t")
                tmp_list1 = whole[6].split(";")
                geneby = whole[7]
                gene_list = list()
                symbol_list = list()
                desc_list = list()
                if len(tmp_list1) >= 2:
                    print("现在大于或者等于2")
                    if "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是1")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 and (
                            "ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是2")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and (
                            "intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是3")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and (
                            "ncRNA_splicing" not in tmp_list1) and "intergenic" and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是4")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif (
                            "splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是5")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and (
                            "intergenic" not in tmp_list1) and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是6")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and (
                            "ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and (
                            "ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是7")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and ("ncRNA_splicing" not in tmp_list1) and (
                            "intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是8")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1(
                            "ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是9")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and (
                            "intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是10")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            elif tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and (
                            "ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是11")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and ("ncRNA_splicing" not in tmp_list1) and (
                            "intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是12")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and (
                            "ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and (
                            "ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是13")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "intergenic" or "ncRNA_intergenic":
                                tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                                tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                                gene_list.append(tmp_intergenic_list)
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and (
                            "intergenic" not in tmp_list1) and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是14")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and (not "ncRNA_splicing" in tmp_list1) and (
                    not "intergenic" in tmp_list1) and (not "ncRNA_intergenic" in tmp_list1):
                        print("现在的情况是15")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in
                                                             geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
                                    tmp_splicing_list = [x[0] for x in
                                                         tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                                    gene_list.append(tmp_splicing_list)
                                elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                    tmp_other_list = geneby_splicing.split(",")
                                    gene_list.append(tmp_other_list)
                                elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                    for m in re.split(r"[(,]", geneby_splicing):
                                        if ":" not in m:
                                            gene_list.append([m])
                                else:
                                    pass
                            else:
                                tmp_other_list = geneby.split(",")[i]
                                gene_list.append(tmp_other_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")

                    else:
                        print("现在的情况是16")
                        tmp_other_list = re.split(r"[;,]", geneby)
                        gene_list.append(tmp_other_list)
                        print(gene_list)
                        for i in gene_list:
                            for j in i:
                                if j in dict_id_symbol.keys():
                                    tmp_symbol = dict_id_symbol[j]
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = dict_id_desc[j]
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                                else:
                                    tmp_symbol = "-"
                                    symbol_list.append(tmp_symbol)
                                    symbol_list_str = ";;".join(symbol_list)
                                    tmp_desc = "-"
                                    desc_list.append(tmp_desc)
                                    desc_list_str = ";;".join(desc_list)
                        w1.write(
                            "\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                                whole[8:]) + "\n")
                else:
                    print("现在tem_list是1")
                    if tmp_list1[0] == "splicing" or tmp_list1[0] == "ncRNA_splicing":
                        if ")," in geneby:
                            tmp_splicing_list_pre = [x.split("(") for x in geneby.split("),")]  # 获取对应位置的splicing的信息
                            tmp_splicing_list = [x[0] for x in
                                                 tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
                            gene_list.append(tmp_splicing_list)
                        elif ")," not in geneby and ")" not in geneby:
                            tmp_other_list = geneby.split(",")
                            gene_list.append(tmp_other_list)
                        elif ")," not in geneby and ")" in geneby:
                            for m in re.split(r"[(,]", geneby):
                                if ":" not in m:
                                    gene_list.append([m])
                        else:
                            pass

                    elif tmp_list1[0] == "intergenic" or tmp_list1[0] == "ncRNA_intergenic":
                        tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(",")]
                        tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                        gene_list.append(tmp_intergenic_list)
                    else:
                        tmp_other_list = geneby.split(",")
                        gene_list.append(tmp_other_list)
                    symbol_list_str = "-"
                    desc_list_str = "-"
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    if len(gene_list) == 0:
                        print
                        "gene_list is empty\t".join(whole[0:8])

                    w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(
                        whole[8:]) + "\n")

    def check_gtf(self, gtf):
        ret = True
        for line in open(gtf):
            if line[0] != '#':
                items = line.strip().split('\t')
                if len(items) == 9:
                    if items[2].lower() == 'cds' and items[7] not in ('0', '1', '2'):
                        ret = False
                        break
        return ret


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.annovar",
            "instant": False,
            "options": dict(
                #ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                #bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/bamlist_new",
                #call_type="sentieon",
                #ref_fasta=test_dir+"/"+"Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                #output=None,
                ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna/output_vcf",
                ref_genome="customer_mode",
                combine_vcf=False,
                des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt",
                des_type="type1",
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/biomart/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.biomart",
                #des_type="type3",
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/gtf/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.gtf",
                #ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/dna/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0_genomic.fna",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190716/Snp_tsg_34853_2442_3325/CallSnpIndel/VcfFilterGatk/output/final.pop.indel.filter.recode.vcf",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190617/Snp_tsg_33912_1820_7969/SnpRna/output_vcf"
                input_file="/mnt/ilustre/users/sanger-dev/workspace/20190716/Snp_tsg_34853_2442_3325/CallSnpIndel/VcfFilterGatk/output/final.vcf"
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190621/Snp_tsg_34423_8856_3754/CallSnpIndel/GvcfTypingV2/output/pop.variant.vcf"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()