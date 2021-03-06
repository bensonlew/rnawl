#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
from mbio.packages.prok_rna.snp_anno import snp_anno
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re


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
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 输入文件
            {"name": "combine_vcf", "type": "bool", "default": False},  # 输入文件
            {"name": "id2name", "type": "string"}
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
                raise OptionError("请输入VCF格式文件", code = "35002201")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '5G'

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

    def dict(self):
        """
        使用picard对参考基因组构建字典
        """
        if 'path' in self.option("ref_fasta").prop.keys():
            ref_fasta = self.option("ref_fasta").prop["path"]
        else:
            ref_fasta = self.option("ref_fasta")
        dict_name = os.path.dirname(ref_fasta) + "/" + ".".join(os.path.basename(ref_fasta).split(".")[:-1]) + ".dict"
        cmd = "program/sun_jdk1.8.0/bin/java -jar {}picard.jar CreateSequenceDictionary R={} O={}"\
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
            self.set_error("构建dict过程error", code = "35002202")

    def combine_vcf(self):
        samples_option = ''
        # 因为samtools产生的就是一个vcf,在module的参数,不能作为glob的输入
        if  os.path.isdir(self.option("input_file").prop["path"]):
            vcf_files = glob.glob('{}/*.vcf'.format(self.option("input_file").prop["path"]))
            for vf in sorted(vcf_files):
                samples_option += ' --variant {}'.format(vf)
        else:
            vcf_files = self.option("input_file").prop["path"]
            samples_option += ' --variant {}'.format(vcf_files)

        cmd = '{}java -jar {} -R {} -T CombineVariants {} -o Combine_Variants.vcf -genotypeMergeOptions UNIQUIFY'\
            .format(self.java_path, self.gatk_path, self.option("ref_fasta").prop["path"], samples_option)
        command = self.add_command("combine_variants", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行combine_variants结束")
        else:
            self.set_error("运行combine_variants出错", code = "35002203")

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
        cmd = "{}gtfToGenePred -genePredExt {} {}.genes_refGene.tmp.txt"\
            .format(self.gtfToGenePred_path, self.ref_gtf, self.option("ref_genome"))
        command = self.add_command("gtftogenepred", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gtfToGenePred完成!")
            self.logger.info("awk \'{print 1\"\\t\"$0}\'" + " {}.genes_refGene.tmp.txt  > {}_refGene.tmp.txt"
                             .format(self.option("ref_genome"), self.option("ref_genome")))
            self.logger.info("awk -F \"\t\" \'{if ($4 == \"+\" || $4 == \"-\") {print $0}}\'" + " {}_refGene.tmp.txt > {}_refGene.txt"
                             .format(self.option("ref_genome"), self.option("ref_genome")))
            try:
                subprocess.check_output("awk \'{print 1\"\\t\"$0}\'" + " {}.genes_refGene.tmp.txt  > {}_refGene.tmp.txt"
                                        .format(self.option("ref_genome"), self.option("ref_genome")), shell=True)
                subprocess.check_output("awk -F \"\t\" \'{if ($4 == \"+\" || $4 == \"-\") {print $0}}\'" + " {}_refGene.tmp.txt > {}_refGene.txt"
                                        .format(self.option("ref_genome"), self.option("ref_genome")), shell=True)
                self.logger.info("提取gtfToGenePred结果信息完成")
                return True
            except subprocess.CalledProcessError:
                self.logger.info("提取gtfToGenePred结果信息出错")
                return False
        else:
            self.set_error("运行gtfToGenePred出错", code = "35002204")

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
            self.set_error("运行retrieve_seq_from_fasta出错", code = "35002205")

    def convert2annovar(self):
        cmd = "{}perl {}convert2annovar.pl -format vcf4 {} > clean.snp.avinput"\
            .format(self.perl_full_path, self.annovar_path, self.vcf_path)
        cmd = "{}perl {}convert2annovar.pl -format vcf4old {} > clean.snp.avinput"\
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
        cmd = "{}perl {}annotate_variation.pl -buildver {} -dbtype refGene clean.snp.avinput ./geneomedb"\
            .format(self.perl_path, self.annovar_path, self.option("ref_genome"), self.option("input_file"))
        command = self.add_command("annotate_variation", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行annotate_variation完成!")
        else:
            self.set_error("运行annotate_variation出错", code = "35002206")

    def set_output(self):
        self.logger.info("snp annotation")
        snp_anno("./clean.snp.avinput.variant_function", "./clean.snp.avinput.exonic_variant_function", self.vcf_path, "./snp_anno.xls")
        self.logger.info("snp annotation done")
        self.logger.info("set ouptput")
        self.logger.info(self.output_dir + "/snp_anno.xls")
        if os.path.exists(self.output_dir + "/snp_anno.xls"):
            os.remove(self.output_dir + "/snp_anno.xls")
        os.link(self.work_dir + "/snp_anno.xls", self.output_dir + "/snp_anno.xls")
        # os.link(self.work_dir + "/snp_annotation.xls", self.output_dir + "/snp_annotation.xls")
        self.logger.info("set ouptput done")

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
        self.desc(self.option("id2name"))
        self.add_name_des(self.work_dir + "/snp_anno.xls", self.work_dir + "/id_name_desc.txt")
        self.end()

    def desc(self, id2name):
        """
        获取蛋白名或描述信息
        """
        df1 = pd.read_table(id2name, header=None, sep="\t")
        df1.iloc[:, [6,5,9]].to_csv(self.work_dir + "/id_name_desc.txt", index=False, sep="\t", header=False)

    def add_name_des(self, snp_anno_new, id_name_desc):
        with open(snp_anno_new) as f1, open(id_name_desc) as f2, open(self.work_dir + "/snp_annotation.xls", "w") as w1:
            dict_id_symbol = collections.defaultdict(dict)
            dict_id_desc = collections.defaultdict(dict)
            # _ = f2.readline()
            for line in f2:
                try:
                    gene_id, symbol, des = line.strip().split("\t")  # 此处一定要strip, 要不然没个des后面都有一个换行符
                    dict_id_symbol[gene_id] = symbol
                    dict_id_desc[gene_id] = des
                except:
                    gene_id, symbol = line.strip().split('\t')
                    dict_id_symbol[gene_id] = symbol
                    dict_id_desc[gene_id] = 'no description'
            f1_col_names = f1.readline().strip().split("\t")
            w1.write("\t".join(f1_col_names[0:8]) + "\t" + "Gene name" + "\t" + "Gene description" + "\t" + "\t".join(f1_col_names[8:]) + "\n")
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
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是2")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and ("intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是3")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and ("ncRNA_splicing" not in tmp_list1) and "intergenic" and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是4")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是5")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and "ncRNA_splicing" in tmp_list1 and ("intergenic" not in tmp_list1) and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是6")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and ("ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是7")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and ("ncRNA_splicing" not in tmp_list1) and ("intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是8")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and "intergenic" in tmp_list1 ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是9")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and ("intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
                        print("现在的情况是10")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and ("ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and "ncRNA_intergenic" in tmp_list1:
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and ("ncRNA_splicing" not in tmp_list1) and ("intergenic" not in tmp_list1) and "ncRNA_intergenic" in tmp_list1:
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and ("ncRNA_splicing" not in tmp_list1) and "intergenic" in tmp_list1 and ("ncRNA_intergenic" not in tmp_list1):
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif ("splicing" not in tmp_list1) and "ncRNA_splicing" in tmp_list1 and ("intergenic" not in tmp_list1) and ("ncRNA_intergenic" not in tmp_list1):
                        print("现在的情况是14")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

                    elif "splicing" in tmp_list1 and (not "ncRNA_splicing" in tmp_list1) and (not "intergenic" in tmp_list1) and (not "ncRNA_intergenic" in tmp_list1):
                        print("现在的情况是15")
                        for i in range(len(tmp_list1)):
                            if tmp_list1[i] == "splicing" or "ncRNA_splicing":
                                geneby_splicing = geneby.split(";")[i]  # 获取对应位置的splicing的信息
                                if ")," in geneby_splicing:
                                    tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]  # 获取对应位置的splicing的信息
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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")

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
                        w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")
                else:
                    print("现在tem_list是1")
                    if tmp_list1[0] == "splicing" or tmp_list1[0] == "ncRNA_splicing":
                        if ")," in geneby:
                            tmp_splicing_list_pre = [x.split("(") for x in geneby.split("),")]  # 获取对应位置的splicing的信息
                            tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]  # 这个得到的list的长度也不一定是1，所以获取gene_name的时候还是会出现多个,并且可能重复
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
                    w1.write("\t".join(whole[0:8]) + "\t" + symbol_list_str + "\t" + desc_list_str + "\t" + "\t".join(whole[8:]) + "\n")


if __name__ == '__main__':
    snp_anno_new = "/mnt/ilustre/users/sanger-dev/sg-users/litangjian/extra1_20180810/extra1"
    id_name_desc = "/mnt/ilustre/users/sanger-dev/sg-users/litangjian/id_name_desc.txt"
    add_name_des(snp_anno_new, id_name_desc)

