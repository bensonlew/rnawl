# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import ConfigParser
import glob
import os
import re
import shutil
import subprocess
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from skbio.parse.sequences import parse_fasta


class KnownMirnaAgent(Agent):
    """
    已知miRNA鉴定
    """

    def __init__(self, parent):
        super(KnownMirnaAgent, self).__init__(parent)
        options = [
            {"name": "species", "type": "string", "default": ""},  # 具体物种,多选时分号分隔
            {"name": "organism_list", "type": "infile", 'format': "small_rna.common"},  # known mirna鉴定的物种列表
            {"name": "mismatch", "type": "int", "default": 0},  # 允许错误匹配
            {"name": "dir", "type": "string", "default": "dir"},  # 输出目录后缀
            {"name": "clean_fa", "type": "infile", "format": "small_rna.fasta"},  # 质控后的序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "miRNA_count", "type": "outfile", "format": "small_rna.common"},  # 已知miRNA表达count表
            {"name": "miRNA_norm", "type": "outfile", "format": "small_rna.common"},  # 已知miRNA表达norm表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完已知miRNA的过滤文件
            {"name": "mrd", "type": "outfile", "format": "small_rna.common"},  # 鉴定完已知miRNA的过滤文件
            {"name": "mirdeep2_version", "type": "string", "default": "0.0.5"},  # 0.0.5 | 0.1.3
            {"name": "database", "type": "string", "default": "mirbase"},  # 参考miRNA数据库
            {"name": "ref_mature", "type": "outfile", "format": "small_rna.fasta"},  # 参考mature序列
            {"name": "ref_hairpin", "type": "outfile", "format": "small_rna.fasta"},  # 参考mature序列
        ]
        self.add_option(options)
        self.step.add_steps("known_mirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.known_mirna.start()
        self.step.update()

    def stepfinish(self):
        self.step.known_mirna.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if self.option("species").lower() == "all_animal" and self.option("database").lower() == "pmiren":
            raise OptionError("PmiREN数据库暂不支持全动物物种分析")
        if self.option("species").lower() == "auto" and not self.option("organism_list").is_set:
            raise OptionError("必须指定物种列表")
        if not self.option("clean_fa").is_set:
            raise OptionError("必须提供质控后的FASTA文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        with open(self.option("clean_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("质控后的序列文件格式不对")

        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "已知miRNA鉴定结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r"pdfs", "", "已知miRNA二级结构文件目录"],
            ["known_mirna_count.xls", "xls", "已知miRNA定量count表"],
            ["known_mirna_norm.xls", "xls", "已知miRNA定量norm表"],
            ["filtered.fa", "", "过滤FASTA文件"],
        ])
        super(KnownMirnaAgent, self).end()


class KnownMirnaTool(Tool):
    def __init__(self, config):
        super(KnownMirnaTool, self).__init__(config)
        if self.option("mirdeep2_version") == "0.0.5":
            self.quantifier = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
            self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        elif self.option("mirdeep2_version") == "0.1.3":
            self.quantifier = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2_0.1.3/mirdeep2-0.1.3/bin/"
            self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.python = 'miniconda2/bin/'
        self.perl = '/program/perl-5.24.0/bin/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/'
        self.parse_arf = self.config.PACKAGE_DIR + "/small_rna/parse_arf.pl"
        self.get_exp = self.config.PACKAGE_DIR + "/small_rna/get_mireap_exp.pl"
        self.merge_mirna = self.config.PACKAGE_DIR + "/small_rna/merge_mirna.py"
        self.mature_mirna = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'mature_mirna']
        self.pre_mirna = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'pre_mirna']
        self.mature2hairpin = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'mature2hairpin']
        self.mirbase = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')]['mirbase']
        self.pmiren = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')]['pmiren']
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.quantifier)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.rnafold)

    def run(self):
        """
        运行
        :return:
        """
        super(KnownMirnaTool, self).run()
        self.known_mirna()
        self.get_tpm()
        self.set_output()
        self.end()

    def known_mirna(self):
        """
        已知miRNA鉴定
        """
        self.logger.info("开始进行已知miRNA鉴定")
        mismatch = self.option("mismatch")
        if self.option("species").lower() == "auto":
            self.get_mirna()
            mature_fa = os.path.join(self.work_dir, "auto_mature.fa")
            hairpin_fa = os.path.join(self.work_dir, "auto_hairpin.fa")
            clean_fa = self.option("clean_fa").prop["path"]
            cmd = "{}perl {}quantifier.pl -p {} -m {} -g {} -r {} -y {}".format(self.perl, self.quantifier,
                                                                                hairpin_fa, mature_fa, mismatch,
                                                                                clean_fa, "auto")
            command = self.add_command("known_mirna_auto", cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("运行quantify脚本完成")
            else:
                self.set_error("运行quantify脚本出错")
            mrd = self.work_dir + "/expression_analyses/expression_analyses_auto/miRBase.mrd"
            filter_fa = os.path.join(self.work_dir, "filtered.fa")
            self.parse_mrd(mrd, clean_fa, filter_fa)
        else:
            species = self.option("species").split(",")
            species_list = list()
            k = 0
            for i in species:
                i = i.strip()
                if self.option("database").lower() == "pmiren":
                    if len(i) >= 3:
                        try:
                            organism = self.pmiren.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                else:
                    if len(i) >= 3:
                        try:
                            organism = self.mirbase.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                self.logger.info("specie: {}".format(i))
                dir = i
                mature_fa = os.path.join(self.work_dir, i + "_mature.fa")
                hairpin_fa = os.path.join(self.work_dir, i + "_hairpin.fa")
                if self.option("database").lower() == 'mirbase':
                    mature_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(i, i)
                    hairpin_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.hairpin.fa".format(i, i)
                else:
                    mature_i = self.config.SOFTWARE_DIR + "/database/PmiREN/species/{}/{}_mature.fa".format(i, i)
                    hairpin_i = self.config.SOFTWARE_DIR + "/database/PmiREN/species/{}/{}_hairpin.fa".format(i, i)
                if os.path.exists(mature_fa):
                    os.remove(mature_fa)
                if os.path.exists(hairpin_fa):
                    os.remove(hairpin_fa)
                if os.path.exists(mature_i) and os.path.exists(hairpin_i) and i not in species_list:
                    k += 1
                    species_list.append(i)
                    os.link(mature_i, mature_fa)
                    os.link(hairpin_i, hairpin_fa)
                else:
                    self.logger.info("{}或{}不存在".format(mature_i, hairpin_i))
                    continue
                if k == 1:
                    clean_fa = self.option("clean_fa").prop["path"]
                    cmd = "{}perl {}quantifier.pl -p {} -m {} -g {} -r {} -y {}".format(self.perl, self.quantifier,
                                                                                        hairpin_fa, mature_fa, mismatch,
                                                                                        clean_fa, dir)
                    command = self.add_command("known_mirna_{}".format(i), cmd).run()
                    self.wait()
                    if command.return_code == 0:
                        self.logger.info("运行quantify脚本完成")
                    else:
                        self.set_error("运行quantify脚本出错")
                    mrd = self.work_dir + "/expression_analyses/expression_analyses_" + i + "/miRBase.mrd"
                    filter_fa = os.path.join(self.work_dir, "filtered.fa")
                    self.parse_mrd(mrd, clean_fa, filter_fa)
                else:
                    clean_fa = os.path.join(self.work_dir, "filtered.fa")
                    cmd = "{}perl {}quantifier.pl -p {} -m {} -g {} -r {} -y {}".format(self.perl, self.quantifier,
                                                                                        hairpin_fa, mature_fa, mismatch,
                                                                                        clean_fa, dir)
                    command = self.add_command("known_mirna_{}".format(i), cmd).run()
                    self.wait()
                    if command.return_code == 0:
                        self.logger.info("运行quantify脚本完成")
                    else:
                        self.set_error("运行quantify脚本出错")
                    mrd = self.work_dir + "/expression_analyses/expression_analyses_" + i + "/miRBase.mrd"
                    input_fa = os.path.join(self.work_dir, "filtered_tmp.fa")
                    filter_fa = os.path.join(self.work_dir, "filtered.fa")
                    subprocess.call('mv {} {}'.format(clean_fa, input_fa), shell=True)
                    self.parse_mrd(mrd, input_fa, filter_fa)
            self.option("species", ",".join(species_list))

    def get_mature_fa(self, mirna_count):
        mature_files = glob.glob(self.work_dir + "/*_mature.fa")
        precursor_files = glob.glob(self.work_dir + "/*_hairpin.fa")
        mature_fa = self.work_dir + "/mature.fa"
        precursor_fa = self.work_dir + "/hairpin.fa"
        expression_analyses_total = os.path.join(self.work_dir, "expression_analyses/expression_analyses_total")
        total_csv = os.path.join(expression_analyses_total, "miRNAs_expressed_all_samples_total.csv")
        mature_id = dict()
        precursor_id = dict()
        with open(mirna_count, "r") as f1:
            f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                mature_id[items[0]] = 1
        with open(total_csv, "r") as f1:
            f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                precursor_id[items[2]] = 1
        seq_id = {}
        seq_sequence = {}
        for file in mature_files:
            with open(file, "r") as f:
                id = ""
                sequence = ""
                for line in f:
                    if line.startswith(">"):
                        id = line.strip().split(" ")[0].split(">")[1]
                        seq_id[id] = line.strip()
                        seq_sequence[id] = ""
                    else:
                        if id != "":
                            sequence = line.strip()
                            seq_sequence[id] += sequence
        with open(mature_fa, "w") as w:
            for id in seq_id.keys():
                if id in mature_id.keys():
                    w.write(seq_id[id] + "\n")
                    w.write(seq_sequence[id] + "\n")
        seq_id = {}
        seq_sequence = {}
        for file in precursor_files:
            with open(file, "r") as f:
                id = ""
                sequence = ""
                for line in f:
                    if line.startswith(">"):
                        id = line.strip().split(" ")[0].split(">")[1]
                        seq_id[id] = line.strip()
                        seq_sequence[id] = ""
                    else:
                        if id != "":
                            sequence = line.strip()
                            seq_sequence[id] += sequence
        with open(precursor_fa, "w") as w:
            for id in seq_id.keys():
                if id in precursor_id.keys():
                    w.write(seq_id[id] + "\n")
                    w.write(seq_sequence[id] + "\n")
        self.logger.info("已知mirna序列提取完毕")

    def parse_mrd(self, mrd, input_fa, output_fa):
        """
        对定量的结果进行处理
        """
        self.logger.info("开始生成过滤后的fasta文件")
        a = dict()
        with open(mrd, "r") as f1:
            for line in f1:
                if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                    a[line.strip().split()[0]] = 1
        with open(output_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.logger.info("fasta过滤文件处理完毕")

    def parse_config(self, file=None, section=None, name=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config.get(section, name)

    def convert_pdf_to_png(self, num, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + olds + ' ' + news
        self.logger.info(cmd)
        command = self.add_command('convert_{}'.format(num), cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            pass
        else:
            self.set_error("PDF转换PNG出错!")

    def get_tpm(self):
        """
        生成均一化的表达量tpm表
        """
        species = self.option("species").split(",")
        self.logger.info(species)
        mirna_exp_files = []
        mrd_files = []
        arf_files = []
        expression_analyses_total = os.path.join(self.work_dir, "expression_analyses/expression_analyses_total")
        if os.path.exists(expression_analyses_total):
            shutil.rmtree(expression_analyses_total)
        os.mkdir(expression_analyses_total)
        for specie in species:
            specie = specie.strip().lower()
            mirna_exp_files.append(self.work_dir + "/miRNAs_expressed_all_samples_{}.csv".format(specie))
            mrd_files.append(
                self.work_dir + "/expression_analyses/expression_analyses_{}/".format(specie) + "miRBase.mrd")
            arf_files.append(self.work_dir + "/expression_analyses/expression_analyses_{}/".format(
                specie) + "{}_mature_mapped.arf".format(specie))
        # 汇总所有表达量表，在植物的情况下用于novel miRNA预测
        mirna_total = os.path.join(expression_analyses_total, "miRNAs_expressed_all_samples_total.csv")
        with open(mirna_total, "w") as w:
            num = 0
            for file in mirna_exp_files:
                num += 1
                with open(file, "r") as exp:
                    if num <= 1:
                        for line in exp:
                            w.write(line)
                    else:
                        header = exp.readline()
                        for line in exp:
                            w.write(line)
        # 整理合并后的mrd文件
        mrd = os.path.join(expression_analyses_total, "miRBase.mrd")
        with open(mrd, "w") as w:
            for file in mrd_files:
                self.logger.info(file)
                with open(file, "r") as r:
                    lines = r.readlines()
                    w.writelines(lines)
        # 整理合并后的arf文件
        arf = os.path.join(expression_analyses_total, "mature_mapped.arf")
        with open(arf, "w") as w:
            for file in arf_files:
                self.logger.info(file)
                with open(file, "r") as r:
                    lines = r.readlines()
                    w.writelines(lines)
        cmd = "{}perl {} -exp {} -mrd {} -config {}".format(self.perl, self.get_exp, mirna_total, mrd,
                                                            self.option("config").prop["path"])
        command = self.add_command("get_tpm", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行get_mireap_exp脚本完成")
        else:
            self.set_error("运行get_mireap_exp脚本出错")

        # 生成mirna详情表
        mirna_detail = os.path.join(self.work_dir, "known_mirna_detail")
        mirna_count = os.path.join(self.work_dir, "known_miR_count.xls")
        self.known_mirna_pre = dict()
        with open(mirna_count, "r") as f, open(mirna_detail, "w") as w:
            header = f.readline()
            array = ["miRNA_name", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position"]
            w.write("\t".join(array) + "\n")
            for line in f:
                items = line.strip().split("\t")
                try:
                    miRNA_seq = self.mature_mirna.find_one({"_id": items[0]})["seq"]
                except:
                    miRNA_seq = "-"
                try:
                    pres = self.mature2hairpin.find({"mature": items[0]})
                    pres = list(pres)
                    if len(list(pres)) == 0:
                        array = [items[0], miRNA_seq, str(len(miRNA_seq)), "-", "-", "-", "-"]
                        w.write("\t".join(array) + "\n")
                    else:
                        for pre in pres:
                            pre_id = pre["hairpin"]
                            self.known_mirna_pre.update({pre_id: pre_id})
                            try:
                                pre_seq = self.pre_mirna.find_one({"_id": pre_id})["seq"]
                            except:
                                pre_seq = "-"
                            try:
                                pre_position = self.pre_mirna.find_one({"_id": pre_id})["position"]
                            except:
                                pre_position = "-"
                            array = [items[0], miRNA_seq, str(len(miRNA_seq)), pre_id, pre_seq, str(len(pre_seq)),
                                     pre_position]
                            w.write("\t".join(array) + "\n")
                except:
                    array = [items[0], miRNA_seq, str(len(miRNA_seq)), "-", "-", "-", "-"]
                    w.write("\t".join(array) + "\n")

    def get_mirna(self):
        self.logger.info("获取物种列表合并去冗余的mirna序列")
        priority_list = self.option("organism_list").prop["path"]
        if self.option("database").lower() == "mirbase":
            mature2hairpin = self.config.SOFTWARE_DIR + "/database/mirbase/mature2hairpin.txt"
        else:
            mature2hairpin = self.config.SOFTWARE_DIR + "/database/PmiREN/mature2hairpin.txt"
        cmd = "{}python {} -priority_list {} -mature2hairpin {} -database {} -prefix {}".format(self.python,
                                                                                                self.merge_mirna,
                                                                                                priority_list,
                                                                                                mature2hairpin,
                                                                                                self.option(
                                                                                                    "database").lower(),
                                                                                                "auto")
        command = self.add_command("get_mirna", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行{}完成".format(cmd))
        else:
            self.set_error("运行{}出错".format(cmd))

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果目录")
        mirna_count = os.path.join(self.work_dir, "known_miR_count.xls")
        mirna_norm = os.path.join(self.work_dir, "known_miR_norm.xls")
        mirna_detail = os.path.join(self.work_dir, "known_mirna_detail")
        expression_analyses_total = os.path.join(self.work_dir, "expression_analyses/expression_analyses_total")
        mrd = os.path.join(expression_analyses_total, "miRBase.mrd")
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_count.xls"))
        os.link(mirna_count, os.path.join(self.output_dir, "known_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_norm.xls"))
        os.link(mirna_norm, os.path.join(self.output_dir, "known_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "known_mirna_detail.xls"))
        if os.path.exists(os.path.join(self.output_dir, "known_mirna.mrd")):
            os.remove(os.path.join(self.output_dir, "known_mirna.mrd"))
        os.link(mrd, os.path.join(self.output_dir, "known_mirna.mrd"))
        self.get_mature_fa(mirna_count)
        mature_fa = self.work_dir + "/mature.fa"
        hairpin_fa = self.work_dir + "/hairpin.fa"
        if os.path.exists(os.path.join(self.output_dir, "mature.fa")):
            os.remove(os.path.join(self.output_dir, "mature.fa"))
        os.link(mature_fa, os.path.join(self.output_dir, "mature.fa"))
        if os.path.exists(os.path.join(self.output_dir, "hairpin.fa")):
            os.remove(os.path.join(self.output_dir, "hairpin.fa"))
        os.link(hairpin_fa, os.path.join(self.output_dir, "hairpin.fa"))
        if os.path.exists(self.output_dir + "/structure_pdf"):
            shutil.rmtree(self.output_dir + "/structure_pdf")
        os.mkdir(self.output_dir + "/structure_pdf")
        pdf_dirs = glob.glob(self.work_dir + "/pdfs_*")
        for pdf_dir in pdf_dirs:
            if os.path.exists(pdf_dir):
                for file in os.listdir(pdf_dir):
                    pre_id = file.strip().split(".pdf")[0]
                    if pre_id in self.known_mirna_pre:
                        os.link(os.path.join(pdf_dir, file), self.output_dir + "/structure_pdf/" + file)
        clean_fa = os.path.join(self.work_dir, "filtered.fa")
        if os.path.exists(os.path.join(self.output_dir, "filtered.fa")):
            os.remove(os.path.join(self.output_dir, "filtered.fa"))
        os.link(clean_fa, os.path.join(self.output_dir, "filtered.fa"))
        self.option("filter_fa", os.path.join(self.output_dir, "filtered.fa"))
        self.option("miRNA_count", os.path.join(self.output_dir, "known_mirna_count.xls"))
        self.option("miRNA_norm", os.path.join(self.output_dir, "known_mirna_norm.xls"))
        self.option("mrd", os.path.join(self.output_dir, "known_mirna.mrd"))
        if self.option("species").lower() == "all_plant":
            self.option("ref_mature", os.path.join(self.work_dir, "all_plant_mature.fa"))
            self.option("ref_hairpin", os.path.join(self.work_dir, "all_plant_hairpin.fa"))
        elif self.option("species").lower() == "all_animal":
            self.option("ref_mature", os.path.join(self.work_dir, "all_animal_mature.fa"))
            self.option("ref_hairpin", os.path.join(self.work_dir, "all_animal_hairpin.fa"))
        elif self.option("species").lower() == "auto":
            self.option("ref_mature", os.path.join(self.work_dir, "auto_mature.fa"))
            self.option("ref_hairpin", os.path.join(self.work_dir, "auto_hairpin.fa"))
        else:
            ref_mature = self.work_dir + "/ref_mature.fa"
            ref_hairpin = self.work_dir + "/ref_hairpin.fa"
            if os.path.exists(ref_mature):
                os.remove(ref_mature)
            if os.path.exists(ref_hairpin):
                os.remove(ref_hairpin)
            mature_fas = glob.glob(self.work_dir + "/*_mature.fa")
            hairpin_fas = glob.glob(self.work_dir + "/*_hairpin.fa")

            with open(ref_mature, "w") as w:
                for file in sorted(mature_fas):
                    with open(file, "r") as f:
                        for line in f:
                            w.write(line)
            with open(ref_hairpin, "w") as w:
                for file in sorted(hairpin_fas):
                    with open(file, "r") as f:
                        for line in f:
                            w.write(line)
            self.option("ref_mature", ref_mature)
            self.option("ref_hairpin", ref_hairpin)
        self.logger.info("设置结果完成")


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "KnownMirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.known_mirna",
            "instant": False,
            "options": dict(
                species="ath,cme,csa",
                clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/uniq.fasta",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                #organism_list="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/orgnism_list",
                database="pmiren"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
