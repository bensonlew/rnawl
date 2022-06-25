# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import shutil
import subprocess
import glob
from Bio import SeqIO
import pandas as pd
from skbio.parse.sequences import parse_fasta


class Mirdeep2Agent(Agent):
    """
    新miRNA预测
    """
    def __init__(self, parent):
        super(Mirdeep2Agent, self).__init__(parent)
        options = [
            {"name": "species", "type": "string", "default": ""},  # 具体物种
            {"name": "ref_mature", "type": "string", "default": ""},  # 参考mature序列
            {"name": "ref_hairpin", "type": "string", "default": ""},  # 参考mature序列
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "mapping_arf", "type": "infile", "format": "small_rna.common"},  # clean reads mapping到参考基因组的结果
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {"name": "mirdeep2_version", "type": "string", "default": "0.1.3"},  # 0.0.5 | 0.1.3
            {"name": "database", "type": "string", "default": "mirbase"},  # 参考miRNA数据库
        ]
        self.add_option(options)
        self.step.add_steps("novel_mirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 30

    def stepstart(self):
        self.step.novel_mirna.start()
        self.step.update()

    def stepfinish(self):
        self.step.novel_mirna.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("input_fa").is_set:
            raise OptionError("必须提供质控后的FASTA文件")
        if not self.option("mapping_arf").is_set:
            raise OptionError("必须提供clean reads mapping到参考基因组的结果文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        with open(self.option("input_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("质控后的序列文件格式不对")
        if self.option("species") in ["Auto", "All"]:
            self.option("species", self.option("species").lower())
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(Mirdeep2Agent, self).end()


class Mirdeep2Tool(Tool):
    def __init__(self, config):
        super(Mirdeep2Tool, self).__init__(config)
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/rna/miniconda3/bin/'
        self.python = "miniconda2/bin/python"
        self.mirdeep = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2_0.1.3/mirdeep2-0.1.3/bin/"
        self.mirdeep_lib = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2_0.1.3/mirdeep2-0.1.3/lib/perl5/"
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/ViennaRNA-2.4.14/bin/'
        self.squid = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/squid-1.9g/'
        self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
        self.quantifier = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2_0.1.3/mirdeep2-0.1.3/bin/"
        self.get_exp = self.config.PACKAGE_DIR + "/small_rna/get_mireap_exp.pl"
        self.bowtie_bash = '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PERL5LIB=self.mirdeep_lib)
        self.set_environ(PATH=self.samtools)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.mirdeep)
        self.set_environ(PATH=self.randfold)
        self.set_environ(PATH=self.squid)
        self.set_environ(PATH=self.rnafold)
        self.set_environ(PATH=self.quantifier)
        if self.option("mirdeep2_version") == "0.0.5":
            self.RNA = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/lib/'
            self.mirdeep = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
            self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Progs/'
            self.lib_rna = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/arch/auto/RNA/'
            self.squid = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/squid-1.9g/'
            self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
            self.quantifier = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
            self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
            self.set_environ(LD_LIBRARY_PATH=self.lib_rna)
            self.set_environ(PATH=self.mirdeep)
            self.set_environ(PATH=self.randfold)
            self.set_environ(PATH=self.squid)
            self.set_environ(PATH=self.rnafold)
            self.set_environ(PERL5LIB=self.RNA)
            self.set_environ(PATH=self.quantifier)

    def run(self):
        """
        运行
        :return:
        """
        super(Mirdeep2Tool, self).run()
        self.get_mirna()
        self.run_remove_white_space_in_id()
        self.run_replace_nonatcgun()
        self.run_mirdeep2()
        if self.mature_fa != "none" and self.hairpin_fa != "none":
            self.parse_mirdeep_result()
        else:
            self.quantify_novel()
            self.get_novel_detail()
        self.end()

    def quantify_novel(self):
        mature_fa = glob.glob(self.work_dir + "/mirna_results*/novel_mature_*.fa")[0]
        hairpin_fa = glob.glob(self.work_dir + "/mirna_results*/novel_pres_*.fa")[0]
        clean_fa = self.option("input_fa").prop["path"]
        cmd = "{}perl {}quantifier.pl -p {} -m {} -g {} -r {} -y {} -d".format(self.perl, self.quantifier,
                                                                            hairpin_fa, mature_fa, 1,
                                                                            clean_fa, "novel")
        command = self.add_command("quantify_novel_mirna", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行quantify脚本完成")
        else:
            self.set_error("运行quantify脚本出错")
        mrd = self.work_dir + "/expression_analyses/expression_analyses_novel/miRBase.mrd"
        exp = self.work_dir + "/miRNAs_expressed_all_samples_novel.csv"
        filter_fa = os.path.join(self.work_dir, "filtered.fa")
        self.parse_mrd(mrd, clean_fa, filter_fa)
        count = self.work_dir + "/novel_miR_count.xls"
        norm = self.work_dir + "/novel_miR_norm.xls"
        cmd = "{}perl {} -exp {} -mrd {} -config {} -count {} -norm {}".format(self.perl, self.get_exp, exp, mrd,
                                                            self.option("config").prop["path"], count, norm)
        command = self.add_command("get_tpm", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行get_mireap_exp脚本完成")
        else:
            self.set_error("运行get_mireap_exp脚本出错")

    def get_novel_detail(self):
        mature_fa = glob.glob(self.work_dir + "/mirna_results*/novel_mature_*.fa")[0]
        hairpin_fa = glob.glob(self.work_dir + "/mirna_results*/novel_pres_*.fa")[0]
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail")
        result_bed = glob.glob(self.work_dir + "/result_*.bed")
        precursors_str = glob.glob(self.work_dir + "/mirdeep_runs/run_*/tmp/precursors.str")[0]
        novel_precursors_str = os.path.join(self.work_dir, "precursors.str")
        novel_mature_fa = {}
        novel_precursor_fa = {}
        precursors_str_dict = {}
        for seq_record in SeqIO.parse(mature_fa, "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(hairpin_fa, "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(precursors_str, "r") as f:
            seq_id = {}
            id = ""
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    id = line.strip().split(" ")[0].split(">")[1]
                    seq_id[id] = line.strip()
                    precursors_str_dict[id] = ""
                else:
                    if id != "":
                        precursors_str_dict[id] += line
        with open(result_bed[0], "r") as f1, open(mirna_detail, "w") as w1, open(novel_precursors_str, "w") as w2:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "score"]
            w1.write("\t".join(head_list) + "\n")
            for line in f1:
                items = line.strip().split("\t")
                if re.match(r'novel', items[3]):
                    id = items[3].split(":")[1]
                    pos = items[0] + "(" + items[5] + "):" + items[6] + "-" + items[7]
                    score = items[4]
                    if id in novel_mature_fa.keys():
                        mature_seq = novel_mature_fa[id]
                        precursor_seq = novel_precursor_fa[id]
                        w1.write("\t".join(
                            [id, str(mature_seq), str(len(mature_seq)), id, str(precursor_seq), str(len(precursor_seq)),
                             pos, score]) + "\n")
                        w2.write(">" + id + "\n" + precursors_str_dict[id])

        mirna_count = os.path.join(self.work_dir, "novel_miR_count.xls")
        mirna_norm = os.path.join(self.work_dir, "novel_miR_norm.xls")
        mrd = glob.glob(self.work_dir + "/mirdeep_runs/run_*/output.mrd")[0]
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(mirna_count, os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(mirna_norm, os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna.mrd")):
            os.remove(os.path.join(self.output_dir, "novel_mirna.mrd"))
        os.link(mrd, os.path.join(self.output_dir, "novel_mirna.mrd"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(novel_precursors_str, os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(mature_fa, os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(hairpin_fa, os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(self.output_dir + "/structure_pdf"):
            shutil.rmtree(self.output_dir + "/structure_pdf")
        os.mkdir(self.output_dir + "/structure_pdf")
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        pdf_files = glob.glob(self.work_dir + "/pdfs_*/*.pdf")
        for file in pdf_files:
            base_name = os.path.basename(file).split(".pdf")[0]
            if base_name in novel_mature_fa.keys():
                if not os.path.exists(self.output_dir + "/structure_pdf/" + os.path.basename(file)):
                    os.link(file, self.output_dir + "/structure_pdf/" + os.path.basename(file))

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
        self.option("filter_fa", output_fa)

    def get_mirna(self):
        if self.option("database").lower() == "none":
            self.mature_fa = "none"
            self.hairpin_fa = "none"
            self.other_mature_fa = "none"
        elif self.option("species") in ['auto', 'all']:
            if self.option("ref_mature") and self.option("ref_hairpin"):
                self.mature_fa = self.option("ref_mature")
                self.hairpin_fa = self.option("ref_hairpin")
                self.other_mature_fa = "none"
        else:
            self.logger.info("species: {}".format(self.option("species")))
            species = self.option("species").split(",")
            self.mature_fa = os.path.join(self.work_dir, species[0] + "_mature.fa")
            self.hairpin_fa = os.path.join(self.work_dir, species[0] + "_hairpin.fa")
            if self.option("database").lower() == "mirbase":
                mature_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(species[0],
                                                                                                         species[0])
                hairpin_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.hairpin.fa".format(species[0],
                                                                                                           species[0])
            elif self.option("database").lower() == "pmiren":
                mature_i = self.config.SOFTWARE_DIR + "/database/PmiREN/species/{}/{}_mature.fa".format(species[0],
                                                                                                        species[0])
                hairpin_i = self.config.SOFTWARE_DIR + "/database/PmiREN/species/{}/{}_hairpin.fa".format(species[0],
                                                                                                          species[0])
            if os.path.exists(self.mature_fa):
                os.remove(self.mature_fa)
            if os.path.exists(self.hairpin_fa):
                os.remove(self.hairpin_fa)
            os.link(mature_i, self.mature_fa)
            os.link(hairpin_i, self.hairpin_fa)

            self.other_mature_fa = os.path.join(self.work_dir, "other_mature.fa")
            ## 预测新miRNA时，把其他所有的物种作为参考
            with open(self.other_mature_fa, "w") as w1:
                mature_fa_file = []
                for i in species:
                    if self.option("database").lower() == "mirbase":
                        mature_fa_file.extend(
                            [self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(i, i)])
                    elif self.option("database").lower() == "pmiren":
                        mature_fa_file.extend(
                            [self.config.SOFTWARE_DIR + "/database/PmiREN/species/{}/{}_mature.fa".format(i, i)])
                    self.logger.info(mature_fa_file)
                if self.option("database").lower() == "mirbase":
                    other_fa_file = glob.glob(self.config.SOFTWARE_DIR + "/database/mirbase/species/*/*.mature.fa")
                elif self.option("database").lower() == "pmiren":
                    other_fa_file = glob.glob(self.config.SOFTWARE_DIR + "/database/PmiREN/species/*/*_mature.fa")
                for file in other_fa_file:
                    if file not in mature_fa_file:
                        with open(file, "r") as r1:
                            lines = r1.readlines()
                            w1.writelines(lines)

    def parse_mirdeep_result(self):
        self.logger.info("处理mirdeep2预测结果")
        file = os.path.basename(glob.glob(self.work_dir + "/miRNAs_expressed_all_samples_*.csv")[0])
        id = re.search("miRNAs_expressed_all_samples_(.*).csv", file).group(1)
        cmd = "{}perl {}convent_mirdeep_result_modify.pl -id {} -o mirdeep -config {}".format(self.perl, self.mirdeep,
                                                                                              id, self.option(
                "config").prop["path"])
        command = self.add_command("parse_mirdeep_result", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("parse_mirdeep_result完成")
        else:
            self.set_error("parse_mirdeep_result出错")

        infile1 = self.work_dir + "/mirdeep/predicted_mature_sumary_table.xls"
        with open(infile1, "r") as f, open(self.work_dir + "/novel_miR_count.xls", "w") as w:
            for line in f:
                if re.search(r'#', line):
                    print(line)
                    if re.search(r'provisional_id', line):
                        items = line.strip().split("\t")
                        w.write("miRNA_ID" + "\t" + "\t".join(items[2:]) + "\n")
                    elif re.search(r'miRBase_miRNAs', line):
                        break
                    else:
                        continue
                elif line.startswith("total"):
                    continue
                else:
                    items = line.strip().split("\t")
                    w.write(items[0] + "\t" + "\t".join(items[2:]) + "\n")
        infile2 = self.work_dir + "/mirdeep/novo_mature_normalization.xls"
        cmd2 = """sed 's/#provisional_id/miRNA_ID/' %s > novel_miR_norm.xls""" % (infile2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info('生成novel_miR_norm文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成novel_miR_norm文件失败')
            self.set_error("生成novel_miR_norm文件失败")

        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail")
        result_bed = glob.glob(self.work_dir + "/result_*.bed")
        precursors_str = glob.glob(self.work_dir + "/mirdeep_runs/run_*/tmp/precursors.str")[0]
        novel_precursors_str = os.path.join(self.work_dir, "precursors.str")
        novel_mature_fa = {}
        novel_precursor_fa = {}
        precursors_str_dict = {}
        for seq_record in SeqIO.parse(self.work_dir + "/mirdeep/novo_mature_seq.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mirdeep/novo_precursor_seq.fa", "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(precursors_str, "r") as f:
            seq_id = {}
            id = ""
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    id = line.strip().split(" ")[0].split(">")[1]
                    seq_id[id] = line.strip()
                    precursors_str_dict[id] = ""
                else:
                    if id != "":
                        precursors_str_dict[id] += line
        with open(result_bed[0], "r") as f1, open(mirna_detail, "w") as w1, open(novel_precursors_str, "w") as w2:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "score"]
            w1.write("\t".join(head_list) + "\n")
            for line in f1:
                items = line.strip().split("\t")
                if re.match(r'novel', items[3]):
                    id = items[3].split(":")[1]
                    pos = items[0] + "(" + items[5] + "):" + items[6] + "-" + items[7]
                    score = items[4]
                    if id in novel_mature_fa.keys():
                        mature_seq = novel_mature_fa[id]
                        precursor_seq = novel_precursor_fa[id]
                        w1.write("\t".join(
                            [id, str(mature_seq), str(len(mature_seq)), id, str(precursor_seq), str(len(precursor_seq)),
                             pos, score]) + "\n")
                        w2.write(">" + id + "\n" + precursors_str_dict[id])

        mirna_count = os.path.join(self.work_dir, "novel_miR_count.xls")
        mirna_norm = os.path.join(self.work_dir, "novel_miR_norm.xls")
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail")
        mrd = glob.glob(self.work_dir + "/mirdeep_runs/run_*/output.mrd")[0]
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(mirna_count, os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(mirna_norm, os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna.mrd")):
            os.remove(os.path.join(self.output_dir, "novel_mirna.mrd"))
        os.link(mrd, os.path.join(self.output_dir, "novel_mirna.mrd"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(novel_precursors_str, os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mirdeep/novo_mature_seq.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mirdeep/novo_precursor_seq.fa",
                os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(self.output_dir + "/structure_pdf"):
            shutil.rmtree(self.output_dir + "/structure_pdf")
        os.mkdir(self.output_dir + "/structure_pdf")
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        pdf_files = glob.glob(self.work_dir + "/mirdeep/pdfs_*/*.pdf")
        for file in pdf_files:
            base_name = os.path.basename(file).split(".pdf")[0]
            if base_name in novel_mature_fa.keys():
                os.link(file, self.output_dir + "/structure_pdf/" + os.path.basename(file))
        ## 整理鉴定完novel miRNA剩余的reads
        mrd = os.path.join(self.output_dir, "novel_mirna.mrd")
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(mrd, "r") as f1:
            for line in f1:
                if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                    a[line.strip().split()[0]] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)

    def run_remove_white_space_in_id(self):
        input_genome = self.option("input_genome").prop["path"]
        self.parse_genome = os.path.join(self.work_dir, os.path.basename(input_genome))
        cmd = "perl {}remove_white_space_in_id.pl {} > {}".format(self.mirdeep, input_genome, self.parse_genome)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('remove_white_space_in_id计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('remove_white_space_in_id计算失败')
            self.set_error("remove_white_space_in_id运行失败")

    def run_replace_nonatcgun(self):
        temp_mature_fa = os.path.join(self.work_dir, "temp_mature.fa")
        temp_hairpin_fa = os.path.join(self.work_dir, "temp_hairpin.fa")
        temp_other_mature_fa = os.path.join(self.work_dir, "temp_other_mature.fa")
        temp_genome = self.parse_genome + "_re"
        if self.mature_fa != "none":
            cmd1 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.mature_fa, temp_mature_fa)
            self.logger.info(cmd1)
            try:
                subprocess.check_output(cmd1, shell=True)
                self.logger.info('replace_nonatcgun计算完成')
            except subprocess.CalledProcessError:
                self.set_error("replace_nonatcgun运行失败")
        else:
            temp_mature_fa = "none"
        if self.hairpin_fa != "none":
            cmd2 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.hairpin_fa, temp_hairpin_fa)
            self.logger.info(cmd2)
            try:
                subprocess.check_output(cmd2, shell=True)
                self.logger.info('replace_nonatcgun计算完成')
            except subprocess.CalledProcessError:
                self.set_error("replace_nonatcgun运行失败")
        else:
            temp_hairpin_fa = "none"
        if self.other_mature_fa != "none":
            cmd3 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.other_mature_fa,
                                                                    temp_other_mature_fa)
            self.logger.info(cmd3)
            try:
                subprocess.check_output(cmd3, shell=True)
                self.logger.info('replace_nonatcgun计算完成')
            except subprocess.CalledProcessError:
                self.set_error("replace_nonatcgun运行失败")
        else:
            temp_other_mature_fa = "none"
        cmd4 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.parse_genome, temp_genome)
        self.logger.info(cmd4)
        try:
            subprocess.check_output(cmd4, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        os.rename(temp_genome, self.parse_genome)
        self.mature_fa = temp_mature_fa
        self.hairpin_fa = temp_hairpin_fa
        self.other_mature_fa = temp_other_mature_fa

    def run_mirdeep2(self):
        """
        动物novel miRNA预测
        """
        tmp_files = glob.glob(self.work_dir + "/miRNAs_expressed_all_samples_*.csv") + glob.glob(
            self.work_dir + "/result_*.bed") + glob.glob(self.work_dir + "/mirdeep_runs/run_*/tmp/precursors.str")
        for file in tmp_files:
            os.remove(file)
        self.logger.info("开始进行novel miRNA预测")
        input_fa = self.option("input_fa").prop["path"]
        mapping_arf = self.option("mapping_arf").prop["path"]
        cmd = "{}perl {}miRDeep2_modify.pl {} {} {} {} {} {}".format(self.perl, self.mirdeep, input_fa,
                                                                     self.parse_genome, mapping_arf, self.mature_fa,
                                                                     self.other_mature_fa, self.hairpin_fa)
        command = self.add_command("mirdeep2", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mirdeep脚本完成")
        else:
            self.set_error("运行mirdeep脚本出错")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Mirdeep2_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.mirdeep2",
            "instant": False,
            "options": dict(
                species="all",
                database='pmiren',
                ref_mature="none",
                ref_hairpin="none",
                #ref_mature="/mnt/ilustre/users/sanger-dev/app/database/mirbase/species/all_plant/all_plant_mature.fa",
                #ref_hairpin="/mnt/ilustre/users/sanger-dev/app/database/mirbase/species/all_plant/all_plant_hairpin.fa",
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/filtered.fa",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/reads_vs_genome.arf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
