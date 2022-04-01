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
            {"name": "mismatch", "type": "int", "default": 0},  # 允许错误匹配
            {"name": "dir", "type": "string", "default": "dir"},  # 输出目录后缀
            {"name": "clean_fa", "type": "infile", "format": "small_rna.fasta"},  # 质控后的序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "miRNA_count", "type": "outfile", "format": "small_rna.common"},  # 已知miRNA表达count表
            {"name": "miRNA_norm", "type": "outfile", "format": "small_rna.common"},  # 已知miRNA表达norm表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完已知miRNA的过滤文件
            {"name": "mrd", "type": "outfile", "format": "small_rna.common"},  # 鉴定完已知miRNA的过滤文件
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
        if self.option("species").lower() == "null":
            raise OptionError("必须指定具体物种")
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
        self.quantifier = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        self.perl = '/program/perl-5.24.0/bin/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/'
        self.parse_arf = self.config.PACKAGE_DIR + "/small_rna/parse_arf.pl"
        self.get_exp = self.config.PACKAGE_DIR + "/small_rna/get_mireap_exp.pl"
        self.mature_mirna = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'mature_mirna']
        self.pre_mirna = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'pre_mirna']
        self.mature2hairpin = \
        Config().get_mongo_client(mtype='small_rna', ref=True)[Config().get_mongo_dbname('small_rna', ref=True)][
            'mature2hairpin']
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
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
        species = self.option("species").split(",")
        k = 0
        for i in species:
            dir = i
            k += 1
            mature_fa = os.path.join(self.work_dir, i + "_mature.fa")
            hairpin_fa = os.path.join(self.work_dir, i + "_hairpin.fa")
            mature_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(i, i)
            hairpin_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.hairpin.fa".format(i, i)
            if os.path.exists(mature_fa):
                os.remove(mature_fa)
            if os.path.exists(hairpin_fa):
                os.remove(hairpin_fa)
            os.link(mature_i, mature_fa)
            os.link(hairpin_i, hairpin_fa)
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

    def get_mature_fa(self, mirna_count):
        mature_files = glob.glob(self.work_dir + "/*_mature.fa")
        precursor_files = glob.glob(self.work_dir + "/*_hairpin.fa")
        mature_fa = self.work_dir + "/mature.fa"
        precursor_fa = self.work_dir + "/hairpin.fa"
        total_csv = self.work_dir + "/miRNAs_expressed_all_samples_total.csv"
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
        for specie in species:
            mirna_exp_files.append(self.work_dir + "/miRNAs_expressed_all_samples_{}.csv".format(specie))
            mrd_files.append(
                self.work_dir + "/expression_analyses/expression_analyses_{}/".format(specie) + "miRBase.mrd")
        ## 汇总所有表达量表，在植物的情况下用于novel miRNA预测
        mirna_total = os.path.join(self.work_dir, "miRNAs_expressed_all_samples_total.csv")
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
        ## 整理合并后的mrd文件
        mrd = os.path.join(self.work_dir, "miRBase.mrd")
        with open(mrd, "w") as w:
            for file in mrd_files:
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

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果目录")
        mirna_count = os.path.join(self.work_dir, "known_miR_count.xls")
        mirna_norm = os.path.join(self.work_dir, "known_miR_norm.xls")
        mirna_detail = os.path.join(self.work_dir, "known_mirna_detail")
        mrd = os.path.join(self.work_dir, "miRBase.mrd")
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_count.xls"))
        os.link(mirna_count, os.path.join(self.output_dir, "known_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_norm.xls"))
        os.link(mirna_norm, os.path.join(self.output_dir, "known_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "known_mirna_detail.xls"))
        if os.path.exists(os.path.join(self.output_dir, "miRBase.mrd")):
            os.remove(os.path.join(self.output_dir, "miRBase.mrd"))
        os.link(mrd, os.path.join(self.output_dir, "miRBase.mrd"))
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
        species = self.option("species").split(",")
        for specie in species:
            pdf_dir = self.work_dir + "/pdfs_{}".format(specie)
            if os.path.exists(pdf_dir):
                for file in os.listdir(pdf_dir):
                    pre_id = file.strip().split(".pdf")[0]
                    if pre_id in self.known_mirna_pre:
                        os.link(os.path.join(pdf_dir, file), self.output_dir + "/structure_pdf/" + file)
        ## 去掉PDF转换为PNG步骤
        '''
        self.logger.info("开始PDF转换为PNG")
        pdf_paths = [os.path.join(self.output_dir + '/structure_pdf', pdf) for pdf in
                     os.listdir(self.output_dir + '/structure_pdf') if re.match(r'^\S+\.pdf$', pdf.strip())]
        self.logger.info(pdf_paths)
        if len(pdf_paths) < 1:
            self.logger.info('未生成pdf文件：路径为 %s' % self.output_dir + '/structure_pdf')
        else:
            for i in range(0, len(pdf_paths)):
                pdf_paths_new = pdf_paths[i]
                re_string = re.compile('.pdf$')
                png_paths = re_string.sub(".png", pdf_paths_new)
                self.convert_pdf_to_png(i+1, pdf_paths_new,png_paths)
        '''
        clean_fa = os.path.join(self.work_dir, "filtered.fa")
        if os.path.exists(os.path.join(self.output_dir, "filtered.fa")):
            os.remove(os.path.join(self.output_dir, "filtered.fa"))
        os.link(clean_fa, os.path.join(self.output_dir, "filtered.fa"))
        self.option("filter_fa", os.path.join(self.output_dir, "filtered.fa"))
        self.option("miRNA_count", os.path.join(self.output_dir, "known_mirna_count.xls"))
        self.option("miRNA_norm", os.path.join(self.output_dir, "known_mirna_norm.xls"))
        self.option("mrd", os.path.join(self.output_dir, "miRBase.mrd"))
        self.logger.info("设置结果完成")


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "KnownMirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.known_mirna",
            "instant": False,
            "options": dict(
                species="dre",
                clean_fa="/mnt/ilustre/users/isanger/sg-users/shicaiping/miRNA/quantifier_test_animal/rfam_trimed.fa",
                config="/mnt/ilustre/users/isanger/sg-users/shicaiping/miRNA/quantifier_test_animal/Uniq.cfg.ini"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
