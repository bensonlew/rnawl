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
from skbio.parse.sequences import parse_fasta


class Mirdp2Agent(Agent):
    """
    新miRNA预测
    """

    def __init__(self, parent):
        super(Mirdp2Agent, self).__init__(parent)
        options = [
            {"name": "category", "type": "string", "default": ""},  # 物种分类，Animal or Plant
            {"name": "species", "type": "string", "default": ""},  # 具体物种
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "mapping_arf", "type": "infile", "format": "small_rna.common"},  # clean reads mapping到参考基因组的结果
            {"name": "mrd", "type": "infile", "format": "small_rna.common"},  # 已知miRNA定量结果文件
            {"name": "known_miRNA_express", "type": "infile", "format": "small_rna.common"},  # 已知miRNA定量结果
            {"name": "clean_fa", "type": "infile", "format": "small_rna.fasta"},  # 用来mapping到参考基因组的fasta序列
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {'default': '', 'type': 'string', 'name': 'index'},  # bowtie索引
            {'default': 'mireap', 'type': 'string', 'name': 'software'},  # 植物预测miRNA所用软件
        ]
        self.add_option(options)
        self.step.add_steps("novel_mirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

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
        if self.option("category").lower() not in ["animal", "plant"]:
            raise OptionError("物种分类只能为Animal/Plant")
        if self.option("species") == "null":
            raise OptionError("必须指定具体物种")
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

        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '30G'

    def end(self):
        super(Mirdp2Agent, self).end()


class Mirdp2Tool(Tool):
    def __init__(self, config):
        super(Mirdp2Tool, self).__init__(config)
        self.mirdp2_bash = "bioinfo/miRNA/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline.bash"
        self.mirdp2_scripts = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miRDP2-v1.1.4/scripts/"
        self.mirdeep = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
        self.mireap = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mireap_0.2/bin/"
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Progs/'
        self.parse_arf = self.config.PACKAGE_DIR + "/small_rna/"
        self.FFW1 = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mireap_0.2/lib/'
        self.RNA = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/lib/'
        self.lib_rna = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/arch/auto/RNA/'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        self.squid = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/squid-1.9g/'
        self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
        self.set_environ(LD_LIBRARY_PATH=self.lib_rna)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        # self.set_environ(PATH=self.perl_old)
        self.set_environ(PERL5LIB=self.FFW1)
        self.set_environ(PERL5LIB=self.RNA)
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.set_environ(PATH=self.parse_arf)
        self.set_environ(PATH=self.mirdeep)
        self.set_environ(PATH=self.randfold)
        self.set_environ(PATH=self.squid)
        self.set_environ(PATH=self.mireap)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.rnafold)
        self.set_environ(PATH=self.mirdp2_bash)
        self.set_environ(PATH=self.mirdp2_scripts)

    def run(self):
        """
        运行
        :return:
        """
        super(Mirdp2Tool, self).run()
        self.get_mirna()
        self.run_mirdp2()
        self.parse_mirdp2_results()
        self.set_mirdp2_results()
        self.end()

    def get_mirna(self):
        species = self.option("species").split(",")
        self.mature_fa = os.path.join(self.work_dir, species[0] + "_mature.fa")
        self.hairpin_fa = os.path.join(self.work_dir, species[0] + "_hairpin.fa")
        mature_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(species[0], species[0])
        hairpin_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.hairpin.fa".format(species[0],
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
                mature_fa_file.extend(
                    [self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(i, i)])
                self.logger.info(mature_fa_file)
            other_fa_file = glob.glob(self.config.SOFTWARE_DIR + "/database/mirbase/species/*/*.mature.fa")
            for file in other_fa_file:
                if file not in mature_fa_file:
                    with open(file, "r") as r1:
                        lines = r1.readlines()
                        w1.writelines(lines)

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

    def parse_mapping_file(self):
        self.logger.info("开始对mapping arf文件进行处理")
        cmd = "perl {}parse_arf.pl -rfam {} -mrd {} -arf {}".format(self.parse_arf,
                                                                    self.option("clean_fa").prop["path"],
                                                                    self.option("mrd").prop["path"],
                                                                    self.option("mapping_arf").prop["path"])
        command = self.add_command("parse_arf", cmd)
        command.software_dir = "/usr/bin/"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mapping arf文件处理完成")
        else:
            self.set_error("对mapping arf文件处理出错")

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
        cmd1 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.mature_fa, temp_mature_fa)
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        cmd2 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.hairpin_fa, temp_hairpin_fa)
        self.logger.info(cmd2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        cmd3 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.other_mature_fa,
                                                                temp_other_mature_fa)
        self.logger.info(cmd3)
        try:
            subprocess.check_output(cmd3, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        self.mature_fa = temp_mature_fa
        self.hairpin_fa = temp_hairpin_fa
        self.other_mature_fa = temp_other_mature_fa

    def run_mirdp2(self):
        """
        植物novel miRNA预测
        """
        input_genome = self.option("input_genome").prop["path"]
        mirdp2_output = self.work_dir + "/mirdp2"
        os.mkdir(mirdp2_output)
        cmd = "{} -g {} -x {} -i {} -o {}".format(self.mirdp2_bash, input_genome, self.option("index"),
                                                            self.option("input_fa").prop["path"], mirdp2_output)
        command = self.add_command("mirdp2", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mirdp2脚本完成")
        else:
            self.set_error("运行mirdp2脚本出错")

    def parse_mirdp2_results(self):
        self.logger.info("开始对mirdp2结果文件进行处理")
        prediction_file = glob.glob(self.work_dir + "/mirdp2/*/*filter_P_prediction")[0]
        mirdp2_result = os.path.join(self.work_dir, "mirdp2_result")
        os.mkdir(mirdp2_result)
        cmd = "{}perl {}parse_mirdp2_result.pl -exp {} -prediction {} -config {} -o {}".format(
            self.perl, self.mirdp2_scripts, self.option("known_miRNA_express").prop["path"],
            prediction_file, self.option("config").prop["path"], mirdp2_result)
        command = self.add_command("parse_mirdp2_results", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mirdp2结果文件处理完成")
        else:
            self.set_error("对mirdp2结果文件处理出错")

        with open(prediction_file, "r") as f:
            for line in f:
                pre_miRNA = str(line.strip().split("\t")[3])[0:12]
                ps = pre_miRNA + "_ss.ps"
                pdf = pre_miRNA + ".pdf"
                cmd1 = "ps2pdf {} {}".format(ps, pdf)
                try:
                    subprocess.check_output(cmd1, shell=True)
                    self.logger.info('ps2pdf转换完成')
                except:
                    self.set_error("ps2pdf转换失败")

    def set_mirdp2_results(self):
        self.logger.info("开始设置mirdp2结果目录")
        ## novel miRNA二级结构图片
        if os.path.exists(os.path.join(self.output_dir, "structure_pdf")):
            shutil.rmtree(os.path.join(self.output_dir, "structure_pdf"))
        if os.path.exists(os.path.join(self.work_dir, "structure_ps")):
            shutil.rmtree(os.path.join(self.work_dir, "structure_ps"))
        os.mkdir(os.path.join(self.output_dir, "structure_pdf"))
        os.mkdir(os.path.join(self.work_dir, "structure_ps"))
        pdf = glob.glob(self.work_dir + "/*.pdf")
        ps = glob.glob(self.work_dir + "/*_ss.ps")
        for file in pdf:
            shutil.move(file, os.path.join(self.output_dir, "structure_pdf"))
        for file in ps:
            shutil.move(file, os.path.join(self.work_dir, "structure_ps"))

        novel_mature_fa = {}
        novel_precursor_fa = {}
        mirna_info = self.work_dir + "/mirdp2_result/novel_miR_mature_infor.xls"
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail.xls")
        for seq_record in SeqIO.parse(self.work_dir + "/mirdp2_result/novel_miR_mature.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mirdp2_result/novel_miR_pre.fa", "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(mirna_info, "r") as f1, open(mirna_detail, "w") as w1:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "energy"]
            w1.write("\t".join(head_list) + "\n")
            headline = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                mat_id = items[0]
                pre_id = items[5]
                pos = items[6] + "(" + items[9] + "):" + items[7] + "-" + items[8]
                energy = items[11]
                mature_seq = novel_mature_fa[mat_id]
                precursor_seq = novel_precursor_fa[pre_id]
                w1.write("\t".join(
                    [mat_id, str(mature_seq), str(len(mature_seq)), pre_id, str(precursor_seq), str(len(precursor_seq)),
                     pos, energy]) + "\n")

        prediction_file = glob.glob(self.work_dir + "/mirdp2/*/*filter_P_prediction")[0]
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(prediction_file, "r") as f1:
            for line in f1:
                reads = line.strip().split("\t")[11].split(";")
                for read in reads:
                    if read not in a:
                        a[read] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.tr", os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_mature.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.fa", os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_count.xls", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_norm.xls", os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))

    def rna_fold(self):
        self.logger.info("开始进行二级结构预测")
        pdfs = self.work_dir + "/mireap/pdfs"
        if os.path.exists(pdfs):
            shutil.rmtree(pdfs)
        else:
            os.mkdir(pdfs)
        cmd = "{}RNAfold < {} > {}".format(self.rnafold, os.path.join(self.work_dir, "mireap/novel_miR_pre.fa"),
                                           os.path.join(pdfs, "NovmiR.str"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('rnafold计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('rnafold计算失败')
            self.set_error("rnafold运行失败")

    def ps_2_pdf(self):
        ps = glob.glob(self.work_dir + "/*.ps")
        for file in ps:
            file1 = file + ".pdf"
            cmd1 = "ps2pdf {} {}".format(file, file1)
            try:
                subprocess.check_output(cmd1, shell=True)
                self.logger.info('ps2pdf转换完成')
            except:
                self.set_error("ps2pdf转换失败")

    def parse_config(self, file=None, section=None, name=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config.get(section, name)


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "NovelMirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.mirdp2",
            "instant": False,
            "options": dict(
                # species="sly",
                # category="Animal",
                # input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel.fa",
                # input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/filtered.fa",
                # config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/qc_file.config",
                # mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/reads_vs_genome.arf",
                # mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRBase.mrd",
                # known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRNAs_expressed_all_samples_total.csv",
                # index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel_index",
                # species="sly",
                # category="Plant",
                # software="mireap",
                # input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel.fa",
                # input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/filtered.fa",
                # clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/uniq.fasta",
                # config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/qc_file.config",
                # mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/reads_vs_genome.arf",
                # mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRBase.mrd",
                # known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRNAs_expressed_all_samples_total.csv",
                # index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel_index",
                species="ath",
                category="Plant",
                software="mirdp2",
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/uniq.fasta",
                clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/uniq.fasta",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/reads_vs_genome.arf",
                mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/miRBase.mrd",
                known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/miRNAs_expressed_all_samples_total.csv",
                index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
