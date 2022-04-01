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


class Mirdp2Agent(Agent):
    """
    新miRNA预测
    """

    def __init__(self, parent):
        super(Mirdp2Agent, self).__init__(parent)
        options = [
            {"name": "other", "type": "infile", "format": "small_rna.fasta"},  # 其他mature序列
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "known_miRNA_express", "type": "infile", "format": "small_rna.common"},  # 已知miRNA定量结果
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {"name": "method", "type": "string", "default": ""},  # mirdeep2, mireap, mirdp2, mir_prefer
            {"name": "mirdeep2_version", "type": "string", "default": "0.1.3"},  # 0.0.5 | 0.1.3
            {'default': '', 'type': 'string', 'name': 'index'},  # bowtie索引
            {'default': 1, 'type': 'float', 'name': 'rpm'},  # mirdeep-p2参数，用于测序数据过滤
            {'default': 1, 'type': 'int', 'name': 'mismatch'},  # mirdeep-p2参数，比对允许的错配
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
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        with open(self.option("input_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("质控后的序列文件格式不对")
        if not self.option("index"):
            raise OptionError("需要提供索引文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        file_size = float(os.path.getsize(self.option('input_fa').prop['path'])) / 1024 / 1024
        memory = int(file_size/50) + 20
        self._memory = '{}G'.format(memory)

    def end(self):
        super(Mirdp2Agent, self).end()


class Mirdp2Tool(Tool):
    def __init__(self, config):
        super(Mirdp2Tool, self).__init__(config)
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/rna/miniconda3/bin/'
        self.bowtie_bash = "bioinfo/align/bowtie-1.2.3-linux-x86_64/bowtie"
        self.mirdp2_bash = "bioinfo/miRNA/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline_modify.bash"
        self.mirdp2_scripts = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miRDP2-v1.1.4/scripts/"
        self.python = "program/Python/bin/python"
        self.bowtie_bash = '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/ViennaRNA-2.4.14/bin/'
        self.squid = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/squid-1.9g/'
        self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
        self.get_exp = self.config.PACKAGE_DIR + "/small_rna/get_mireap_exp.pl"
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.samtools)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.randfold)
        self.set_environ(PATH=self.squid)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.rnafold)


    def run(self):
        """
        运行
        :return:
        """
        super(Mirdp2Tool, self).run()
        mirdp2_output = self.work_dir + "/mirdp2"
        if not os.path.exists(os.path.join(mirdp2_output, "filtered", "filtered_filter_P_prediction.bed")):
            self.run_mirdp2()
        self.parse_mirdp2_results()
        self.set_mirdp2_results()
        self.end()

    def run_mirdp2(self):
        """
        miRDP2 novel miRNA identification
        """
        input_genome = self.option("input_genome").prop["path"]
        mirdp2_output = self.work_dir + "/mirdp2"
        if os.path.exists(mirdp2_output):
            shutil.rmtree(mirdp2_output)
        os.mkdir(mirdp2_output)
        cmd = "{} -g {} -x {} -i {} -o {} -R {} -M {}".format(self.mirdp2_bash, input_genome, self.option("index"),
                                                              self.option("input_fa").prop["path"], mirdp2_output,
                                                              self.option("rpm"), self.option("mismatch"))

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
        if os.path.exists(mirdp2_result):
            shutil.rmtree(mirdp2_result)
        os.mkdir(mirdp2_result)
        cmd = "{}perl {}parse_mirdp2_result.pl -prediction {} -config {} -o {}".format(
            self.perl, self.mirdp2_scripts, prediction_file, self.option("config").prop["path"], mirdp2_result)
        command = self.add_command("parse_mirdp2_results", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mirdp2结果文件处理完成")
        else:
            self.set_error("对mirdp2结果文件处理出错")

        with open(prediction_file, "r") as f:
            for line in f:
                pre_miRNA = str(line.strip().split("\t")[3])
                ps = pre_miRNA + "_ss.ps"
                pdf = pre_miRNA + ".pdf"
                if os.path.exists(self.work_dir + "/" + ps):
                    cmd1 = "ps2pdf {} {}".format(ps, pdf)
                else:
                    continue
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
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.tr",
                os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_mature.fa",
                os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.fa",
                os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_count.xls",
                os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_norm.xls",
                os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Mirdp2_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.mirdp2",
            "instant": False,
            "options": dict(
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/filtered.fa",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                rpm=1,
                mismatch=1
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
