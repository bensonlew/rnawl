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


class MirPreferAgent(Agent):
    """
    新miRNA预测
    """

    def __init__(self, parent):
        super(MirPreferAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {'default': '', 'type': 'string', 'name': 'index'},  # bowtie索引
            {'default': 1, 'type': 'int', 'name': 'mismatch'},  # 比对允许的错配
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
        self._memory = '10G'

    def end(self):
        super(MirPreferAgent, self).end()


class MirPreferTool(Tool):
    def __init__(self, config):
        super(MirPreferTool, self).__init__(config)
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/rna/miniconda3/bin/'
        self.bowtie_bash = "bioinfo/align/bowtie-1.2.3-linux-x86_64/bowtie"
        self.python = "miniconda2/bin/python"
        self.mir_prefer = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miR-PREFeR-master/"
        self.mir_prefer_config = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miR-PREFeR-master/example/config.example"
        self.bowtie_bash = '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.2.3-linux-x86_64/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/ViennaRNA-2.4.14/bin/'
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.samtools)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.rnafold)

    def run(self):
        """
        运行
        :return:
        """
        super(MirPreferTool, self).run()
        self.split_fasta()
        self.bowtie_align_reads()
        self.parse_config_file()
        self.run_mir_prefer()
        self.rna_fold(self.work_dir + "/mir_prefer/pdfs",
                      os.path.join(self.work_dir, "mir_prefer/novel_miRNA.precursor.fa"))
        self.ps_2_pdf()
        self.parse_mir_prefer_results()
        self.count_2_tpm()
        self.end()

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

    def split_fasta(self):
        self.logger.info("开始序列拆分")
        input_fa = self.option("input_fa").prop["path"]
        for seq_record in SeqIO.parse(input_fa, "fasta"):
            sample = re.search(r'(S\d+)_\d+_x\d+', seq_record.id).group(1)
            sample_fa = os.path.join(self.work_dir, sample + ".fa")
            with open(sample_fa, "a") as w:
                w.write(str(">" + seq_record.id + "\n" + seq_record.seq + "\n"))

    def bowtie_align_reads(self):
        self.logger.info("开始序列比对")
        self.samples = list()
        config_file = self.option("config").prop["path"]
        config = ConfigParser.ConfigParser()
        config.read(config_file)
        for item in config.items("NAME"):
            self.samples.append(item[0].upper())
        for sample in self.samples:
            cmd = "{}/bowtie -a -v {} {} -p 10 -f {} -S {}".format(self.bowtie_bash, self.option("mismatch"),
                                                                   self.option("index"),
                                                                   os.path.join(self.work_dir, sample + ".fa"),
                                                                   os.path.join(self.work_dir, sample + ".sam"))
            command_name = "bowtie_align_{}".format(sample.lower())
            command = self.add_command(command_name, cmd)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("bowtie比对成功")
            else:
                self.set_error("bowtie比对出错")

    def parse_config_file(self):
        output_dir = os.path.join(self.work_dir, "mir_prefer")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(os.path.join(self.work_dir, "mir_prefer"))
        os.link(self.mir_prefer_config, os.path.join(self.work_dir, os.path.basename(self.mir_prefer_config)))
        align_sam = list()
        for sam in sorted(glob.glob(os.path.join(self.work_dir + "/*.sam"))):
            align_sam.append(sam)
        with open("config.example", "r") as f, open("config.txt", "w") as w:
            for line in f:
                if line.startswith("#"):
                    continue
                elif line.startswith("PIPELINE_PATH"):
                    w.write("PIPELINE_PATH = {}/scripts".format(self.mir_prefer) + "\n")
                elif line.startswith("FASTA_FILE"):
                    w.write("FASTA_FILE =  {}".format(self.option("input_genome").prop['path']) + "\n")
                elif line.startswith("ALIGNMENT_FILE = "):
                    w.write("ALIGNMENT_FILE = {}".format(",".join(align_sam)) + "\n")
                elif line.startswith("PRECURSOR_LEN"):
                    w.write("PRECURSOR_LEN = 300" + "\n")
                elif line.startswith("READS_DEPTH_CUTOFF"):
                    w.write("READS_DEPTH_CUTOFF = 20" + "\n")
                elif line.startswith("NUM_OF_CORE"):
                    w.write("NUM_OF_CORE = 10" + "\n")
                elif line.startswith("OUTFOLDER"):
                    w.write("OUTFOLDER = {}".format(output_dir) + "\n")
                elif line.startswith("NAME_PREFIX"):
                    w.write("NAME_PREFIX = novel" + "\n")
                elif line.startswith("MAX_GAP"):
                    w.write("MAX_GAP = 300" + "\n")
                elif line.startswith("MIN_MATURE_LEN"):
                    w.write("MIN_MATURE_LEN = 18" + "\n")
                elif line.startswith("MAX_MATURE_LEN"):
                    w.write("MAX_MATURE_LEN = 24" + "\n")
                elif line.startswith("ALLOW_NO_STAR_EXPRESSION"):
                    w.write("ALLOW_NO_STAR_EXPRESSION = Y" + "\n")
                elif line.startswith("ALLOW_3NT_OVERHANG"):
                    w.write("ALLOW_3NT_OVERHANG = N" + "\n")
                elif line.startswith("CHECKPOINT_SIZE"):
                    w.write("CHECKPOINT_SIZE = 30" + "\n")
                else:
                    continue

    def run_mir_prefer(self):
        cmd = "{} {}/miR_PREFeR_modify.py -L -k -d pipeline {}".format(self.python, self.mir_prefer,
                                                                       os.path.join(self.work_dir, "config.txt"))
        command = self.add_command("mir_prefer", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mir_prefer脚本完成")
        else:
            self.set_error("运行mir_prefer脚本出错")

    def parse_mir_prefer_results(self):
        detail_csv = glob.glob(self.work_dir + "/mir_prefer/*.detail.csv")[0]
        novel_mirna_count = self.output_dir + "/novel_mirna_count.xls"
        novel_mirna_detail = self.output_dir + "/novel_mirna_detail.xls"
        novel_mature_seq = self.output_dir + "/novel_mature_seq.fa"
        novel_precursor_seq = self.output_dir + "/novel_precursor_seq.fa"
        novel_precursors_str = self.output_dir + "/novel_precursors.str"
        structure = self.work_dir + "/mir_prefer/pdfs/NovmiR.str"
        precursor_energy = dict()
        precursor_map = dict()
        with open(structure, "r") as f:
            i = 0
            for line in f:
                i += 1
                if i % 3 == 1:
                    id = ""
                    pos = ""
                    energy = ""
                if line.startswith(">"):
                    pos = line.strip().split(" ")[0].replace(">", "").replace(":", "_")
                    id = line.strip().split(" ")[2]
                if i % 3 == 0:
                    #  energy = re.match(r'\((.*)\)', line.strip().split(" ")[1:]).group(1)
                    #  有时候会出现特殊情况，如( -14.60)
                    energy = line.strip().split(" (")[1].split(")")[0].strip()
                    if float(energy) > 0:
                        energy = None
                if id and pos and energy and id not in precursor_energy:
                    precursor_energy[id] = energy
                    precursor_map[pos] = id
        with open(detail_csv, "r") as f, open(novel_mirna_detail, "w") as w1, open(novel_mirna_count, "w") as w2, open(
                novel_mature_seq, "w") as w3, open(novel_precursor_seq, "w") as w4, open(novel_precursors_str,
                                                                                         "w") as w5:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "energy"]
            w1.write("\t".join(head_list) + "\n")
            head = f.readline()
            sample_list = list()
            samples = head.strip("\n").split(", ")[9].split(",")
            for sam in samples:
                if sam:
                    self.logger.info(sam)
                    try:
                        sample = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=sam)
                    except:
                        sample = ""
                    self.logger.info(sample)
                    if sample and sample not in sample_list:
                        sample_list.append(sample)
            w2.write("miRNA\t" + "\t".join(sample_list) + "\n")
            f.readline()
            for line in f:
                items = line.strip("\n").split(",")
                counts = items[9:]
                sample_count = list()
                if len(counts) % 4 == 0:
                    for i, count in enumerate(counts):
                        if i % 4 == 1:
                            sample_count.append(int(count))
                miRNA_id = items[0].strip()
                seq_id = items[1].strip()
                precursor_start = items[2].strip()
                precursor_end = items[3].strip()
                strand = items[4].strip()
                precursor_seq = items[5].strip()
                precursor_str = items[6].strip()
                mature_seq = items[7].strip()
                position = str(seq_id) + "(" + strand + "):" + str(precursor_start) + "-" + str(precursor_end)
                w1.write("\t".join([miRNA_id.replace("-precursor", ""), mature_seq, str(
                    len(mature_seq)), miRNA_id, precursor_seq, str(len(precursor_seq)), position, str(
                    precursor_energy[miRNA_id])]) + "\n")
                w2.write(
                    miRNA_id.replace("-precursor", "") + "\t" + "\t".join(map(lambda x: str(x), sample_count)) + "\n")
                w3.write(">" + miRNA_id.replace("-precursor", "") + "\n" + mature_seq + "\n")
                w4.write(">" + miRNA_id + "\n" + precursor_seq + "\n")
                w5.write(">" + miRNA_id + "\n" + precursor_seq + "\n" + precursor_str + "\n")
        if os.path.exists(self.output_dir + "/structure_pdf"):
            shutil.rmtree(self.output_dir + "/structure_pdf")
        os.mkdir(self.output_dir + "/structure_pdf")
        pdf_files = glob.glob(self.work_dir + "/*_ss.ps.pdf")
        for file in pdf_files:
            base_name = precursor_map[os.path.basename(file).split("_ss.ps.pdf")[0]]
            os.link(file, self.output_dir + "/structure_pdf/" + base_name + ".pdf")
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))

    def count_2_tpm(self):
        novel_mirna_count = self.output_dir + "/novel_mirna_count.xls"
        novel_mirna_norm = self.output_dir + "/novel_mirna_norm.xls"
        novel_mirna_count_pd = pd.read_table(novel_mirna_count, index_col=0)
        for i in novel_mirna_count_pd.columns:
            sum_i = novel_mirna_count_pd[i].sum()
            novel_mirna_count_pd[i] = novel_mirna_count_pd[i].apply(lambda x: round(float(x*1000000/sum_i), 2))
        novel_mirna_count_pd.to_csv(novel_mirna_norm, sep="\t")

    def rna_fold(self, output_dir, input):
        self.logger.info("开始进行二级结构预测")
        pdfs = output_dir
        if os.path.exists(pdfs):
            shutil.rmtree(pdfs)
        else:
            os.mkdir(pdfs)
        cmd = "{}RNAfold < {} > {}".format(self.rnafold, input, os.path.join(pdfs, "NovmiR.str"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('rnafold计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('rnafold计算失败')
            self.set_error("rnafold运行失败")

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
            "id": "MirPrefer_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.mir_prefer",
            "instant": False,
            "options": dict(
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/filtered.fa",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                mismatch=1
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
