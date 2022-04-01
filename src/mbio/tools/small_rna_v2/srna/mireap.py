# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest
import shutil
import subprocess
import glob
from Bio import SeqIO
from skbio.parse.sequences import parse_fasta


class MireapAgent(Agent):
    """
        新miRNA预测
        """
    def __init__(self, parent):
        super(MireapAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "mapping_arf", "type": "infile", "format": "small_rna.common"},  # clean reads mapping到参考基因组的结果
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {"name": "min_count", "type": "int", "default": 10},  # min novel miRNA counts
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
        if not self.option("input_fa").is_set:
            raise OptionError("必须提供质控后的FASTA文件")
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
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MireapAgent, self).end()


class MireapTool(Tool):
    def __init__(self, config):
        super(MireapTool, self).__init__(config)
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.mireap = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mireap_0.2/bin/"
        self.parse_arf = self.config.PACKAGE_DIR + "/small_rna/"
        self.FFW1 = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mireap_0.2/lib/'
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/rna/miniconda3/bin/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2_0.1.3/ViennaRNA-2.4.14/bin/'
        self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
        self.lib_rna = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/arch/auto/RNA/'
        self.RNA = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/lib/'
        self.python = "program/Python/bin/python"
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(LD_LIBRARY_PATH=self.lib_rna)
        self.set_environ(PATH=self.mireap)
        self.set_environ(PATH=self.samtools)
        self.set_environ(PERL5LIB=self.FFW1)
        self.set_environ(PATH=self.parse_arf)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.randfold)
        self.set_environ(PATH=self.rnafold)
        self.set_environ(PERL5LIB=self.RNA)

    def run(self):
        super(MireapTool, self).run()
        self.get_mapfile()
        self.run_mireap()
        self.parse_mireap_results()
        self.rna_fold(self.work_dir + "/mireap/pdfs", os.path.join(self.work_dir, "mireap/novel_miR_pre.fa"))
        self.ps_2_pdf()
        self.set_mireap_results()
        self.end()

    def get_mapfile(self):
        self.logger.info("开始对mapping arf文件进行处理")
        seq_ids = set()
        arf = self.option("mapping_arf").prop["path"]
        input_fa = self.option("input_fa").prop["path"]
        filtered_fa = self.work_dir + "/filtered.fa"
        map = self.work_dir + "/map.txt"
        with open(filtered_fa, "w") as w:
            for seq_record in SeqIO.parse(input_fa, "fasta"):
                seq_ids.add(seq_record.id)
                short_id = seq_record.id.split("_x")[0]
                count = seq_record.id.split("_x")[1]
                w.write(">" + short_id + " " + str(count) + "\n")
                w.write(str(seq_record.seq) + "\n")
        with open(arf, "r") as f, open(map, "w") as w:
            for line in f:
                items = line.strip().split("\t")
                seq_id = items[0]
                count = int(seq_id.split("_x")[1])
                seq_short = seq_id.split("_x")[0]
                if seq_id in seq_ids and count >= self.option("min_count"):
                    map_info = "{}\t{}\t{}\t{}\n".format(items[5], items[7], items[8], items[10])
                    w.write(seq_short + "\t" + map_info)

    def parse_mireap_results(self):
        self.logger.info("开始对mireap结果文件进行处理")
        cmd = "{}perl {}parse_mireap_result_v2.pl -aln {} -gff {} -config {} -o {}".format(
            self.perl, self.mireap, os.path.join(self.work_dir, "mireap-Nov.aln"),
            os.path.join(self.work_dir, "mireap-Nov.gff"), self.option("config").prop["path"],
            os.path.join(self.work_dir, "mireap"))
        command = self.add_command("parse_mireap_results", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mireap结果文件处理完成")
        else:
            self.set_error("对mireap结果文件处理出错")

    def run_mireap(self):
        self.logger.info("开始进行novel miRNA预测")
        map = self.work_dir + "/map.txt"
        input_genome = self.option("input_genome").prop["path"]
        cmd = "perl {}mireap.pl -i {} -m {} -r {} -t Nov".format(self.mireap, self.work_dir + "/filtered.fa",
                                                                 map, input_genome)
        command = self.add_command("mireap", cmd)
        command.software_dir = "/usr/bin/"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mireap脚本完成")
        else:
            self.set_error("运行mireap脚本出错")

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

    def set_mireap_results(self):
        self.logger.info("开始设置mireap结果目录")
        ## novel miRNA二级结构图片
        if os.path.exists(os.path.join(self.output_dir, "structure_pdf")):
            shutil.rmtree(os.path.join(self.output_dir, "structure_pdf"))
        os.mkdir(os.path.join(self.output_dir, "structure_pdf"))
        pdf = glob.glob(self.work_dir + "/*.pdf")
        for file in pdf:
            os.link(file, os.path.join(self.output_dir, "structure_pdf", os.path.basename(file)))

        novel_mature_fa = {}
        novel_precursor_fa = {}
        mirna_info = self.work_dir + "/mireap/novel_miR_mature_infor.xls"
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail.xls")
        for seq_record in SeqIO.parse(self.work_dir + "/mireap/novel_miR_mature.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mireap/novel_miR_pre.fa", "fasta"):
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

        mrd = os.path.join(self.work_dir, "mireap-Nov.aln")
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(mrd, "r") as f1:
            for line in f1:
                if len(line.strip().split(" ")) >= 3:
                    if re.compile(r'S.*_\d+\s+').match(line.strip().split(" ")[1]):
                        mirna_id = line.strip().split(" ")[1] + "_" + line.strip().split(" ")[2]
                        a[mirna_id] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(self.work_dir + "/mireap/novel_miR_pre.tr", os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mireap/novel_miR_mature.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mireap/novel_miR_pre.fa", os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(self.work_dir + "/mireap/novel_miR_count.xls", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(self.work_dir + "/mireap/novel_miR_norm.xls", os.path.join(self.output_dir, "novel_mirna_norm.xls"))
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
            "id": "Mireap_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.mireap",
            "instant": False,
            "options": dict(
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
