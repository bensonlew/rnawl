# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last modified by shicaiping at 20180517

import os
import shutil
import re
import subprocess
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.express.single_sample import *
import unittest

class TranscriptAbstractAgent(Agent):
    """
    提取参考基因组最长序列，作为基因注释的输入文件
    """
    def __init__(self, parent):
        super(TranscriptAbstractAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因组fasta文件
            {"name": "ref_genome_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因组gtf文件
            {"name": "ref_genome_gff", "type": "infile", "format": "gene_structure.gff3"},  # 参考基因组gff文件
            {"name": "trans_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 输出做注释的转录本序列
            {"name": "gene_file", "type": "outfile", "format": "ref_rna_v2.gene_list"},  # 输出最长转录本
            {"name": "trans2gene", "type": "outfile", "format": "ref_rna_v2.common"},  # gene2trans
            {"name": "length_file", "type": "outfile", "format": "annotation.cog.cog_list"}  # 输出注释转录本序列的长度
        ]
        self.add_option(options)
        self.step.add_steps("Transcript")
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Transcript.start()
        self.step.update()

    def step_end(self):
        self.step.Transcript.finish()
        self.step.update()

    def check_option(self):
        if not self.option("ref_genome").is_set:
            raise OptionError("请设置参考基因组custom文件", code = "33708501")
        if self.option("ref_genome_gtf").is_set or self.option("ref_genome_gff").is_set:
            pass
        else:
            raise OptionError("请设置参考基因组gtf文件或gff文件", code = "33708502")

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(TranscriptAbstractAgent, self).end()

class TranscriptAbstractTool(Tool):
    def __init__(self, config):
        super(TranscriptAbstractTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"
        self.long_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/"
        self.python_path = "program/Python/bin/"
        self.length = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/fastalength"

    def run_gffread(self):
        old_fasta = self.option("ref_genome").prop["path"]  # 统一按自定义方式传参考基因组
        fasta = self.work_dir + "/" + os.path.basename(self.option("ref_genome").prop["path"])
        if os.path.exists(fasta):
            os.remove(fasta)
        os.link(old_fasta, fasta)
        gtf = self.option("ref_genome_gtf").prop["path"]
        cmd = "{}gffread {} -g {} -w exons.fa".format(self.gffread_path, gtf, fasta)
        self.logger.info("开始运行cufflinks的gffread，合成、提取exons")
        command = self.add_command("gffread", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("exons提取完成")
            self.option('trans_fa', self.work_dir + '/exons.fa')
        else:
            self.set_error("运exons提取出错", code = "33708503")

    def run_exons_length(self):
        """exons的长度"""
        exon_path = os.path.join(self.work_dir, "exons.fa")
        length_path = self.output_dir + "/exons_length.txt"
        cmd = "{} {} > {}".format(self.length, exon_path, length_path)
        self.logger.info("开始提取exons的长度")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("提取exons的长度完成")
        except subprocess.CalledProcessError:
            self.set_error("提取exons的长度失败", code = "33708504")
        self.option("length_file", length_path)

    def run_long_transcript(self):
        exon_path = os.path.join(self.work_dir, "exons.fa")
        cmd = "{}python {}annotation_longest.py -i {}".format(self.python_path, self.long_path, exon_path)
        self.logger.info("提取最长序列")
        command = self.add_command("the_longest", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运提取最长序列完成")
        else:
            self.set_error("提取最长序列出错", code = "33708505")
        output1 = os.path.join(self.work_dir, "exons.fa")
        if os.path.exists(self.output_dir + "/exons.fa"):
            os.remove(self.output_dir + "/exons.fa")
        os.link(output1, self.output_dir + "/exons.fa")
        output2 = os.path.join(self.work_dir, "the_longest_exons.fa")
        if os.path.exists(self.output_dir + "/the_longest_exons.fa"):
            os.remove(self.output_dir + "/the_longest_exons.fa")
        os.link(output2, self.output_dir + "/the_longest_exons.fa")

    def get_gene_list(self):
        output_path = self.work_dir + "/output/the_longest_exons.fa"
        gene_list_path = self.work_dir + '/output/gene_list.txt'
        gene_lists = []
        with open(output_path, 'rb') as f, open(gene_list_path, 'wb') as w:
            lines = f.readlines()
            for line in lines:
                m = re.match(r">(.+) gene=(.+)", line)
                if m:
                    trans_name = m.group(1)
                    if trans_name not in gene_lists:
                        w.write(trans_name + '\n')
                        gene_lists.append(trans_name)
                else:
                    n = re.match(r">(.+) transcript:(.+)", line)
                    if n:
                        trans_name = n.group(1)
                        if trans_name not in gene_lists:
                            w.write(trans_name + '\n')
                            gene_lists.append(trans_name)
        self.option('gene_file', gene_list_path)

    def gene2trans(self):
        trans2gene_tmp = os.path.join(self.work_dir, 'trans2gene_tmp')
        trans2gene = os.path.join(self.work_dir, 'trans2gene')
        gtf(self.option('ref_genome_gtf').prop['path'], trans2gene_tmp)
        with open(trans2gene_tmp, "r") as r, open(trans2gene, "w") as w:
            for line in r:
                line = line.strip().split()
                w.write(line[1] + "\t" + line[0] + "\n")
        self.option('trans2gene', self.work_dir + "/trans2gene")
        self.logger.info("提取gene_transcript成功！")

    def run(self):
        super(TranscriptAbstractTool, self).run()
        self.run_gffread()
        self.gene2trans()
        #self.run_exons_length()
        #self.run_long_transcript()
        #self.get_gene_list()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data5 = {
            "id": "NewTrnasFa_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "medical_transcriptome.transcript_abstract",
            "instant": False,
            "options": dict(
                ref_genome_gtf = "/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/RefrnaAssemble/output/NewTranscripts/ref_and_new.gtf",
                ref_genome = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # ref_genome_gtf="/mnt/ilustre/users/sanger-dev/workspace/20180606/Refrna_tsg_30361/output/assembly/NewTranscripts/ref_and_new.gtf",
                # ref_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
            )
        }
        wsheet = Sheet(data=data5)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()