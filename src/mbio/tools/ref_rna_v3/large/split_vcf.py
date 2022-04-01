# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
import shutil
# import pandas as pd
__author__ = 'shicaiping'


class SplitVcfAgent(Agent):
    """
    split_fasta description
    """
    def __init__(self, parent):
        super(SplitVcfAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'vcf', 'type': 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('10')

    def end(self):
        super(SplitVcfAgent, self).end()


class SplitVcfTool(Tool):
    """
    split_fasta description
    """
    def __init__(self, config):
        self.vcftools = "/bioinfo/rna/vcftools_0.1.13/bin/vcftools"
        super(SplitVcfTool, self).__init__(config)

    def split_vcf(self):
        splits = list()
        chromosomes = list()
        size = 0
        genome = self.option("fasta").prop['path']
        vcf = self.option("vcf").prop['path']

        def check_chromosomes(chromosomes):
            f = []
            add_chrs = []
            les = 0
            for i in chromosomes:
                if (les + len(i)) > 131000:
                    les = 0
                    add_chrs.append(i)
                    f.append(add_chrs)
                    add_chrs=[]
                else:
                    les += len(i)
                    add_chrs.append(i)
            if add_chrs:
                f.append(add_chrs)
            return f


        # add_splits = check_chromosomes()
        for seq_record in SeqIO.parse(genome, 'fasta'):
            if float(len(seq_record.seq))/1024/1024 >= 50:
                splits.append(seq_record.id)
            else:
                if (size + float(len(seq_record.seq))/1024/1024) >= 50:
                    size = 0
                    chromosomes.append(seq_record.id)
                    new_chrs = check_chromosomes(chromosomes)
                    for chrs in new_chrs:
                        splits.append(chrs)
                    chromosomes = list()
                else:
                    size += float(len(seq_record.seq))/1024/1024
                    chromosomes.append(seq_record.id)

        if chromosomes:
            splits.append(chromosomes)
        if splits:
            with open(os.path.join(self.work_dir, "splits_file"),"w") as w:
                for s in splits:
                    try:
                        w.write(s + "\n")
                    except:
                        w.write("|".join(s) + "\n")
        if chromosomes:
            with open(os.path.join(self.work_dir, "chromosomes_file"),"w") as w:
                for s in chromosomes:
                    w.write(s+"\n")
        if not splits and len(chromosomes) >= 3:
            self.logger.info("染色体长度较短，按照染色体个数分成3份！")
            a = len(chromosomes)/3
            splits.append(chromosomes[0:a])
            splits.append(chromosomes[a:2 * a])
            splits.append(chromosomes[2 * a:])
        if splits:
            for i, chrs in enumerate(splits):
                vcf_output = os.path.join(self.output_dir, os.path.basename(vcf).replace(".vcf", "_{}.vcf".format(i)))
                if isinstance(chrs, list):
                    chr = ''
                    for tmp_chr in chrs:
                        chr += '--chr {} '.format(tmp_chr)
                else:
                    chr = '--chr {} '.format(chrs)
                cmd = "{} --vcf {} {} --recode --recode-INFO-all --out {}".format(self.vcftools, vcf, chr, vcf_output)
                command = self.add_command("extract_{}_vcf".format(i), cmd).run()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("提取染色体：{} vcf文件成功".format(i))
                else:
                    self.set_error("提取染色体：{} vcf文件失败".format(i))
        else:
            self.logger.info("染色体长度不符合拆分条件，按照不拆分运行！")
            os.link(vcf, os.path.join(self.output_dir, os.path.basename(vcf)))

    def run(self):
        super(SplitVcfTool, self).run()
        self.split_vcf()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SplitVcf" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.large.split_vcf",
            "instant": False,
            "options": dict(
                fasta="/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/plants/Brassica_napus/genoscope/dna/Brassica_napus_v4.1.chromosomes.fa",
                vcf="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/CallSnpIndelSentieon/VcfFilterGatk/output/final.vcf"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
