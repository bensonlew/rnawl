# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# version 1.0

import os
import re
import glob
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.gene_structure.gff3 import Gff3File
import unittest


class IntronExonAgent(Agent):
    def __init__(self, parent):
        super(IntronExonAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考序列文件
            {"name": "gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考注释文件
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},  # 参考注释文件
            {"name": "exon", "type": "outfile", "format": "gene_structure.bed"},# 外显子bed注释文件
            {"name": "intron", "type": "outfile", "format": "gene_structure.bed"},# 内含子bed注释文件
            {"name": "exon_fa", "type": "outfile", "format": "small_rna.fasta"},# 外显子序列
            {"name": "intron_fa", "type": "outfile", "format": "small_rna.fasta"},# 内含子序列
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option("gtf").is_set:
            raise OptionError("必须设置参数gtf")
        if not self.option("gtf").is_set:
            if not self.option("gff").is_set:
                raise OptionError("gff和gtf中必须有一个作为参数传入")

    def set_resource(self):
        self._cpu = 4
        self._memory = '{}G'.format(int(os.path.getsize(self.option("input_genome").prop['path']) / 1024.0 ** 3 + 20))

    def end(self):
        super(IntronExonAgent, self).end()

class IntronExonTool(Tool):
    def __init__(self, config):
        super(IntronExonTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.ref_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.bedtools = "/bioinfo/rna/bedtools2-master/bin/"

    def exon_bed(self):
        self.logger.info("正在提取exon注释bed文件")
        if self.option("gff").is_set:
            origin_gff_path = self.option("gff").prop["path"]
            gff_name = os.path.split(origin_gff_path)[1]
            self.logger.info("gff的名称为{}".format(gff_name))
            new_gff_path = os.path.join(self.work_dir, gff_name)
            if os.path.exists(new_gff_path):
                os.remove(new_gff_path)
            os.link(origin_gff_path, new_gff_path)
            gff = Gff3File()
            gff.set_path(new_gff_path)
            gff.set_gtf_file(new_gff_path + ".gtf")
            if os.path.exists(new_gff_path + ".gtf"):
                os.remove(new_gff_path + ".gtf")
            gff.set_gffread_path(self.config.SOFTWARE_DIR + "/bioinfo/rna/cufflinks-2.2.1/gffread")
            gff.to_gtf()
            if os.path.exists(new_gff_path + ".gtf.bed"):
                os.remove(new_gff_path + ".gtf.bed")
            gff.set_gtf2bed_path(self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/gtf2bed.py")
            gff.gtf_to_bed(new_gff_path + ".gtf")
        else:
            new_gtf = self.work_dir + "/" + os.path.basename(self.option("gtf").prop["path"])
            if os.path.exists(new_gtf):
                os.remove(new_gtf)
            os.link(self.option("gtf").prop["path"], new_gtf)
            self.option("gtf").set_path(new_gtf)
            self.option("gtf").to_bed()
        self.logger.info("提取exon注释bed文件完成")

    def intron_bed(self):
        '''
        Extract Intron regions from bed file (must be 12-column).  output is 6-column Tab separated bed file, each row represents one intron
        '''
        self.logger.info("正在提取intron注释bed文件")
        exon_bed = self.option("gtf").prop["path"] + ".bed"
        intron_bed = self.option("gtf").prop["path"] + ".intron.bed"
        with open (exon_bed, "r") as f, open (intron_bed, "w") as w:
            for line in f:
                try:
                    if line.startswith('#'):continue
                    if line.startswith('track'):continue
                    if line.startswith('browser'):continue
                    fields = line.split("\t")
                    chrom = fields[0]
                    tx_start = int(fields[1])
                    tx_end = int(fields[2])
                    gene_name = fields[3]
                    strand = fields[5].replace(" ","_")
                    cds_start = int(fields[6])
                    cds_end = int(fields[7])
                    if int(fields[9] ==1):continue
                    exon_starts = map(int, fields[11].rstrip( ',\n' ).split( ',' ))
                    exon_starts = map((lambda x: x + tx_start ), exon_starts)
                    exon_ends = map(int, fields[10].rstrip( ',\n' ).split( ',' ))
                    exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
                    intron_start = exon_ends[:-1]
                    intron_end=exon_starts[1:]
                    if (strand == '-'):
                        intron_num=len(intron_start)
                        for st,end in zip(intron_start,intron_end):
                            w.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + gene_name + "_intron_" + str(intron_num) + "\t0\t" + strand + '\n')
                            intron_num -= 1
                    else:
                        intron_num = 1
                        for st,end in zip(intron_start,intron_end):
                            w.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + gene_name + "_intron_" + str(intron_num) + "\t0\t" + strand + '\n')
                            intron_num += 1
                except:
                    self.logger.info("提取intron注释bed文件失败")
        self.logger.info("提取intron注释bed文件成功")

    def get_exon(self):
        exon = self.option("gtf").prop["path"] + ".bed"
        cmd = self.bedtools + "bedtools getfasta -fi %s -bed %s -fo %s -name" % (self.genome_fasta, exon, self.work_dir + "/" + self.ref_name + ".exon.fa")
        command = self.add_command("get_exon", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("exon序列提取完成")
        elif command.return_code == -11:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("exon序列提取出错!")

    def get_intron(self):
        intron = self.option("gtf").prop["path"] + ".intron.bed"
        cmd = self.bedtools + "bedtools getfasta -fi %s -bed %s -fo %s -name" % (self.genome_fasta, intron, self.work_dir + "/" + self.ref_name + ".intron.fa")
        command = self.add_command("get_intron", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("intron序列提取完成")
        else:
            self.set_error("intron序列提取出错!")

    def set_output(self):
        intron_bed_path = self.option("gtf").prop["path"] + ".intron.bed"
        intron_bed = self.output_dir + "/" + self.ref_name + ".intron.bed"
        if os.path.exists(intron_bed):
            os.remove(intron_bed)
        os.link(intron_bed_path, intron_bed)
        exon_bed_path = self.option("gtf").prop["path"] + ".bed"
        exon_bed = self.output_dir + "/" + self.ref_name + ".exon.bed"
        if os.path.exists(exon_bed):
            os.remove(exon_bed)
        os.link(exon_bed_path, exon_bed)
        intron_path = self.work_dir + "/" + self.ref_name + ".intron.fa"
        intron = self.output_dir + "/" + self.ref_name + ".intron.fa"
        if os.path.exists(intron):
            os.remove(intron)
        os.link(intron_path, intron)
        exon_path = self.work_dir + "/" + self.ref_name + ".exon.fa"
        exon = self.output_dir + "/" + self.ref_name + ".exon.fa"
        if os.path.exists(exon):
            os.remove(exon)
        os.link(exon_path, exon)
        self.option('intron_fa', intron)
        self.option('exon_fa', exon)
        self.option("intron", intron_bed)
        self.option("exon", exon_bed)

    def run(self):
        super(IntronExonTool, self).run()
        self.exon_bed()
        self.intron_bed()
        self.get_exon()
        self.get_intron()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "IntronExon_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.intron_exon",
            "instant": False,
            "options": dict(
                input_genome="/mnt/ilustre/users/sanger-test/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                gtf="/mnt/ilustre/users/sanger-test/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()