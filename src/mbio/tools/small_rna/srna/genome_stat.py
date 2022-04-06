# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import ConfigParser
import os
import subprocess
import time
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from skbio.parse.sequences import parse_fasta


class GenomeStatAgent(Agent):
    """
    repeat比对结果统计
    """

    def __init__(self, parent):
        super(GenomeStatAgent, self).__init__(parent)
        options = [
            {"name": "intron_bed", "type": "infile", "format": "gene_structure.bed"},  # intron bed文件
            {"name": "exon_bed", "type": "infile", "format": "gene_structure.bed"},  # exon bed文件
            {"name": "mapping_arf", "type": "infile", "format": "small_rna.common"},  # clean reads mapping到参考基因组的结果
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完已知miRNA的过滤文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
        ]
        self.add_option(options)
        self.step.add_steps("genome_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.genome_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.genome_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("intron_bed").is_set:
            raise OptionError("必须提供intron bed文件")
        if not self.option("exon_bed").is_set:
            raise OptionError("必须提供exon bed文件")
        if not self.option("mapping_arf").is_set:
            raise OptionError("必须提供mapping arf文件")
        if not self.option("query").is_set:
            raise OptionError("必须提供输入FASTA文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
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
        super(GenomeStatAgent, self).end()


class GenomeStatTool(Tool):
    def __init__(self, config):
        super(GenomeStatTool, self).__init__(config)
        self.python = "miniconda2/bin/python"
        self.intersectBed = self.config.SOFTWARE_DIR + "/bioinfo/seq/bedtools-2.25.0/bin/intersectBed"
        intersect_bed_v2170 = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/miniconda3/bin/intersectBed')
        if os.path.isfile(intersect_bed_v2170):
            self.intersectBed = intersect_bed_v2170

    def run(self):
        super(GenomeStatTool, self).run()
        self.parse_arf()
        self.intergenic_bed()
        self.intergenic_stat()
        self.intron_bed()
        self.intron_stat()
        self.exon_bed()
        self.exon_stat()
        self.end()

    def parse_arf(self):
        self.logger.info("将mapping arf文件转换成bed文件")
        start = time.time()
        mapping_arf = self.option("mapping_arf").prop["path"]
        mapping_bed = self.output_dir + "/genome_mapping.bed"
        input_fa = self.option("query").prop["path"]
        self.in_fa = {}
        for seq_id, seq_sequence in parse_fasta(input_fa):
            self.in_fa[seq_id] = seq_sequence
        with open(mapping_arf, "r") as f, open(mapping_bed, "w") as w1:
            for line in f:
                items = line.strip().split("\t")
                if items[0] in self.in_fa:
                    w1.write(
                        items[5] + "\t" + items[7] + "\t" + items[8] + "\t" + items[0] + "\t0\t" + items[10] + "\n")
        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

    def intergenic_bed(self):
        self.logger.info("开始计算mapping上基因间区区域的reads")
        start = time.time()
        mapping_bed = self.output_dir + "/genome_mapping.bed"
        exon_bed = self.option("exon_bed").prop["path"]

        intergenic_mapping_bed = self.output_dir + "/intergenic_mapping.bed"
        cmd = """%s -a %s -b %s -v> %s""" % (self.intersectBed, mapping_bed, exon_bed, intergenic_mapping_bed)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成intergenic_mapping_bed文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成intergenic_mapping_bed文件失败')
            self.set_error("生成intergenic_mapping_bed文件失败")
        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

    def intergenic_stat(self):
        self.logger.info("开始对序列比对intergenic区域进行统计")
        start = time.time()
        intergenic_mapping_bed = self.output_dir + "/intergenic_mapping.bed"
        intergenic_mapping_out = self.output_dir + "/intergenic_mapping.xls"
        intergenic_seq_file = self.work_dir + "/intergenic_mapping_seq.txt"
        with open(intergenic_mapping_bed, "r") as f1, open(intergenic_mapping_out, "w") as w1:
            intergenic_seq_id = {}
            self.intergenic_fa = {}
            samples = {}
            intergenic_stat = {}
            head = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                seq_id = items[3]
                if seq_id in self.in_fa:
                    if seq_id not in self.intergenic_fa:
                        self.intergenic_fa[seq_id] = self.in_fa[seq_id]
                    self.in_fa.pop(seq_id)
                    sample = seq_id.split("_")[0]
                    num = int(seq_id.split("_x")[1])
                    if sample not in samples:
                        samples[sample] = 1
                    if seq_id not in intergenic_seq_id:
                        intergenic_seq_id[seq_id] = 1
                        if intergenic_stat.has_key(sample):
                            intergenic_stat[sample] += num
                        else:
                            intergenic_stat[sample] = num
                    else:
                        pass
                else:
                    pass
            w1.write("Samples\tIntergenic\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not intergenic_stat.has_key(key):
                    intergenic_stat[key] = 0
                w1.write("{}\t{}\n".format(sample_name, intergenic_stat[key]))
            with open(intergenic_seq_file, "w") as w:
                for key in intergenic_seq_id:
                    w.write("{}\n".format(key))
        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

    def intron_bed(self):
        self.logger.info("开始计算mapping上intron区域的reads")
        start = time.time()
        mapping_bed = self.output_dir + "/genome_mapping.bed"
        intron_bed = self.option("intron_bed").prop["path"]
        if os.path.getsize(intron_bed) == 0:
            intron_mapping_bed = self.output_dir + "/intron_mapping.bed"
            cmd = "touch %s" % (intron_mapping_bed)
            subprocess.check_output(cmd, shell=True)
        else:
            intron_mapping_bed = self.output_dir + "/intron_mapping.bed"
            cmd = """%s -a %s -b %s -wo> %s""" % (self.intersectBed, mapping_bed, intron_bed, intron_mapping_bed)
            self.logger.info(cmd)
            try:
                subprocess.check_output(cmd, shell=True)
                self.logger.info('生成intron_mapping_bed文件成功')
            except subprocess.CalledProcessError:
                self.logger.info('生成intron_mapping_bed文件失败')
                self.set_error("生成intron_mapping_bed文件失败")
            end = time.time()
            self.logger.info("运行时间为{}s".format(end - start))

    def exon_bed(self):
        self.logger.info("开始计算mapping上exon区域的reads")
        start = time.time()
        mapping_bed = self.output_dir + "/genome_mapping.bed"
        exon_bed = self.option("exon_bed").prop["path"]
        exon_mapping_bed = self.output_dir + "/exon_mapping.bed"
        cmd = """%s -a %s -b %s -wo> %s""" % (self.intersectBed, mapping_bed, exon_bed, exon_mapping_bed)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成exon_mapping_bed文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成exon_mapping_bed文件失败')
            self.set_error("生成exon_mapping_bed文件失败")
        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

    def intron_stat(self):
        self.logger.info("开始对序列比对intron区域进行统计")
        start = time.time()
        intron_mapping_bed = self.output_dir + "/intron_mapping.bed"
        intron_summary_out = self.output_dir + "/intron_summary.xls"
        intron_seq_file = self.work_dir + "/intron_mapping_seq.txt"
        with open(intron_mapping_bed, "r") as f1, open(intron_summary_out, "w") as w1:
            intron_seq_id = {}
            samples = {}
            intron_stat = {}
            head = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                seq_id = items[3]
                overlap_ratio = float(int(items[-1]) / (int(items[2]) - int(items[1])))
                if overlap_ratio > 0.9 and items[5] == items[11]:
                    if seq_id in self.in_fa:
                        self.in_fa.pop(seq_id)
                        sample = seq_id.split("_")[0]
                        num = int(seq_id.split("_x")[1])
                        if sample not in samples:
                            samples[sample] = 1
                        if seq_id not in intron_seq_id:
                            intron_seq_id[seq_id] = 1
                            if intron_stat.has_key(sample):
                                intron_stat[sample] += num
                            else:
                                intron_stat[sample] = num
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            w1.write("Samples\tIntron\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not intron_stat.has_key(key):
                    intron_stat[key] = 0
                w1.write("{}\t{}\n".format(sample_name, intron_stat[key]))
            with open(intron_seq_file, "w") as w:
                for key in intron_seq_id:
                    w.write("{}\n".format(key))

        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

    def exon_stat(self):
        self.logger.info("开始对序列比对exon区域进行统计")
        start = time.time()
        exon_mapping_bed = self.output_dir + "/exon_mapping.bed"
        exon_summary_out = self.output_dir + "/exon_summary.xls"
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        exon_seq_file = self.work_dir + "/exon_mapping_seq.txt"
        with open(exon_mapping_bed, "r") as f1, open(exon_summary_out, "w") as w1:
            exon_seq_id = {}
            samples = {}
            exon_stat = {}
            head = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                seq_id = items[3]
                overlap_ratio = float(int(items[-1]) / (int(items[2]) - int(items[1])))
                if overlap_ratio > 0.9 and items[5] == items[11]:
                    if seq_id in self.in_fa:
                        self.in_fa.pop(seq_id)
                        sample = seq_id.split("_")[0]
                        num = int(seq_id.split("_x")[1])
                        if sample not in samples:
                            samples[sample] = 1
                        if seq_id not in exon_seq_id:
                            exon_seq_id[seq_id] = 1
                            if exon_stat.has_key(sample):
                                exon_stat[sample] += num
                            else:
                                exon_stat[sample] = num
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            w1.write("Samples\tExon\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not exon_stat.has_key(key):
                    exon_stat[key] = 0
                w1.write("{}\t{}\n".format(sample_name, exon_stat[key]))

        with open(filter_fa, "w") as w:
            for key in self.in_fa:
                w.write('>%s\n%s\n' % (key, self.in_fa[key]))
            for key in self.intergenic_fa:
                w.write('>%s\n%s\n' % (key, self.intergenic_fa[key]))
        with open(exon_seq_file, "w") as w:
            for key in exon_seq_id:
                w.write("{}\n".format(key))
        self.option('filter_fa', filter_fa)
        end = time.time()
        self.logger.info("运行时间为{}s".format(end - start))

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
            "id": "GenomeStat_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.genome_stat",
            "instant": False,
            "options": dict(
                intron_bed="/mnt/lustre/users/sanger/workspace/20190529/Smallrna_i-sanger_177340/Srna/IntronExon/output/Populus_trichocarpa.intron.bed",
                exon_bed="/mnt/lustre/users/sanger/workspace/20190529/Smallrna_i-sanger_177340/Srna/IntronExon/output/Populus_trichocarpa.exon.bed",
                mapping_arf="/mnt/lustre/users/sanger/workspace/20190529/Smallrna_i-sanger_177340/MapperAndStat/reads_vs_genome.arf",
                config="/mnt/lustre/users/sanger/workspace/20190529/Smallrna_i-sanger_177340/MirnaQc/output/clean_data/qc_file.config",
                query="/mnt/lustre/users/sanger/workspace/20190529/Smallrna_i-sanger_177340/Srna/RepeatStat/output/filtered.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
