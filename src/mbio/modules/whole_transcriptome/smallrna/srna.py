# -*- coding: utf-8 -*-
# __author__ = 'caiping.shi'

import os
import re
import time
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import shutil
import pandas as pd


class SrnaModule(Module):
    """
    miRNA质控的module,用于miRNA质控，多个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    """
    def __init__(self, work_id):
        super(SrnaModule, self).__init__(work_id)
        options = [
            {"name": "category", "type": "string", "default": ""}, # 物种分类，Animal or Plant
            {"name": "species", "type": "string", "default": ""}, # 具体物种
            {"name": "reference", "type": "infile", "format": "small_rna.fasta"}, # 参考基因组fasta文件
            {"name": "clean_fa", "type": "infile", "format": "small_rna.fasta"},  # 质控后的序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"}, # 鉴定完已知miRNA的过滤文件
            {"name": "gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考注释文件
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},  # 参考注释文件
            {"name": "repeat_result", "type": "string", "default": ""},
            {"name": "arf", "type": "infile", "format": "small_rna.common"},# clean reads mapping到参考基因组的结果
            {"name": "repeat", "type": "bool", "default": False}, ## 是否需要重新跑repeatmasker, True跳过，False需要重新跑
            {"name": "repeat_fa", "type": "infile", "format": "small_rna.fasta"}, ## repeatmasker注释出来的fasta文件
            {"name": "repeat_gff", "type": "infile", "format": "gene_structure.gff3"}, ## repeatmasker注释出来的gff文件
            {"name": "total_mirna_count", "type": "outfile", "format": "small_rna.common"}, # 合并的mirna count结果表
            {"name": "total_mirna_norm", "type": "outfile", "format": "small_rna.common"}, # 合并的mirna norm结果表
            {"name": "list", "type": "string", "default": ""}, # 样本列表，用来决定样本展示顺序
            {"name": "qc_output", "type": "string", "default": ""}, # 质控结果目录
            {"name": "extract_intron_exon", "type": "bool", "default": True}, # 是否直接从reads_vs_genome.arf结果中提取mapping到内含子和外显子区间的reads
            {"name": "mismatch", "type": "int", "default": 0},  # 允许错误匹配
            {"name": "bowtie_rfam", "type": "bool", "default": True},  # 是否采用bowtie比对到rfam
            {"name": "bowtie_repeat", "type": "bool", "default": True},  # 是否采用bowtie比对到repeat
        ]
        self.add_option(options)
        self.step.add_steps("srna")

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if self.option("category").lower() not in ["animal", "plant"]:
            raise OptionError("物种分类只能为Animal/Plant")
        if self.option("species") == "null":
            raise OptionError("必须指定具体物种")
        if not self.option("clean_fa").is_set:
            raise OptionError("必须设置质控后的文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        with open (self.option("clean_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("质控后的序列文件格式不对")
        if not self.option("gtf").is_set:
            if not self.option("gff").is_set:
                raise OptionError("gff和gtf中必须有一个作为参数传入")

    def run_known_mirna(self):
        options = {
            "species": self.option("species"),
            "clean_fa": self.option("clean_fa"),
            "config": self.option("config"),
            "mismatch": self.option("mismatch"),
        }
        self.known_mirna.set_options(options)
        if not self.option("extract_intron_exon"):
            self.known_mirna.on('end', self.run_blast_rfam)
        else:
            self.known_mirna.on('end', self.run_intron_exon)
        self.known_mirna.on('end', self.set_output, "known_mirna")
        self.known_mirna.on('start', self.set_step, {'start': self.step.known_mirna})
        self.known_mirna.on('end', self.set_step, {'end': self.step.known_mirna})
        self.known_mirna.run()

    def run_bowtie_rfam(self):
        options = {
            "file_type": "fasta",
            "seq_type": "se",
            "software": "bowtie2",
        }
        if not self.option("extract_intron_exon"):
            options.update({"se": self.known_mirna.option("filter_fa").prop["path"]})
        else:
            options.update({"se": self.genome_stat.option("filter_fa").prop["path"]})
        self.logger.info(options)
        self.bowtie_rfam.set_options(options)
        self.bowtie_rfam.on('end', self.run_rfam_stat)
        self.bowtie_rfam.on('start', self.set_step, {'start': self.step.bowtie_rfam})
        self.bowtie_rfam.on('end', self.set_step, {'end': self.step.bowtie_rfam})
        self.bowtie_rfam.run()

    def run_blast_rfam(self):
        options = {
            "database": "rfam",
            "lines": 50000,
            "evalue": 1e-3,
        }
        if not self.option("extract_intron_exon"):
            options.update({"query": self.known_mirna.option("filter_fa")})
        else:
            options.update({"query": self.genome_stat.option("filter_fa")})
        self.logger.info(options)
        self.blast_rfam.set_options(options)
        self.blast_rfam.on('end', self.run_xml2table)
        self.blast_rfam.on('start', self.set_step, {'start': self.step.blast_rfam})
        self.blast_rfam.on('end', self.set_step, {'end': self.step.blast_rfam})
        self.blast_rfam.run()

    def run_xml2table(self):
        options = {
            "blastout": self.blast_rfam.option("outxml"),
        }
        self.xml2table.set_options(options)
        self.xml2table.on('end', self.set_output, "rfam_blast")
        self.xml2table.on('end', self.run_rfam_stat)
        self.xml2table.on('start', self.set_step, {'start': self.step.xml2table})
        self.xml2table.on('end', self.set_step, {'end': self.step.xml2table})
        self.xml2table.run()

    def run_rfam_stat(self):
        options = {
            "config": self.option("config"),
        }
        if not self.option("extract_intron_exon"):
            options.update({"query": self.known_mirna.option("filter_fa")})
        else:
            options.update({"query": self.genome_stat.option("filter_fa")})
        if self.option("bowtie_rfam") == True:
            options.update({"bowtie_sam": os.path.join(self.bowtie_rfam.output_dir, "bowtie_rfam.sam")})
        else:
            options.update({"blast_table": self.xml2table.option("blast_table")})
        self.rfam_stat.set_options(options)
        self.rfam_stat.on('end', self.set_output, "rfam_stat")
        if self.option("repeat"):
            if self.option("bowtie_repeat") == True:
                self.rfam_stat.on('end', self.run_bowtie_repeat)
            else:
                self.rfam_stat.on('end', self.run_blast_repeat)
        else:
            self.rfam_stat.on('end', self.run_repeatmasker)
        self.rfam_stat.on('start', self.set_step, {'start': self.step.rfam_stat})
        self.rfam_stat.on('end', self.set_step, {'end': self.step.rfam_stat})
        self.rfam_stat.run()

    def run_repeatmasker(self):
        options = {
            "input_genome": self.option("reference"),
            "repeat_result": self.option("repeat_result"),
        }
        self.repeatmasker.set_options(options)
        self.repeatmasker.on('end', self.set_output, "repeatmasker")
        if self.option("bowtie_repeat") == True:
            self.repeatmasker.on('end', self.run_bowtie_repeat)
        else:
            self.repeatmasker.on('end', self.run_blast_repeat)
        self.repeatmasker.on('start', self.set_step, {'start': self.step.repeatmasker})
        self.repeatmasker.on('end', self.set_step, {'end': self.step.repeatmasker})
        self.repeatmasker.run()

    def run_bowtie_repeat(self):
        options = {
            "file_type": "fasta",
            "seq_type": "se",
            "software": "bowtie2",
            "se": self.rfam_stat.option("filter_fa").prop["path"],
        }
        if self.option("repeat"):
            options.update({"repeat": self.option("repeat_fa").prop["path"]})
        else:
            options.update({"repeat": self.repeatmasker.option("repeat").prop["path"]})
        self.logger.info(options)
        self.bowtie_repeat.set_options(options)
        self.bowtie_repeat.on('end', self.run_repeat_stat)
        self.bowtie_repeat.on('start', self.set_step, {'start': self.step.bowtie_repeat})
        self.bowtie_repeat.on('end', self.set_step, {'end': self.step.bowtie_repeat})
        self.bowtie_repeat.run()

    def run_blast_repeat(self):
        if self.option("repeat"):
            options = {
                "query": self.rfam_stat.option("filter_fa"),
                "database": "repeat",
                "reference": self.option("repeat_fa"),
            }
        else:
            options = {
                "query": self.rfam_stat.option("filter_fa"),
                "database": "repeat",
                "reference": self.repeatmasker.option("repeat"),
            }
        self.blast_repeat.set_options(options)
        self.blast_repeat.on('end', self.run_xml2table_1)
        self.blast_repeat.on('start', self.set_step, {'start': self.step.blast_repeat})
        self.blast_repeat.on('end', self.set_step, {'end': self.step.blast_repeat})
        self.blast_repeat.run()

    def run_xml2table_1(self):
        options = {
            "blastout": self.blast_repeat.option("outxml"),
        }
        self.xml2table_1.set_options(options)
        self.xml2table_1.on('end', self.set_output, "repeat_blast")
        self.xml2table_1.on('end', self.run_repeat_stat)
        self.xml2table_1.on('start', self.set_step, {'start': self.step.xml2table_1})
        self.xml2table_1.on('end', self.set_step, {'end': self.step.xml2table_1})
        self.xml2table_1.run()

    def run_repeat_stat(self):
        if not self.option("repeat"):
            options = {
                "query": self.rfam_stat.option("filter_fa"),
                "config": self.option("config"),
                "gff": self.repeatmasker.option("gff"),
            }
        else:
            options = {
                "query": self.rfam_stat.option("filter_fa"),
                "config": self.option("config"),
                "gff": self.option("repeat_gff"),
            }
        if self.option("bowtie_repeat") == True:
            options.update({"bowtie_sam": os.path.join(self.bowtie_repeat.output_dir, "bowtie_repeat.sam")})
        else:
            options.update({"blast_table": self.xml2table_1.option("blast_table")})
        self.repeat_stat.set_options(options)
        self.repeat_stat.on('end', self.set_output, "repeat_stat")
        if not self.option("extract_intron_exon"):
            self.repeat_stat.on('end', self.run_intron_exon)
        else:
            self.repeat_stat.on('end', self.run_novel_mirna)
        self.repeat_stat.on('start', self.set_step, {'start': self.step.repeat_stat})
        self.repeat_stat.on('end', self.set_step, {'end': self.step.repeat_stat})
        self.repeat_stat.run()

    def run_intron_exon(self):
        options = {
            "input_genome": self.option("reference"),
            "gtf": self.option("gtf"),
        }
        self.intron_exon.set_options(options)
        self.intron_exon.on('end', self.set_output, "intron_exon")
        if not self.option("extract_intron_exon"):
            self.intron_exon.on('end', self.run_blast_exon)
        else:
            self.intron_exon.on('end', self.run_genome_stat)
        self.intron_exon.on('start', self.set_step, {'start': self.step.intron_exon})
        self.intron_exon.on('end', self.set_step, {'end': self.step.intron_exon})
        self.intron_exon.run()

    def run_genome_stat(self):
        options = {
            "intron_bed": self.intron_exon.option("intron"),
            "exon_bed": self.intron_exon.option("exon"),
            "mapping_arf": self.option("arf"),
            "query": self.known_mirna.option("filter_fa"),
            "config": self.option("config"),
        }
        self.genome_stat.set_options(options)
        self.genome_stat.on('end', self.set_output, "genome_stat")
        if self.option("bowtie_rfam") == True:
            self.genome_stat.on('end', self.run_bowtie_rfam)
        else:
            self.genome_stat.on('end', self.run_blast_rfam)
        self.genome_stat.on('start', self.set_step, {'start': self.step.genome_stat})
        self.genome_stat.on('end', self.set_step, {'end': self.step.genome_stat})
        self.genome_stat.run()

    def run_blast_exon(self):
        options = {
            "query": self.repeat_stat.option("filter_fa"),
            "database": "exon",
            "reference": self.intron_exon.option("exon_fa"),
        }
        self.blast_exon.set_options(options)
        self.blast_exon.on('end', self.run_xml2table_2)
        self.blast_exon.on('start', self.set_step, {'start': self.step.blast_exon})
        self.blast_exon.on('end', self.set_step, {'end': self.step.blast_exon})
        self.blast_exon.run()

    def run_xml2table_2(self):
        options = {
            "blastout": self.blast_exon.option("outxml"),
        }
        self.xml2table_2.set_options(options)
        self.xml2table_2.on('end', self.set_output, "exon_blast")
        self.xml2table_2.on('end', self.run_exon_stat)
        self.xml2table_2.on('start', self.set_step, {'start': self.step.xml2table_2})
        self.xml2table_2.on('end', self.set_step, {'end': self.step.xml2table_2})
        self.xml2table_2.run()

    def run_exon_stat(self):
        options = {
            "blast_table": self.xml2table_2.option("blast_table"),
            "query": self.repeat_stat.option("filter_fa"),
            "config": self.option("config"),
        }
        self.exon_stat.set_options(options)
        self.exon_stat.on('end', self.set_output, "exon_stat")
        self.exon_stat.on('end', self.run_blast_intron)
        self.exon_stat.on('start', self.set_step, {'start': self.step.exon_stat})
        self.exon_stat.on('end', self.set_step, {'end': self.step.exon_stat})
        self.exon_stat.run()

    def run_blast_intron(self):
        options = {
            "query": self.exon_stat.option("filter_fa"),
            "database": "intron",
            "reference": self.intron_exon.option("intron_fa"),
        }
        self.blast_intron.set_options(options)
        self.blast_intron.on('end', self.run_xml2table_3)
        self.blast_intron.on('start', self.set_step, {'start': self.step.blast_intron})
        self.blast_intron.on('end', self.set_step, {'end': self.step.blast_intron})
        self.blast_intron.run()

    def run_xml2table_3(self):
        options = {
            "blastout": self.blast_intron.option("outxml"),
        }
        self.xml2table_3.set_options(options)
        self.xml2table_3.on('end', self.set_output, "intron_blast")
        self.xml2table_3.on('end', self.run_intron_stat)
        self.xml2table_3.on('start', self.set_step, {'start': self.step.xml2table_3})
        self.xml2table_3.on('end', self.set_step, {'end': self.step.xml2table_3})
        self.xml2table_3.run()

    def run_intron_stat(self):
        options = {
            "blast_table": self.xml2table_3.option("blast_table"),
            "query": self.exon_stat.option("filter_fa"),
            "config": self.option("config"),
        }
        self.intron_stat.set_options(options)
        self.intron_stat.on('end', self.set_output, "intron_stat")
        self.intron_stat.on('end', self.run_novel_mirna)
        self.intron_stat.on('start', self.set_step, {'start': self.step.intron_stat})
        self.intron_stat.on('end', self.set_step, {'end': self.step.intron_stat})
        self.intron_stat.run()

    def run_novel_mirna(self):
        options = {
            "species": self.option("species"),
            "category": self.option("category"),
            "input_genome": self.option("reference"),
            "config": self.option("config"),
            "mapping_arf": self.option("arf"),
            "mrd": self.known_mirna.option("mrd"),
        }
        if self.option("category").lower() == "plant":
            options.update({
                "known_miRNA_express": os.path.join(self.known_mirna.work_dir,"miRNAs_expressed_all_samples_total.csv"),
                "clean_fa": self.option("clean_fa")
            })
        if not self.option("extract_intron_exon"):
            options.update({"input_fa": self.intron_stat.option("filter_fa")})
        else:
            options.update({"input_fa": self.repeat_stat.option("filter_fa")})
        self.novel_mirna.set_options(options)
        self.novel_mirna.on('end', self.set_output, "novel_mirna")
        self.novel_mirna.on('end', self.run_srna_stat)
        self.novel_mirna.on('start', self.set_step, {'start': self.step.novel_mirna})
        self.novel_mirna.on('end', self.set_step, {'end': self.step.novel_mirna})
        self.novel_mirna.run()

    def run_srna_stat(self):
        options = {
            "known_mirna_count": self.known_mirna.option("miRNA_count").prop["path"],
            "known_mrd": self.known_mirna.output_dir + "/miRBase.mrd",
            "novel_mirna_count": self.novel_mirna.option("novel_count").prop["path"],
            "rfam_summary": self.rfam_stat.output_dir + "/rfam_summary.xls",
            "repeat_summary": self.repeat_stat.output_dir + "/repeat_summary.xls",
            "input_fa": self.option("clean_fa"),
            "config": self.option("config"),
            "category": self.option("category"),
            "list": self.option("list"),
            "qc_output": self.option("qc_output"),
        }
        if self.option("category").lower() == "animal":
            options.update({'novel_mrd': self.novel_mirna.output_dir + "/miRBase.mrd"})
        else:
            options.update({'novel_mrd': self.novel_mirna.work_dir + "/mireap-Nov.aln"})
        if not self.option("extract_intron_exon"):
            options.update({'exon_summary': self.exon_stat.output_dir + "/exon_summary.xls"})
            options.update({'intron_summary': self.intron_stat.output_dir + "/intron_summary.xls"})
        else:
            options.update({'exon_summary': self.genome_stat.output_dir + "/exon_summary.xls"})
            options.update({'intron_summary': self.genome_stat.output_dir + "/intron_summary.xls"})
        self.srna_stat.set_options(options)
        self.srna_stat.on('end', self.set_output, "srna_stat")
        self.srna_stat.on('end', self.end)
        self.srna_stat.on('start', self.set_step, {'start': self.step.srna_stat})
        self.srna_stat.on('end', self.set_step, {'end': self.step.srna_stat})
        self.srna_stat.run()

    def set_output(self, event):
        # pass
        obj = event["bind_object"]
        if event['data'] == 'known_mirna':
            self.move2outputdir(obj.output_dir, 'known_mirna')
        if event['data'] == 'rfam_blast':
            self.move2outputdir(obj.output_dir, 'rfam_blast')
        if event['data'] == 'rfam_stat':
            self.move2outputdir(obj.output_dir, 'rfam_stat')
        if event['data'] == 'repeatmasker':
            self.move2outputdir(obj.output_dir, 'repeatmasker')
        if event['data'] == 'repeat_blast':
            self.move2outputdir(obj.output_dir, 'repeat_blast')
        if event['data'] == 'repeat_stat':
            self.move2outputdir(obj.output_dir, 'repeat_stat')
        if event['data'] == 'intron_exon':
            self.move2outputdir(obj.output_dir, 'intron_exon')
        if event['data'] == 'exon_blast':
            self.move2outputdir(obj.output_dir, 'exon_blast')
        if event['data'] == 'exon_stat':
            self.move2outputdir(obj.output_dir, 'exon_stat')
        if event['data'] == 'intron_blast':
            self.move2outputdir(obj.output_dir, 'intron_blast')
        if event['data'] == 'intron_stat':
            self.move2outputdir(obj.output_dir, 'intron_stat')
        if event['data'] == 'srna_stat':
            self.move2outputdir(obj.output_dir, 'srna_stat')
        if event['data'] == 'novel_mirna':
            self.move2outputdir(obj.output_dir, 'novel_mirna')
        if event['data'] == 'bowtie_rfam':
            self.move2outputdir(obj.output_dir, 'bowtie_rfam')
        if event['data'] == 'bowtie_repeat':
            self.move2outputdir(obj.output_dir, 'bowtie_repeat')

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)

    def run(self):
        self.known_mirna = self.add_tool("small_rna.srna.known_mirna")
        if self.option("bowtie_rfam") == True:
            self.bowtie_rfam = self.add_tool("small_rna.srna.bowtie_rfam")
        else:
            self.blast_rfam = self.add_module("small_rna.srna.blast")
            self.xml2table = self.add_tool("small_rna.srna.xml2table")
        self.rfam_stat = self.add_tool("small_rna.srna.rfam_stat")
        self.repeatmasker = self.add_module("small_rna.srna.repeatmasker")
        if self.option("bowtie_repeat") == True:
            self.bowtie_repeat = self.add_tool("small_rna.srna.bowtie_repeat")
        else:
            self.blast_repeat = self.add_module("small_rna.srna.blast")
            self.xml2table_1 = self.add_tool("small_rna.srna.xml2table")
        self.repeat_stat = self.add_tool("small_rna.srna.repeat_stat")
        self.intron_exon = self.add_tool("small_rna.srna.intron_exon")
        if not self.option("extract_intron_exon"):
            self.blast_exon = self.add_module("small_rna.srna.blast")
            self.xml2table_2 = self.add_tool("small_rna.srna.xml2table")
            self.exon_stat = self.add_tool("small_rna.srna.exon_stat")
            self.blast_intron = self.add_module("small_rna.srna.blast")
            self.xml2table_3 = self.add_tool("small_rna.srna.xml2table")
            self.intron_stat = self.add_tool("small_rna.srna.intron_stat")
        else:
            self.genome_stat = self.add_tool("small_rna.srna.genome_stat")
        self.novel_mirna = self.add_tool("small_rna.srna.novel_mirna")
        self.srna_stat = self.add_tool("small_rna.srna.srna_stat")
        if not self.option("extract_intron_exon"):
            self.step.add_steps('known_mirna', 'blast_rfam', "xml2table", "rfam_stat",
                                "repeatmasker", "blast_repeat", "xml2table_1", "repeat_stat",
                                "intron_exon", "blast_exon", "xml2table_2", "exon_stat",
                                "blast_intron", "xml2table_3", "intron_stat", "novel_mirna", "srna_stat", "bowtie_rfam", "bowtie_repeat")
        else:
            self.step.add_steps('known_mirna', 'blast_rfam', "xml2table", "rfam_stat",
                                "repeatmasker", "blast_repeat", "xml2table_1", "repeat_stat",
                                "intron_exon", "genome_stat", "novel_mirna", "srna_stat", "bowtie_rfam", "bowtie_repeat")
        self.run_known_mirna()
        super(SrnaModule, self).run()

    def end(self):
        known_mirna_count = os.path.join(self.known_mirna.output_dir, "known_mirna_count.xls")
        novel_mirna_count = os.path.join(self.novel_mirna.output_dir, "novel_mirna_count.xls")
        if os.path.exists(os.path.join(self.output_dir, "known_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "known_mirna_count.xls"))
        os.link(known_mirna_count, os.path.join(self.output_dir, "known_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(novel_mirna_count, os.path.join(self.output_dir, "novel_mirna_count.xls"))
        total_mirna_norm = os.path.join(self.output_dir, "total_mirna_norm.xls")
        total_mirna_count = os.path.join(self.output_dir, "total_mirna_count.xls")
        known_mirna_norm = os.path.join(self.output_dir, "known_mirna_norm.xls")
        novel_mirna_norm = os.path.join(self.output_dir, "novel_mirna_norm.xls")
        known_mirnas = list()
        novel_mirnas = list()
        with open(known_mirna_count, "r") as f1, open(novel_mirna_count, "r") as f2, open(total_mirna_count, "w") as w1:
            head_f2 = f2.readline()
            head_f1 = f1.readline()
            w1.write(head_f1)
            for line in f1:
                known_mirnas.append(line.strip().split("\t")[0])
                w1.write(line)
            for line in f2:
                novel_mirnas.append(line.strip().split("\t")[0])
                w1.write(line)
        total_mirna_count_pd = pd.read_table(total_mirna_count, header=0, index_col=0)
        samples = total_mirna_count_pd.columns
        for sample in samples:
            total_mirna_count_pd[sample] = (
                    total_mirna_count_pd[sample] / total_mirna_count_pd[sample].sum() * 1000000).round(4)
        known_mirna_norm_pd = total_mirna_count_pd.ix[known_mirnas]
        novel_mirna_norm_pd = total_mirna_count_pd.ix[novel_mirnas]
        total_mirna_count_pd.to_csv(total_mirna_norm, sep="\t")
        known_mirna_norm_pd.to_csv(known_mirna_norm, sep="\t")
        novel_mirna_norm_pd.to_csv(novel_mirna_norm, sep="\t")
        os.remove(os.path.join(self.output_dir, "known_mirna", "known_mirna_norm.xls"))
        os.remove(os.path.join(self.output_dir, "novel_mirna", "novel_mirna_norm.xls"))
        shutil.copy(known_mirna_norm, os.path.join(self.output_dir, "known_mirna", "known_mirna_norm.xls"))
        shutil.copy(novel_mirna_norm, os.path.join(self.output_dir, "novel_mirna", "novel_mirna_norm.xls"))
        self.option("total_mirna_count", self.output_dir + "/total_mirna_count.xls")
        self.option("total_mirna_norm", self.output_dir + "/total_mirna_norm.xls")
        super(SrnaModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "srna" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "small_rna.srna.srna",
            "instant": False,
            "options": dict(
                category="Animal",
                species="mmu",
                clean_fa="/mnt/lustre/users/sanger/workspace/20190524/Smallrna_sanger_181227/FastaUniq/uniq.fasta",
                config="/mnt/lustre/users/sanger/workspace/20190524/Smallrna_sanger_181227/MirnaQc/qc_file.config",
                reference="/mnt/lustre/users/sanger/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa",
                arf="/mnt/lustre/users/sanger/workspace/20190524/Smallrna_sanger_181227/MapperAndStat__1/output/reads_vs_genome.arf",
                gtf="/mnt/lustre/users/sanger/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/gtf/Mus_musculus.GRCm38.96.gtf",
                repeat=True,
                extract_intron_exon=True,
                repeat_fa="/mnt/lustre/users/sanger/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/Annotation_v2/repeatmasker/repeatmasker_merge.repeat.fa",
                repeat_gff="/mnt/lustre/users/sanger/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/Annotation_v2/repeatmasker/repeatmasker_merge.SSR.gff",
                list="/mnt/lustre/users/sanger/workspace/20190524/Smallrna_sanger_181227/MirnaQc/output/clean_data/list.txt",
                qc_output="/mnt/lustre/users/sanger/workspace/20190524/Smallrna_sanger_181227/MirnaQc/output"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()