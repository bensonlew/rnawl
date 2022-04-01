# -*- coding:utf-8 -*-
# __author__ = '封一统'
"""原核转录组一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import subprocess
import datetime
import shutil
import re
import time
import gevent
import functools
from biocluster.config import Config
import pandas as pd
from mbio.packages.prok_rna.copy_file import CopyFile
import unittest
from biocluster.wpm.client import *


# 定义用于统计导表时间的装饰器
def time_count(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run ' + func_name + ' at ' + start_time)
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End ' + func_name + ' at ' + end_time)
        print("{}函数执行时间约为{}s".format(func.__name__, end - start))
    return wrapper


class ProkrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(ProkrnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "strand_specific", "type": "bool", "default": True}, # 当为PE测序时，是否为链特异性
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "align_species", "type": "string", "default": "is_ncbi"},  # 选择参考物种是来着NCBI或用户自己上传
            {"name": "species_name", "type":"string", "default": ""}, # 客户研究对象在NCBI里的命名
            {"name": "genome_id", "type":"string", "default": ""}, # 客户研究对象的参考基因组等在NCBI里的accession号
            {"name": "genome_db", "type": "infile", "format": "prok_rna.fasta"},  # 客户自行上传的参考基因组
            {"name": "gff_or_gtf", "type": "string", "default": "gff"},  # 客户自行上传参考基因组时需要选择是上传gtf或者gff文件
            {"name": "gff_or_gtf_file", "type": "infile", "format": "prok_rna.common"},  # 客户自行上传参考基因组时上传的gtf或者gff文件
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"}, # 对照表
            {"name": "mapping_stop", "type": "string", "default": ''}, # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_sample_percent", "type": "string", "default": ""}, # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_ratio", "type": "string", "default": ""}, # 设定比对率小于多少时判定样本不达标
            {"name": "rrna_stop", "type": "string", "default": ''}, # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_sample_percent", "type": "string", "default": ""}, # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_ratio", "type": "string", "default": ""}, # 设定rRNA比率大于多少时判定样本不达标
            {"name": "map_assess_method", "type": "string", "default": "saturation,distribution,coverage,chr_stat"},
            ## 高级参数设置
            # 转录组功能注释
            {"name": "nr_evalue", "type": "string", "default": "1e-5"},
            {"name": "kegg_evalue", "type": "string", "default": "1e-5"},
            {"name": "cog_evalue", "type": "string", "default": "1e-5"},
            {"name": "swissprot_evalue", "type": "string", "default": "1e-5"},
            {"name": "pfam_evalue", "type": "string", "default": "1e-5"},
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},  # NR比对e值
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},  # KEGG注释使用的e值
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-5},  # Swissprot比对使用的e值
            {"name": "cog_blast_evalue", "type": "float", "default": 1e-5},  # COG比对使用的e值
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},  # PFAM比对使用的e值
            # 表达量分析
            {"name": "express_method", "type": "string", "default": "rsem"},# 表达量分析手段: Salmon, Kallisto, RSEM
            {"name": "exp_way", "type": "string", "default": "fpkm"}, # fpkm or tpm
            # 表达差异分析
            {"name": "diff_method", "type": "string", "default": "DESeq2"},# 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  #选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  #Bonferroni,Holm,BH,BY
        ]
        #获取输出目录
        self.workflow_output_tmp = self._sheet.output
        self.logger.info(self.workflow_output_tmp)
        try:
            if re.match(r'tsanger:',self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
                self.workflow_output2 = self.workflow_output_tmp.replace('tsanger:','')
            elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp
                self.workflow_output2 = self.workflow_output_tmp
            else:
                self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
                self.workflow_output2 = self.workflow_output_tmp.replace('sanger:','')
        except:
            self.workflow_output = 'rere...'
            self.workflow_output2 = '/mnt/...'
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())

        #获取数据库中该物种相关文件路径
        if self.option('align_species') == 'is_ncbi':
            self.db = Config().get_mongo_client(mtype="prok_rna", dydb_forbid=True)[Config().get_mongo_dbname("prok_rna", dydb_forbid=True)]
            col = self.db["sg_genome_db"]
            genome_info = col.find_one({"organism_name" : self.option("species_name"), "genome_id" : self.option("genome_id")})
            if not genome_info:
                self.set_error("数据库中不存在该物种注释信息，程序退出", code = "15000101")
        self.db_path = self.config.SOFTWARE_DIR + "/database/Transeq_DB_v1/"
        if self.option('align_species') == 'is_ncbi':
            if not os.path.exists(os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna')):
                self.set_error("服务器内不存在您输入的NCBI accession号的参考基因组文件，请检查", code = "15000102")
            if not os.path.exists(os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff')):
                self.set_error("服务器内不存在您输入的NCBI accession号的参考基因组的gff文件，请检查", code = "15000103")
        else:
            if not self.option('genome_db').is_set or not os.path.exists(self.option('genome_db').prop['path']):
                self.set_error("没有找到您上传的参考基因组", code = "15000104")
            if not self.option('gff_or_gtf_file').is_set or not os.path.exists(self.option('gff_or_gtf_file').prop['path']):
                self.set_error("没有找到您上传的gff或gtf文件", code = "15000105")

        self.ref_genome = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option(
            'align_species') == 'is_ncbi' else self.option('genome_db').prop['path']

        #添加tool/module
        self.rock_index = self.add_tool("prok_rna.rockhopper_index")
        self.rockhopper = self.add_tool("prok_rna.rockhopper")
        self.filecheck = self.add_tool("prok_rna.file_check")
        self.qc = self.add_module("prok_rna.hiseq_qc")
        self.qc_stat_before = self.add_module("prok_rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("prok_rna.hiseq_reads_stat")
        self.mapping = self.add_module("prok_rna.rnaseq_mapping")
        self.annot_filter = self.add_module("prok_rna.annot_filter")
        self.annot_class = self.add_module("prok_rna.annot_class")
        self.annot_mapdb = self.add_module("prok_rna.annot_mapdb")
        self.annot_orfpfam = self.add_module("prok_rna.annot_orfpfam")
        self.map_qc = self.add_module("prok_rna.map_assessment")
        self.snp = self.add_module("prok_rna.sam_rna")
        self.express = self.add_module("prok_rna.quant")
        self.exp_pca = self.add_tool("prok_rna.exp_pca")
        self.exp_corr = self.add_tool("prok_rna.exp_corr")
        self.exp_venn = self.add_tool('prok_rna.exp_venn')
        self.diffexpress = self.add_tool("prok_rna.diffexp")
        self.srna = self.add_module("prok_rna.srna")
        self.promote = self.add_tool('prok_rna.promote')

        #判断流程结束这次采用计数的方法，对于特定的tool，set_output之后这个数会加一，等这个数为7时，流程结束
        # self.end_number = 0
        # if self.option('strand_specific'):
        #     self.to_end = 7
        # else:
        #     self.to_end = 6

        self.final_tools = [self.qc_stat_before, self.map_qc, self.snp, self.exp_pca, self.exp_corr, self.exp_venn, self.promote, self.srna, self.annot_class]
        # self.final_tools.remove(self.snp)

        #一个判断是否通过了检查的参数，通过了为True,不通过为False
        self.pass_judge = True

        # 添加step，显示在页面进度条
        all_steps = ["rock_index", "filecheck", "rna_qc", "mapping", "express", "diffexpress", "snp_rna",
                     "map_qc", "annot_mapdb", "annot_orfpfam", "annot_filter",  "annot_class", "exp_pca",
                     "exp_corr", "exp_venn", "qc_stat_before", "qc_stat_after", "rockhopper", "srna", "promote"
                     ]
        if self.option("group_table").prop["sample_number"] < 3:
            all_steps.remove('exp_pca')
            for step in all_steps:
                self.step.add_steps(step)
        else:
            for step in all_steps:
                self.step.add_steps(step)

    def check_options(self):
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE", code = "15000106")
        try:
            nr_evalue = float(self.option("nr_evalue"))
            cog_evalue = float(self.option("cog_evalue"))
            kegg_evalue = float(self.option("kegg_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范", code = "15000107")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("cog_blast_evalue", cog_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围", code = "15000108")
        if not self.option("cog_blast_evalue") > 0 and not self.option("cog_blast_evalue") < 1:
            raise OptionError("Cog比对的E值超出范围", code = "15000109")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围", code = "15000110")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围", code = "15000111")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围", code = "15000112")
        if self.option('align_species') != 'is_ncbi':
            if self.option('gff_or_gtf_file').prop['path'].endswith('gtf') and self.option('gff_or_gtf') == 'gff':
                raise OptionError("说好的你给我传gff，结果你给我gtf", code = "15000113")
            if self.option('gff_or_gtf_file').prop['path'].endswith('gff') and self.option('gff_or_gtf') == 'gtf':
                raise OptionError("说好的你给我传gtf，结果你给我gff", code = "15000114")
        if self.option('mapping_sample_percent') or self.option('rrna_sample_percent'):
            try:
                mapping_sample_percent = float(self.option("mapping_sample_percent"))
                mapping_ratio = float(self.option("mapping_ratio"))
                rrna_sample_percent = float(self.option("rrna_sample_percent"))
                rrna_ratio = float(self.option("rrna_ratio"))
            except:
                raise OptionError("判断终止那边不要传入非数字形式", code = "15000115")
            else:
                if mapping_ratio and not mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整", code = "15000116")
                if not mapping_ratio and mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整", code = "15000117")
                if rrna_sample_percent and not rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整", code = "15000118")
                if not rrna_sample_percent and rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整", code = "15000119")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        """
        ref-rna workflow run方法
        :return:
        """
        self.rock_index.on('end', self.run_filecheck)
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc_stat_after.on('end', self.run_mapping)
        self.mapping.on('end', self.run_map_assess)
        self.mapping.on('end', self.run_snp)
        self.mapping.on('end', self.run_annot_mapdb)
        self.mapping.on('end', self.run_annot_orfpfam)
        self.mapping.on('end', self.run_promote)
        self.on_rely([self.annot_mapdb,self.annot_orfpfam], self.run_annot_filter)
        self.annot_filter.on('end', self.run_annot_class)
        if self.option("strand_specific"):
            self.mapping.on('end', self.run_rockhopper)
            self.mapping.on('end', self.run_srna)
            self.rockhopper.on('end', self.run_express)
        else:
            self.mapping.on('end', self.run_express)
            self.final_tools.remove(self.srna)
        self.express.on('end', self.run_diffexpress)
        self.diffexpress.on('end', self.run_exp_corr)
        self.diffexpress.on('end', self.run_exp_pca)
        self.diffexpress.on('end', self.run_exp_venn)
        self.on_rely(self.final_tools, self.end)
        self.run_rock_index()
        super(ProkrnaWorkflow, self).run()

    def run_rock_index(self):
        if self.option("align_species") == 'is_ncbi':
            opts = {
                "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                "input_file": os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff'),
                "type": "gff"
            }
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": self.option('gff_or_gtf_file').prop['path'],
                "type": self.option('gff_or_gtf')
            }
        self.rock_index.set_options(opts)
        self.rock_index.on('end', self.set_output, 'rock_index')
        self.rock_index.on('start', self.set_step, {'start': self.step.rock_index})
        self.rock_index.on('end', self.set_step, {'end': self.step.rock_index})
        self.rock_index.run()
        pass

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            "in_gtf": self.rock_index.option('gtf_exon'),
        }
        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        if self.option("control_file").is_set:
            opts.update({'control_file': self.option('control_file')})
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        self.qc.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'quality_score_system': self.option('quality_score_system')
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.rna_qc})
        self.qc.on('end', self.set_step, {'end': self.step.rna_qc})
        self.qc.run()

    def run_qc_stat(self, event):
        if self.option('quality_score_system').lower() == "phred+33" or self.option('quality_score_system').lower() == "phred 33":
            quality = 33
        else:
            quality = 64
        if event['data']:
            self.qc_stat_after.set_options({
                'fastq_dir': self.qc.option('sickle_dir'),
                'fq_type': self.option('fq_type'),
                'quality': quality,
                'rfam':True,
                'dup': True
            })
        else:
            self.qc_stat_before.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality': quality,
            })
        if event['data']:
            self.qc_stat_after.on('end', self.set_output, 'qc_stat_after')
            self.qc_stat_after.on("start", self.set_step, {"start": self.step.qc_stat_after})
            self.qc_stat_after.on("end", self.set_step, {"end": self.step.qc_stat_after})
            self.qc_stat_after.run()
        else:
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
            self.qc_stat_before.on("start", self.set_step, {"start": self.step.qc_stat_before})
            self.qc_stat_before.on("end", self.set_step, {"end": self.step.qc_stat_before})
            self.qc_stat_before.run()

    def run_mapping(self):
        opts = {
            "ref_genome": self.ref_genome,
            "seq_method": self.option("fq_type"),   # PE or SE
            "fastq_dir": self.qc.option("sickle_dir"),
            "ref_gtf": self.rock_index.option("gtf"),
        }
        if self.option("strand_specific"):
            opts.update(
                {
                    "strand_specific":True,
                 }
            )
        self.mapping.set_options(opts)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_rockhopper(self):
        self.logger.info("开始运行rockhopper")
        if self.option("align_species") == 'is_ncbi':
            opts = {
                "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                "input_file": os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff'),
                "type": "gff"
            }
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": self.option('gff_or_gtf_file').prop['path'],
                "type": self.option('gff_or_gtf')
            }
        opts.update({"group_list": self.option('group_table').prop['path'], "trimPairFq": self.qc.option('fq_list').prop['path']})
        self.rockhopper.set_options(opts)
        self.rockhopper.on("end", self.set_output, "rockhopper")
        self.rockhopper.on('start', self.set_step, {'start': self.step.rockhopper})
        self.rockhopper.on('end', self.set_step, {'end': self.step.rockhopper})
        self.rockhopper.run()

    def run_annot_mapdb(self, event):
        self.logger.info("开始运行diamond注释")
        opts = {
            "query" : self.rock_index.option("query"),
            "method": "diamond",
            "nr_db": "bacteria",
            "lines": 5000
        }
        self.annot_mapdb.set_options(opts)
        self.annot_mapdb.on("end", self.set_output, "annot_mapdb")
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.logger.info("开始运行pfam注释")
        opts = {
            "lines": 5000,
            "pep": self.rock_index.option("querypep"),
        }
        self.annot_orfpfam.set_options(opts)
        self.annot_orfpfam.on("end", self.set_output, "annot_orfpfam")
        self.annot_orfpfam.on('start', self.set_step, {'start': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_step, {'end': self.step.annot_orfpfam})
        self.annot_orfpfam.run()

    def run_annot_filter(self):
        options = {
            "blast_nr_xml" : self.annot_mapdb.output_dir + "/nr/blast.xml",
            "blast_eggnog_xml" : self.annot_mapdb.output_dir + "/eggnog/blast.xml",
            "blast_kegg_xml": self.annot_mapdb.output_dir + "/kegg/blast.xml",
            "blast_swissprot_xml" : self.annot_mapdb.output_dir + "/swissprot/blast.xml",
            "pfam_domain" : self.annot_orfpfam.output_dir + "/pfam_domain",
            "blast2go_annot" : self.annot_mapdb.output_dir + "/GO/go_annot.xls",
            'nr_evalue': self.option('nr_blast_evalue'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        }
        self.annot_filter.set_options(options)
        self.annot_filter.on('start', self.set_step, {'start': self.step.annot_filter})
        self.annot_filter.on('end', self.set_step, {'end': self.step.annot_filter})
        self.annot_filter.on('end', self.set_output, "annot_filter")
        self.annot_filter.run()

    def run_annot_class(self):
        filter_dir = self.annot_filter.output_dir
        options = {
            # 'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml': filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': filter_dir + "/pfam/pfam_domain.filter.xls",
            "blast2go_annot": filter_dir + "/go/go_annot.xls.filter.xls",
            "des": self.rock_index.option('ptt').prop['path']
        }
        self.annot_class.set_options(options)
        self.annot_class.on('start', self.set_step, {'start': self.step.annot_class})
        self.annot_class.on('end', self.set_step, {'end': self.step.annot_class})
        self.annot_class.on('end', self.set_output, "annot_class")
        self.annot_class.run()

    def run_express(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option("fq_type") == "PE":
                self.libtype = "rf"
            else:
                self.libtype = "r"
            transcriptome = self.rockhopper.option("genome_fa").prop['path']
        else:
            self.libtype = None
            transcriptome = self.rock_index.option("query")
        opts = {
            "fastq" : self.qc.option("fq_list"),
            "method" : self.option("express_method").lower(),
            "libtype" : self.libtype,
            "transcriptome" : transcriptome,
            "id2name" : self.rock_index.option("ptt").prop["path"],
        }
        self.express.set_options(opts)
        self.express.on("end", self.set_output, "express")
        self.express.on('start', self.set_step, {'start': self.step.express})
        self.express.on('end', self.set_step, {'end': self.step.express})
        self.express.run()

    def deal_express_file(self):
        self.express_deal = self.express.output_dir + "/transcript.{}.matrix.dealed".format(self.option("exp_way").lower())
        exp_df = pd.read_table(self.express.output_dir + "/transcript.{}.matrix".format(self.option("exp_way").lower()), header =0)
        exp_df.drop(['gene_name', 'type'], axis=1, inplace=True)
        exp_df.to_csv(self.express_deal, sep="\t", index=False)

    def run_exp_pca(self):
        # time.sleep(60)
        self.logger.info("开始运行pca")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : self.express_deal
            }
        else:
            opts = {
                "exp" : self.express_deal
            }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.run()

    def run_exp_corr(self):
        # time.sleep(60)
        self.logger.info("开始运行聚类分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : self.express_deal
            }
        else:
            opts = {
                "exp" : self.express_deal
            }
        self.exp_corr.set_options(opts)
        self.exp_corr.on("end", self.set_output, "exp_corr")
        self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.run()

    def run_exp_venn(self):
        # time.sleep(60)
        self.logger.info("开始运行聚类分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "express_matrix" : self.express_deal,
                "group_table" : self.option('group_table')
            }
        else:
            opts = {
                "express_matrix" : self.express_deal,
                "group_table": self.option('group_table')
            }
        self.exp_venn.set_options(opts)
        self.exp_venn.on("end", self.set_output, "exp_venn")
        self.exp_venn.on('start', self.set_step, {'start': self.step.exp_venn})
        self.exp_venn.on('end', self.set_step, {'end': self.step.exp_venn})
        self.exp_venn.run()

    def run_diffexpress(self):
        self.deal_express_file()
        self.logger.info("开始运行基因差异表达分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            count_file = self.express.output_dir + "/transcript.count.matrix"
            fpkm_file = self.express_deal
            opts = {
                "count" : count_file,
                "exp" : fpkm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                "pvalue_padjust" : self.option("pvalue_padjust"),
                "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method")
            }
        else:
            count_file = self.express.output_dir + "/transcript.count.matrix"
            tpm_file = self.express_deal
            opts = {
                "count" : count_file,
                "exp" : tpm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                "pvalue_padjust" : self.option("pvalue_padjust"),
                "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method")
            }
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

    def run_snp(self):
        self.logger.info("开始运行snp步骤")
        opts = {
            "ref_genome_custom": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option('align_species') == 'is_ncbi' else self.option('genome_db').prop['path'],
            "ref_genome":  "customer_mode",
            "ref_gtf": self.rock_index.option('gtf'),
            "fq_list": self.qc.option("fq_list").prop["path"],
            "bamlist": self.mapping.option("bamlist"),
            "id2name": self.rock_index.option("ptt").prop['path'],
        }
        self.snp.set_options(opts)
        self.snp.on("start", self.set_step, {"start": self.step.snp_rna})
        self.snp.on("end", self.set_step, {"end": self.step.snp_rna})
        self.snp.on("end", self.set_output, "snp")
        self.snp.run()

    def run_map_assess(self):
        opts = {
            "bam": self.mapping.option("bam_output"),
            "bed": self.filecheck.option("bed"),
        }
        self.map_qc.set_options(opts)
        self.map_qc.on("start", self.set_step, {"start": self.step.map_qc})
        self.map_qc.on("end", self.set_step, {"end": self.step.map_qc})
        self.map_qc.on("end", self.set_output, "map_qc")
        self.map_qc.run()

    def run_srna(self):
        self.logger.info("开始运行srna")
        if self.option("align_species") == 'is_ncbi':
            opts = {
                "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                "input_file": os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff'),
                "type": "gff"
            }
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": self.option('gff_or_gtf_file').prop['path'],
                "type": self.option('gff_or_gtf')
            }
        opts.update({"group_list": self.option('group_table').prop['path'], "trimPairFq": self.qc.option('fq_list').prop['path']})
        self.srna.set_options(opts)
        self.srna.on("end", self.set_output, "srna")
        self.srna.on('start', self.set_step, {'start': self.step.srna})
        self.srna.on('end', self.set_step, {'end': self.step.srna})
        self.srna.run()

    def run_promote(self):
        self.logger.info("开始运行promote")
        opts = {
            "sequence": self.rock_index.option('query'),
            "assemble": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option('align_species') == 'is_ncbi' else self.option('genome_db').prop['path'],
            "rock_index": self.rock_index.option('ptt'),
            "sample": "all"
        }
        self.promote.set_options(opts)
        self.promote.on("end", self.set_output, "promote")
        self.promote.on('start', self.set_step, {'start': self.step.promote})
        self.promote.on('end', self.set_step, {'end': self.step.promote})
        self.promote.run()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error("需要移动到output目录的文件夹不存在", code = "15000120")
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

    def set_output(self, event):
        # pass
        obj = event["bind_object"]
        if event['data'] == 'qc':
            self.move2outputdir(obj.output_dir, 'QC_stat')
        if event['data'] == 'qc_stat_before':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
            self.logger.info("开始设置qc的输出目录")
        if event['data'] == 'qc_stat_after':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
        if event['data'] == 'mapping':
            self.move2outputdir(obj.output_dir, 'mapping')
            self.judge_pass()
        if event['data'] == 'map_qc':
            self.move2outputdir(obj.output_dir, 'map_qc')
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'snp')
            # self.end_number += 1
        if event['data'] == 'express':
            self.move2outputdir(obj.output_dir, 'express')
        if event['data'] == 'exp_pca':
            self.move2outputdir(obj.output_dir, 'exp_pca')
            # self.end_number += 1
        if event['data'] == 'exp_corr':
            self.move2outputdir(obj.output_dir, 'exp_corr')
            # self.end_number += 1
        if event['data'] == 'exp_venn':
            self.move2outputdir(obj.output_dir, 'exp_venn')
            # self.end_number += 1
        if event['data'] == 'diffexpress':
            self.move2outputdir(obj.output_dir, 'diffexpress')
            # self.end_number += 1
        if event['data'] == 'srna':
            self.move2outputdir(obj.output_dir, 'srna')
            # self.end_number += 1
        if event['data'] == 'annot_mapdb':
            self.move2outputdir(obj.output_dir, 'annot_mapdb')
        if event['data'] == 'annot_orfpfam':
            self.move2outputdir(obj.output_dir, 'annot_orfpfam')
        if event['data'] == 'promote':
            self.move2outputdir(obj.output_dir, 'promote')
            # self.end_number += 1
        if event['data'] == 'annot_class':
            self.move2outputdir(obj.output_dir, 'annot_class')
            # self.end_number += 1
        if event['data'] == 'rockhopper':
            self.move2outputdir(obj.output_dir, 'rockhopper')
        if event['data'] == 'rock_index':
            self.move2outputdir(obj.output_dir, 'rock_index')
        # if self.end_number == self.to_end-1:
        if not self.pass_judge:
            self.end()

    def end(self):
        if not self.pass_judge:
            self.run_api_qc_mapping()
            db = Config().get_mongo_client(mtype="project")[Config().get_mongo_dbname("project")]
            email = db['sg_task_email']
            email.update({"task_id": self.task_id}, {"$set": {"status": "5"}},
                       upsert=True)
            super(ProkrnaWorkflow, self).end()
        else:
            self.run_api()
            ## 更新一系列主表的字段，用于页面交互分析
            self.fa = self.workflow_output + "/Rockhopper/genome.feature.fa" if self.option('strand_specific') else self.workflow_output + "/Sequence_database/cds.fa"
            db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
            col = db["sg_task"]
            col.update({"task_id" : self.task_id}, {"$set": {"rock_index": self.workflow_output + "/Sequence_database/"}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"ref_gtf": self.workflow_output + "/Sequence_database/reshape.gtf"}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"ref_genome": self.ref_genome}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"strand_specific": self.option('strand_specific')}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"genome_id": self.option('genome_id')}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.fa}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"fastq": self.workflow_output + "/QC/cleandata/fq_list.txt"}}, upsert=True)
            col.update({"task_id": self.task_id}, {"$set": {"id2name": self.workflow_output + "/Sequence_database/ptt.bed"}},
                        upsert=True)
            col.update({"task_id": self.task_id}, {"$set": {"srna_files": self.workflow_output2 + "/srna/"}},
                       upsert=True)
            col1 = db["sg_annotation_stat"]
            col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)

            col2 = db["sg_exp"]
            col2.update({"task_id" : self.task_id}, {"$set": {"count_file": self.workflow_output + "/Express/ExpAnnalysis/transcript.count.matrix.xls", "exp_level" : self.exp_level}}, upsert=True)

            col3 = db["sg_specimen"]
            # col3.update({"task_id": self.task_id}, {"$set": {"fq_type": self.option('fq_type').lower()}}, upsert=True)
            col4 = db["sg_srna_fold"]
            col4.update({"task_id": self.task_id}, {"$set": {"result_dir": self.workflow_output2 + "/srna/srna_fold"}}, upsert=True)
            col5 = db["sg_snp"]
            col5.update({"task_id": self.task_id}, {"$set": {"result_dir": self.workflow_output + "/SNP/snp_annotation_detail.xls"}},
                        upsert=True)
            self.merge_annotation_exp_matrix() # 表达量表增加注释信息
            self.merge_annotation_diffexp_matrix() # 差异表达量表增加注释信息
            self.modify_output()

            super(ProkrnaWorkflow, self).end()

    def modify_output(self):
        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"
        #Background
        os.mkdir(target_dir + "/Background")
        genome_stat = self.genome_stat
        os.link(genome_stat, target_dir + "/Background/" + os.path.basename(genome_stat))
        #seq_db
        os.mkdir(target_dir + "/Sequence_database")
        for file in os.listdir(os.path.join(origin_dir, 'rock_index')):
            os.link(os.path.join(origin_dir, 'rock_index', file), os.path.join(target_dir, 'Sequence_database', file))
        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        os.mkdir(target_dir + "/QC")
        os.mkdir(target_dir + "/QC/cleandata")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        for file in os.listdir(origin_dir + "/QC_stat/sickle_dir"):
            file_path = os.path.join(origin_dir + "/QC_stat/sickle_dir", file)
            os.link(file_path, target_dir + "/QC/cleandata/" + file)
        # Align
        # AlignBam
        os.makedirs(target_dir + "/Align/AlignBam")
        for file in os.listdir(origin_dir + "/mapping/bam"):
            file_path = os.path.join(origin_dir + "/mapping/bam", file)
            os.link(file_path, target_dir + "/Align/AlignBam/" + file)
        # AlignStat
        os.makedirs(target_dir + "/Align/AlignStat")
        for file in os.listdir(origin_dir + "/mapping/stat"):
            file_path = os.path.join(origin_dir + "/mapping/stat", file)
            os.link(file_path, target_dir + "/Align/AlignStat/" + file.split(".stat")[0] + "_align_stat.txt")
        # QualityAssessment
        os.makedirs(target_dir + "/Align/QualityAssessment")
        if os.path.exists(origin_dir + "/map_qc/chr_stat"):
            for file in os.listdir(origin_dir + "/map_qc/chr_stat"):
                file_path = os.path.join(origin_dir + "/map_qc/chr_stat", file)
                os.link(file_path, target_dir + "/Align/QualityAssessment/" + file.split("bam_chr_stat.xls")[0] + "chr_distribution.xls")
        if os.path.exists(origin_dir + "/map_qc/coverage"):
            for file in os.listdir(origin_dir + "/map_qc/coverage"):
                file_path = os.path.join(origin_dir + "/map_qc/coverage", file)
                os.link(file_path, target_dir + "/Align/QualityAssessment/" + file.split("coverage_")[1])
        if os.path.exists(origin_dir + "/map_qc/distribution"):
            for file in os.listdir(origin_dir + "/map_qc/distribution"):
                file_path = os.path.join(origin_dir + "/map_qc/distribution", file)
                os.link(file_path, target_dir + "/Align/QualityAssessment/" + file.split(".reads_distribution.txt")[0] + ".region_distribution.xls")
        if os.path.exists(origin_dir + "/map_qc/saturation"):
            files = glob.glob(origin_dir + "/map_qc/saturation/" + "*.eRPKM.xls.cluster_percent.xls")
            for file in files:
                file_name = os.path.basename(file)
                os.link(file, target_dir + "/Align/QualityAssessment/" + file_name.split(".eRPKM.xls.cluster_percent.xls")[0].split("satur_")[1] + ".saturation.xls")
        #Annotation
        os.makedirs(target_dir + "/Annotation")
        CopyFile().linkdir(origin_dir + "/annot_class", target_dir + "/Annotation")
        CopyFile().linkdir(origin_dir + "/annot_mapdb", target_dir + "/Annotation/annot_mapdb")
        CopyFile().linkdir(origin_dir + "/annot_orfpfam", target_dir + "/Annotation/annot_mapdb/pfam")
        if os.path.exists(target_dir + "/Annotation/" + os.path.basename(genome_stat)):
            os.remove(target_dir + "/Annotation/" + os.path.basename(genome_stat))
            os.link(genome_stat, target_dir + "/Annotation/" + os.path.basename(genome_stat))
        else:
            os.link(genome_stat, target_dir + "/Annotation/" + os.path.basename(genome_stat))
        #Expression
        os.mkdir(target_dir + "/Express")
        os.mkdir(target_dir + "/Express/ExpAnnalysis")
        if os.path.exists(self.express.output_dir + "/transcript.count.matrix"):
            trans_count = os.path.join(self.express.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
        if os.path.exists(self.express.output_dir + "/transcript.tpm.matrix"):
            trans_tpm = os.path.join(self.express.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
        if os.path.exists(self.express.output_dir + "/transcript.fpkm.matrix"):
            trans_fpkm = os.path.join(self.express.output_dir + "/transcript.fpkm.matrix")
            os.link(trans_fpkm, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.xls")
        if os.path.exists(self.express.output_dir + "/transcript.tpm.matrix.annot.xls"):
            trans_tpm_anno = os.path.join(self.express.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
        if os.path.exists(self.express.output_dir + "/transcript.fpkm.matrix.annot.xls"):
            trans_fpkm_anno = os.path.join(self.express.output_dir + "/transcript.fpkm.matrix.annot.xls")
            os.link(trans_fpkm_anno, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls")
        os.mkdir(target_dir + "/Express/ExpCorr")
        gene_corr = os.path.join(self.exp_corr.output_dir + "/sample_correlation.xls")
        os.link(gene_corr, target_dir + "/Express/ExpCorr/sample_correlation.xls ")
        if self.option("group_table").prop["sample_number"] > 2:
            os.mkdir(target_dir + "/Express/ExpPCA")
            gene_pca = os.path.join(self.exp_pca.output_dir + "/PCA.xls")
            os.link(gene_pca, target_dir + "/Express/ExpPCA/PCA.xls")
            gene_ratio = os.path.join(self.exp_pca.output_dir + "/Explained_variance_ratio.xls")
            os.link(gene_ratio, target_dir + "/Express/ExpPCA/Explained_variance_ratio.xls")
        # DiffExpress
        os.mkdir(target_dir + "/DiffExpress")
        for file in os.listdir(self.diffexpress.output_dir):
            diffexpress = os.path.join(self.diffexpress.output_dir, file)
            os.link(diffexpress, target_dir + "/DiffExpress/" + file)
        # SNP
        os.makedirs(target_dir + "/SNP")
        CopyFile().linkdir(self.work_dir + "/snp_tmp", target_dir + "/SNP")
        self.merge1(target_dir + "/SNP/data_anno_pre.xls", target_dir + "/SNP/snp_annotation.xls")
        if os.path.exists(target_dir + "/SNP/snp_annotation.xls"):
            os.remove(target_dir + "/SNP/snp_annotation.xls")
        if os.path.exists(target_dir + "/SNP/data_anno_pre.xls"):
            os.remove(target_dir + "/SNP/data_anno_pre.xls")
        if self.option("strand_specific"):
            os.makedirs(target_dir + "/Rockhopper")
            for file in os.listdir(self.rockhopper.output_dir):
                rockhopper = os.path.join(self.rockhopper.output_dir, file)
                os.link(rockhopper, target_dir + "/Rockhopper/" + file)
            os.makedirs(target_dir + "/srna")
            for file in os.listdir(self.srna.output_dir):
                if u'ockhopper' not in file:
                    srna = os.path.join(self.srna.output_dir, file)
                    CopyFile().linkdir(srna, target_dir + "/srna/" + file)
                else:
                    os.makedirs(target_dir + "/srna/srna_predict")
                    srna = os.path.join(self.srna.output_dir, file, 'genome.predicted_RNA.fa')
                    os.link(srna, target_dir + "/srna/srna_predict" + '/genome.predicted_RNA.fa')

        sdir = self.add_upload_dir(target_dir)
        sdir.add_regexp_rules([

            [r"Align/AlignBam/.*\.bam", "", "样本比对bam文件", 1],
            [r"Align/AlignStat/.*_align_stat\.txt", "", "样本比对统计结果表"],
            [r"Align/QualityAssessment/.*\.chr_distribution\.xls", "", "Reads在不同染色体的分布统计表"],
            [r"Align/QualityAssessment/.*\.geneBodyCoverage\.txt", "", "基因覆盖度分布结果文件"],
            [r"Align/QualityAssessment/.*\.region_distribution\.xls", "", "Reads在不同区域的分布统计表"],
            [r"Align/QualityAssessment/.*\.saturation\.xls", "", "测序饱和度分析结果文件"],
            [r"Express/ExpAnnalysis/transcript.count.matrix.xls", "", "transcript count表达定量结果表"],
            [r"Express/ExpAnnalysis/transcript\..*m\.matrix.xls", "", "transcript 表达定量结果表"],
            [r"Express/ExpAnnalysis/transcript\..*m\.matrix.annot.xls", "", "transcript tpm表达定量注释结果表"],
            [r"DiffExpress/.*_vs_.*\.xls", "xls", "差异分析结果表"],
            [r"DiffExpress/.*\.DE\.list", "xls", "差异基因列表"],
            [r"DiffExpress/.*summary\.xls", "xls", "差异统计结果表"],
            [r"DiffExpress/.*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表"],
            [r"DiffExpress/.*_vs_.*\..*annot\.xls", "xls", "差异统计注释结果表"],
            ])
        rel_path = [
            [".", "", "流程分析结果目录"],
            ["Background", "", "项目背景目录"],
            ["Background/genome_stat.xls", "", " 参考基因组注释信息表"],
            ["Sequence_database", "", "序列文件数据库", 1],
            ["QC", "", "测序数据统计与质控结果"],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表"],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表"],
            ["QC/cleandata", "", "质控结果目录", 1],
            ["Align", "", "比对结果目录"],
            ["Align/AlignBam", "", "比对结果bam目录", 1],
            ["Align/AlignStat", "", "比对结果统计目录"],
            ["Align/QualityAssessment", "", "比对结果整体评估目录"],
            ["Annotation/annot_class/anno_stat", "", "CDS注释统计结果目录"],
            ["Annotation/annot_class/anno_stat/all_anno_detail.xls", "", "CDS注释详细信息的表"],
            ["Annotation/annot_class/anno_stat/all_annotation_statistics.xls", "", "CDS注释信息占比统计"],
            ["Annotation/annot_class/anno_stat/venn", "", "CDS在各个数据库注释结果统计VENN图目录"],
            ["Annotation/annot_class/anno_stat/venn/pfam_venn.txt", "", "有与pfam数据库相关信息的CDS"],
            ["Annotation/annot_class/anno_stat/venn/kegg_venn.txt", "", "有与kegg数据库相关信息的CDS"],
            ["Annotation/annot_class/anno_stat/venn/go_venn.txt", "", "有与go数据库相关信息的CDS"],
            ["Annotation/annot_class/anno_stat/venn/cog_venn.txt", "", "有与cog数据库相关信息的CDS"],
            ["Annotation/annot_class/anno_stat/venn/nr_venn.txt", "", "有与nr数据库相关信息的CDS"],
            ["Annotation/annot_class/anno_stat/venn/swissprot_venn.txt", "", "有与swissprot数据库相关信息的CDS"],
            ["Annotation/annot_class/go", "", "鉴定CDS与GO数据库比对结果详情目录"],
            ["Annotation/annot_class/go/query_gos.list", "", "鉴定CDS与GO的对应关系"],
            ["Annotation/annot_class/go/go_detail.xls", "", "鉴定CDS GO注释详情表，包含GO等级信息"],
            ["Annotation/annot_class/go/go12level_statistics.xls", "", "鉴定CDS的GO1、2层级的详细信息"],
            ["Annotation/annot_class/go/go123level_statistics.xls", "", "鉴定CDS的GO1、2、3层级的详细信息"],
            ["Annotation/annot_class/go/go1234level_statistics.xls", "", "鉴定CDS的GO1、2、3、4层级的详细信息"],
            ["Annotation/annot_class/cog", "", "鉴定CDS与COG数据库比对结果详情目录"],
            ["Annotation/annot_class/cog/cog_table.xls", "", "鉴定CDS序列文件与COG比对结果表"],
            ["Annotation/annot_class/cog/cog_summary.xls", "", "COG数据库注释结果汇总"],
            ["Annotation/annot_class/cog/cog_list.xls", "", "鉴定CDS与COG的对应关系"],
            ["Annotation/annot_class/kegg", "", "鉴定CDS与KEGG数据库比对结果详情目录"],
            ["Annotation/annot_class/kegg/pathways.tar.gz", "", "KEGG通路注释结果图片目录"],
            ["Annotation/annot_class/kegg/pathway_table.xls", "", "每行以通路为单位显示通路注释详情表"],
            ["Annotation/annot_class/kegg/kegg_table.xls", "", "每行以CDS为单位显示通路注释详情表"],
            ["Annotation/annot_class/kegg/pid.txt", "", "每个注释通路和映射上去的CDS联系表"],
            ["Annotation/annot_class/blast_xml", "", "CDS比对注释结果目录", 1],
            ["Annotation/annot_class/blast_xml/vs_kegg.xml", "", "鉴定CDS与数据库比对结果文件目录"],
            ["Annotation/annot_class/blast_xml/vs_eggnog.xml", "", "鉴定CDS与EGGNOG数据库比对结果xml形式"],
            ["Annotation/annot_class/blast_xml/pfam_domain", "", "鉴定CDS与pfam数据库比对结果xml形式"],
            ["Annotation/annot_class/blast_xml/vs_nr.xml", "", "鉴定CDS与NR数据库比对结果xml形式"],
            ["Annotation/annot_class/blast_xml/vs_swissprot.xml", "", "鉴定CDS与NR数据库比对结果xml形式"],
            ["Annotation/annot_class/annot_mapdb", "", "数据库原始比对结果(evalue 1e-3)"],
            ["Annotation/annot_class", "", "注释结果目录"],
            ["Express", "", "表达量分析结果目录"],
            ["Express/ExpAnnalysis", "", "表达定量结果目录"],
            ["Express/ExpCorr", "", "样本间相关性分析"],
            ["Express/ExpCorr/sample_correlation.xls ", "", "样本间相关性分析矩阵表"],
            ["Express/ExpPCA", "", "样本间PCA分析"],
            ["Express/ExpPCA/PCA.xls", "", "样本间PCA分析结果表"],
            ["Express/ExpPCA/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表"],
            ["DiffExpress", "", "表达量差异分析"],
            ["SNP", "", "SNP/InDel分析结果文件"],
            ["SNP/snp_annotation_statistics.xls", "", "snp分析结果注释详情表格"],
            ["SNP/snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格"],
            ["SNP/snp_freq_statistics.xls", "", "SNP频率统计结果表格"],
            ["SNP/snp_depth_statistics.xls", "", "SNP深度统计结果表格"],
            ["SNP/snp_position_distribution.xls", "", "SNP不同区域布析结果表格"],
            ["SNP/indel_position_distribution.xls", "", "InDel不同区域布析结果表格"],
            ]
        if self.option("strand_specific"):
            extra_rel_path = [
                ["Rockhopper", "", "Rockhopper分析结果文件"],
                ["Rockhopper/cds.bed", "", "cds的bed文件"],
                ["Rockhopper/cds.fa", "", "cds的fasta文件"],
                ["Rockhopper/cds.faa", "", "cds翻译的氨基酸序列文件"],
                ["Rockhopper/cds.faa", "", "cds翻译的氨基酸序列文件"],
                ["srna", "", "sRNA分析结果文件"],
                ["srna/srna_predict", "", "sRNA预测结果文件"],
                ["srna/srna_annot", "", "sRNA注释分析结果文件"],
                ["srna/srna_fold", "", "sRNA二级结构预测结果文件"],
                ["srna/srna_target", "", "sRNA靶基因预测结果文件"],
            ]
            rel_path += extra_rel_path
        sdir.add_relpath_rules(rel_path)

    def run_api(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.prok_rna')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.stop_timeout_check()
        self.export_genome_info()
        self.export_rock_index()
        self.export_qc()
        self.export_map_assess()
        if self.option("strand_specific") == True:
            self.export_srna()
            self.export_structure(srna = True)
        self.export_structure(promote = True)
        self.export_snp()
        # self.export_promote()
        self.export_annotation()
        self.export_expression()
        self.logger.info("导表完成")

    def run_api_qc_mapping(self):
        self.logger.info('经判断，序列比对或rRNA比例不达标，流程提前终止')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.prok_rna')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.export_qc()

    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("prok_rna.genome_info")
        genome_path = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option('align_species') == 'is_ncbi' else self.option('genome_db').prop['path']
        all_seq = ''
        with open(genome_path, 'r') as fna_r:
            for seq in fna_r.read().split('\n>'):
                all_seq += ''.join(seq.strip().split('\n')[1:])
        self.genome_size = len(all_seq)
        gc = 0
        for i in all_seq:
            if i.lower() == 'g' or i.lower() == 'c':
                gc += 1
        self.gc_content = float(gc)/self.genome_size*100
        self.gc_content = round(self.gc_content, 2)
        # self.gc_content = (all_seq.lower().count('g') + all_seq.lower().count('c'))/self.genome_size*100
        # gc_content = int(gc_content)
        file_path = genome_path
        self.species_name = self.option('species_name')
        self.accession = self.option('genome_id')
        self.api_geno.add_genome_info(file_path=file_path,species_name=self.species_name, accession=self.accession, genome_size=self.genome_size,gc_content=self.gc_content,major=False)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("prok_rna.prok_rna_qc")
        qc_stat_before = self.qc_stat_before.output_dir
        qc_stat_after = self.qc_stat_after.output_dir
        list_txt = os.path.join(self.option('fastq_dir').prop['path'], 'list.txt')
        fq_type = self.option("fq_type").upper()
        self.api_qc.add_samples_info(qc_stat_before, qc_stat_after, list_txt, fq_type=fq_type)
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat"
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat"  # 将qc前导表加于该处
        self.api_qc.add_gragh_info(quality_stat_before, "before")
        # self.api_qc.add_bam_path(self.mapping.output_dir + "/bam")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(self.option("group_table").prop["path"])
        self.logger.info("group_detail为：" + str(self.group_detail))
        self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control_file").prop["path"], self.group_id)
        self.compare_detail = compare_detail
        self.api_qc.add_bam_path(self.workflow_output)
        # map_stat = os.path.join(self.mapping.output_dir, 'stat')
        # self.api_qc.add_bowtie2_mapping_stat(map_stat)
        rrna_stat = os.path.join(self.qc_stat_after.output_dir, 'RfamStat')
        self.api_qc.add_rfam_stat(rrna_stat)

    @time_count
    def export_srna(self):
        self.api_srna = self.api.api("prok_rna.srna")
        predict_dir = os.path.join(self.srna.output_dir, "rockhopper/genome.predicted_RNA.bed.xls")
        fold_dir = os.path.join(self.srna.output_dir, "srna_fold")
        annot_dir = os.path.join(self.srna.output_dir, "srna_annot")
        target_dir = os.path.join(self.srna.output_dir, "srna_target")
        params = self.task_id
        predict_id = self.api_srna.add_main_table('sg_srna_predict', params=params)
        self.api_srna.add_predict_rna_bed(predict_id, predict_dir)
        annot_id = self.api_srna.add_main_table('sg_srna_anno', params=params)
        self.api_srna.add_srna_annot(annot_id, annot_dir)
        fold_id = self.api_srna.add_main_table('sg_srna_fold', params=params)
        self.api_srna.add_srna_fold(fold_id, fold_dir)
        target_id = self.api_srna.add_main_table('sg_srna_target', params=params)
        # if os.path.getsize(os.path.join(target_dir, 'combine_RNAplex_RNAhybrid')) < 2000000:
        #     self.api_srna.add_srna_target(target_id, target_dir)
        # else:
        #     self.logger.info('srna靶基因预测结果太多，为避免导表出错，只导入差异srna的靶基因')
        #     self.api_srna.add_srna_target_filter_bydiff(target_id, target_dir, self.diffexpress.output_dir)
        self.api_srna.add_srna_target(target_id, target_dir)
        # self.api_srna.run(predict_dir, annot_dir, fold_dir, target_dir, params)

    def export_rock_index(self):
        self.api_srna = self.api.api("prok_rna.srna")
        db_path = self.rock_index.output_dir
        genome_path = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option('align_species') == 'is_ncbi' else self.option('genome_db').prop['path']
        self.api_srna.build_seq_database(db_path + '/cds.fa.db.sqlite3', db_path + '/cds.fa', table_name = 'seq_gene', type_ = 'gene')
        self.api_srna.build_seq_database(db_path + '/cds.faa.db.sqlite3', db_path + '/cds.faa', table_name = 'seq_protein', type_ = 'protein')
        self.api_srna.build_seq_database(db_path + '/genome_fna.db.sqlite3', genome_path, table_name='seq_genome',
                                       type_='genome')

    @time_count
    def export_map_assess(self):
        gevent.sleep()
        self.api_map = self.api.api("prok_rna.prok_rna_qc")
        stat_dir = self.mapping.output_dir + "/stat"
        self.api_map.add_bowtie2_mapping_stat(stat_dir)
        analysis = self.option("map_assess_method").split(",")
        if "saturation" in analysis:
            file_path = self.map_qc.output_dir + "/saturation"
            self.api_map.add_rpkm_table(file_path)
        if "coverage" in analysis:
            coverage = self.map_qc.output_dir + "/coverage"
            self.api_map.add_coverage_table(coverage)
        if "distribution" in analysis:
            distribution = self.map_qc.output_dir + "/distribution"
            self.api_map.add_distribution_table(distribution)
        if "chr_stat" in analysis:
            chrom_distribution = self.map_qc.output_dir + "/chr_stat"
            self.api_map.add_chrom_distribution_table(chrom_distribution)

    @time_count
    def export_annotation(self):
        self.api_annotation = self.api.api("prok_rna.prokrna_annotation")
        annot_dir = self.annot_class.output_dir
        self.genome_stat = annot_dir + "/genome_stat.xls"
        with open(self.rock_index.option('query').prop['path'], 'r') as cds_r:
            cds_num = len(cds_r.read().strip().split('\n>'))
        with open(self.genome_stat, 'w') as genome_s:
            genome_s.write('Organism' + '\t' + self.species_name + '\n')
            genome_s.write('Accession' + '\t' + self.accession + '\n')
            genome_s.write('Size' + '\t' + str(self.genome_size) + '\n')
            genome_s.write('GC' + '\t' + str(self.gc_content) + '\n')
            genome_s.write('CDS' + '\t' + str(cds_num) + '\n')
        params = {
            "nr_evalue": str(self.option("nr_blast_evalue")),
            "swissprot_evalue":str(self.option("swissprot_blast_evalue")),
            "cog_evalue": str(self.option("cog_blast_evalue")),
            "kegg_evalue": str(self.option("kegg_blast_evalue")),
            "pfam_evalue": str(self.option("pfam_blast_evalue")),
        }
        self.api_annotation.run_prok_rna(self.genome_stat, annot_dir, params)

    @time_count
    def export_snp(self):
        gevent.sleep()
        self.logger.info("开始进行Snpfinal的导表")
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            submit_location="snp",
            task_type=2,
            method_type="samtools"
        )
        self.api_snp = self.api.api("prok_rna.ref_snp")
        snp_anno = self.snp.output_dir
        df = pd.read_table(snp_anno + '/snp_anno.xls', header=0, sep='\t', low_memory=False)
        df.columns = df.columns.map(lambda x: x.split("bam/")[1] if "bam" in x else x)
        df.to_csv(snp_anno + '/snp_anno', header=True, sep="\t", index=False)
        if os.path.exists(snp_anno + '/snp_anno.xls'):
            os.remove(snp_anno + '/snp_anno.xls')
        os.link(snp_anno + '/snp_anno', snp_anno + '/snp_anno.xls')
        os.remove(snp_anno + '/snp_anno')
        if os.path.exists(self.work_dir + "/snp_tmp"):
            shutil.rmtree(self.work_dir + "/snp_tmp")
        os.mkdir(self.work_dir + "/snp_tmp")
        if os.path.exists(self.work_dir + "/snp_tmp/snp_annotation.xls"):
            os.remove(self.work_dir + "/snp_tmp/snp_annotation.xls")
        os.link(self.snp.work_dir + "/snp_annotation.xls", self.work_dir + "/snp_tmp/snp_annotation.xls")
        self.api_snp.add_snp_main(snp_anno=snp_anno, params=params,
                                  task_id=task_id, method_type="samtools", project_sn=project_sn, new_output=self.work_dir + "/snp_tmp")

    @time_count
    def export_expression(self):
        gevent.sleep()
        all_exp = self.api.api("prok_rna.all_exp")
        group_dict = self.option('group_table').prop['group_dict']
        quant_method = self.option('express_method')
        group_id = self.group_id
        control_id = self.control_id
        task_id = self.task_id
        project_sn = self.project_sn

        # add exp matrix
        ## 还需确认路径信息，因为后续交互需要用到express.output_dir
        exp_output = self.express.output_dir
        self.exp_level = 'mRNA+sRNA' if self.option('strand_specific') else 'mRNA'
        if self.option("exp_way").lower() == "tpm":
            params = dict(
                task_id=task_id,
                submit_location="exp_detail",
                task_type=2,
                method=quant_method,
                exp_type='TPM',
            )
        else:
            params = dict(
                task_id=task_id,
                submit_location="exp_detail",
                task_type=2,
                method=quant_method,
                exp_type='FPKM',
            )
        # print('litangjian')
        # print(params)
        if self.option("exp_way").lower() == "tpm":
            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            trans_exp_id = all_exp.add_exp(
                exp_matrix, quant_method='rsem', group_dict=group_dict, main_id_m=None, main_id_s=None,
                main_id_ms=None, group_id=group_id,
                add_distribution=False, exp_type='TPM', project_sn=project_sn, task_id=task_id,
                lib_type=self.libtype, params=params)
            params_distribution = dict(
                task_id=task_id,
                exp_id=str(trans_exp_id),
                group_dict=group_dict,
                group_id=str(group_id),
                submit_location="expgraph",
                task_type=2,
                exp_level=self.exp_level,
            )
            all_exp.add_distribution(exp_matrix, group_dict, params=params_distribution,
                              quant_method=quant_method, project_sn=project_sn, task_id=task_id, exp_level=self.exp_level)
        else:
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            trans_exp_id = all_exp.add_exp(exp_matrix, quant_method='rsem', group_dict=group_dict, main_id_m=None, main_id_s=None,
                main_id_ms=None, group_id=group_id,
                add_distribution=True, exp_type='FPKM', project_sn=project_sn, task_id=task_id,
                lib_type=self.libtype, params=params)
            params_distribution = dict(
                task_id=task_id,
                exp_id=str(trans_exp_id),
                group_dict=group_dict,
                group_id=str(group_id),
                submit_location="expgraph",
                task_type=2,
                exp_level=self.exp_level,
            )
            all_exp.add_distribution(exp_matrix, group_dict, params=params_distribution,
                              quant_method=quant_method, project_sn=project_sn, task_id=task_id, exp_level=self.exp_level)

        # add gene corr
        corr_output = self.exp_corr.work_dir
        params = dict(
            task_id=task_id,
            submit_location='expcorr',
            task_type=2,
            exp_id=str(trans_exp_id),
            group_id=str(group_id),
            exp_level=self.exp_level,
            group_dict=group_dict,
            scm="complete",
            scd="euclidean",
            # quant_method=quant_method,
            corr_method="pearson",
        )
        all_exp.add_exp_corr2(corr_output, exp_level=self.exp_level, params=params,
                              project_sn=project_sn, task_id=task_id)

        # add gene venn
        venn_output = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
        params = dict(
            task_id=task_id,
            submit_location='expvenn',
            task_type=2,
            exp_id=str(trans_exp_id),
            group_id=str(group_id),
            exp_level=self.exp_level,
            group_dict=group_dict,
            # quant_method=quant_method,
        )
        try:
            all_exp.add_exp_venn(venn_output, quant_method=quant_method, exp_level=self.exp_level, params=params,
                              project_sn=project_sn, task_id=task_id)
        except:
            self.logger.info('venn没有导表')

        # add gene pca
        if self.option("group_table").prop["sample_number"] > 2:
            pca_output = self.exp_pca.output_dir
            params = dict(
                task_id=task_id,
                submit_location="exppca",
                task_type=2,
                exp_id=str(trans_exp_id),
                group_id=str(group_id),
                exp_level=self.exp_level,
                group_dict=group_dict,
                # quant_method=quant_method,
            )
            all_exp.add_exp_pca2(pca_output, exp_level=self.exp_level,
                                 params=params, project_sn=project_sn, task_id=task_id)

        # add gene diff
        diff_output = self.diffexpress.output_dir
        exp_id, exp_level = trans_exp_id, self.exp_level
        diff_method = self.option('diff_method')
        stat_type = self.option('pvalue_padjust')
        params = dict(
            task_id=task_id,
            submit_location="diff_detail",
            task_type=2,
            exp_id=str(exp_id),
            group_id=str(group_id),
            control_id=str(control_id),
            exp_level=exp_level,
            group_dict=group_dict,
            fc=str(self.option('fc')),
            correct_method=self.option('padjust_way'),
            stat_type=stat_type,
            stat_cutoff=self.option('diff_fdr_ci'),
            # quant_method=quant_method,
            diff_method=diff_method,
        )
        all_exp.add_diffexp(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
                            exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
                            project_sn=project_sn, task_id=task_id, params=params,
                            pvalue_padjust=stat_type)

    def export_structure(self, srna=False, promote=False):
        structure = self.api.api("prok_rna.gene_structure")
        if srna:
            structure.add_rock_structure(self.rockhopper.output_dir)
        if promote:
            structure.add_promote_structure(self.promote.output_dir + '/all_promoter_result.xls')

    # 添加注释信息
    def merge_annotation_exp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        exp_output = self.express.output_dir
        annot = os.path.join(self.output_dir, 'annot_class/anno_stat/all_anno_detail.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        trans_annot_pd = all_annot.reset_index().set_index('gene_id')

        if self.option("exp_way").lower() == "tpm":
            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
        else:
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
        trans_pd = pd.read_table(exp_matrix, header=0, index_col=0)
        trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
        trans_out = os.path.join(exp_output, 'transcript.{}.matrix.annot.xls'.format(self.option('exp_way').lower()))
        trans_result.to_csv(trans_out, header=True, index=True, sep='\t')

    def merge_annotation_diffexp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        annot = os.path.join(self.output_dir, 'annot_class/anno_stat/all_anno_detail.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        diff_output = self.diffexpress.output_dir

        duplicate_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls') ## 为了防止流程重跑的时候反复增加注释结果表
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + '/' + '*_vs_*.*.xls')
        for each in target_files:
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_pd, all_annot], axis=1)
            gene_out = each.split('.xls')[0] + '.annot.xls'
            gene_result.to_csv(gene_out, header=True, index=True, sep='\t')

    def merge1(self, x, y):
        df_anno_pre1 = pd.read_table(x, header=0, sep="\t", low_memory=False)
        tmp_list = df_anno_pre1.columns[:-14].tolist()
        tmp_list.append(df_anno_pre1.columns[-1])
        df_anno_pre1_select = df_anno_pre1.loc[:, tmp_list]
        df_anno_pre1_select["index1"] = df_anno_pre1['alt'].apply(str) + df_anno_pre1['anno'].apply(str) + \
        df_anno_pre1['chrom'].apply(str) + df_anno_pre1["end"].apply(str) + df_anno_pre1["start"].apply(str) + df_anno_pre1["ref"].apply(str)
        df_anno_pre2 = pd.read_table(y, header=0, sep="\t", low_memory=False)
        df_anno_pre2_select_pre = df_anno_pre2.loc[:, df_anno_pre2.columns[:13]]
        df_anno_pre2_select_pre.rename(columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno", "END": "End", "START": "Start", "REF": "Ref",
                                                "MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
        df_anno_pre2_select = df_anno_pre2_select_pre[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt", "Total depth", "QUAL",
                                                      "Anno", "MUT type", "MUT info"]]
        df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2['CHROM'].apply(str) + \
                                        df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + df_anno_pre2["REF"].apply(str)
        df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
        df_join_pre.drop(columns=['index1'], inplace=True)
        df_join_pre.to_csv(self.work_dir + "/upload_results/SNP/snp_annotation_detail.xls", sep="\t", index=False)

    def judge_pass(self):
        if self.option('mapping_sample_percent') or self.option('rrna_sample_percent'):
            sample_num = self.option("group_table").prop["sample_number"]
            mapping_files = os.listdir(os.path.join(self.mapping.output_dir, 'stat'))
            if len(mapping_files) != sample_num:
                self.set_error("不是所有的样本都有mapping统计结果", code = "15000121")
            rrna_files = os.listdir(os.path.join(self.qc_stat_after.output_dir, 'RfamStat'))
            if len(rrna_files) != sample_num:
                self.set_error("不是所有的样本都有rRNA比例统计结果", code = "15000122")
            map_failed = 0
            for map in mapping_files:
                with open(os.path.join(self.mapping.output_dir, 'stat',map), 'r') as map_r:
                    for line in map_r.readlines():
                        if u'overall alignment rate' in line:
                            if float(line.strip().split('%')[0]) < float(self.option('mapping_ratio')):
                                map_failed += 1
            if float(map_failed)/sample_num*100 > float(self.option('mapping_sample_percent')):
                self.pass_judge = False
                return

            rrna_failed = 0
            for rrna in rrna_files:
                with open(os.path.join(self.qc_stat_after.output_dir, 'RfamStat', rrna), 'r') as rrna_r:
                    _ = rrna_r.readline()
                    rrna_ratio = rrna_r.readline().strip().split('\t')[1]
                    if float(rrna_ratio) > float(self.option('rrna_ratio')):
                        rrna_failed += 1
            if float(rrna_failed)/sample_num*100 > float(self.option('rrna_sample_percent')):
                self.pass_judge = False
                return



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "prok_rna_workflow" + str(random.randint(1, 10000)),
    #         "type": "workflow",
    #         "name": "prok_rna.prokrna",
    #         "instant": False,
    #         "options": dict(
    #             fq_type="PE",
    #             quality_score_system="phred+33",
    #             species_name="Yersinia enterocolitica subsp. enterocolitica 8081",
    #             genome_id="GCF_000009345.1",
    #             fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/prok_rna_fastq_dir",
    #             group_table='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/group.txt',
    #             control_file='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/control.txt',
    #         )
    #     }
    #     # data['options']['method'] = 'rsem'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     # #
    #     # data['id'] += '1'
    #     # data['options']['method'] = 'salmon'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     #
    #     data['id'] += '_fyt'
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()

    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "prok_rna_workflow_" + id,
            "type": "workflow",
            "name": "prok_rna.prokrna",
            "options": dict(
                    fq_type="PE",
                    quality_score_system="phred+33",
                    species_name="Yersinia enterocolitica subsp. enterocolitica 8081",
                    genome_id="GCF_000009345.1",
                    align_species="is_upload",
                    genome_db="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/GCF_000009345.1.fna",
                    gff_or_gtf="gtf",
                    gff_or_gtf_file="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/reshape.gtf",
                    fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/prok_rna_fastq_dir",
                    group_table='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/group.txt',
                    control_file='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/control.txt',
                    mapping_ratio='40',
                    mapping_sample_percent='90',
                    rrna_ratio='90',
                    rrna_sample_percent='90',
                )
        }

        info = worker.add_task(data)
        print(info)


if __name__ == '__main__':
    unittest.main()
