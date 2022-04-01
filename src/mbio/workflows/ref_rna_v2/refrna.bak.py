# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""有参转录组一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import subprocess
import datetime
import json
import shutil
import re
import time
import gevent
import functools
from biocluster.config import Config
import pandas as pd
from mbio.packages.ref_rna_v2.copy_file import CopyFile


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


class RefrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(RefrnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 单样本 OR 多样本
            {"name": "sample_num", "type": "string", "default": "multiple"},  # 测序类型，single OR multiple
            ## 基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "strand_specific", "type": "bool", "default": False}, # 当为PE测序时，是否为链特异性
            {"name": "strand_dir", "type": "string", "default": "forward"}, # 当链特异性时为True时，forward or reverse
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "taxonmy", "type":"string", "default": "Animal"}, # 物种类别
            {"name": "ref_genome", "type": "string", "default": "Custom"},  # 参考基因组，具体物种名称
            {"name": "genome_version", "type": "string", "default": "Custom"},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": "Custom"},  # 参考基因组注释版本
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"}, # 对照表
            {"name": "fastp", "type": "int", "default": 0},  # 是否选用fastp进行质控
            ## 高级参数设置
            # 分析水平
            {"name": "level", "type": "string", "default": "gene"}, # 分析水平，gene or genetrans
            # 序列比对分析
            {"name": "align_method", "type": "string", "default": "Hisat"},  # 比对方法，Tophat or Hisat
            {"name": "map_assess_method", "type": "string", "default": "saturation,distribution,coverage,chr_stat"},#质量评估分析项
            # 转录本组装
            {"name": "is_assemble", "type": "bool", "default": True}, # 是否进行有参拼接
            {"name": "assemble_method", "type": "string", "default": "stringtie"},# 拼接方法，Cufflinks or Stringtie
            # 转录组功能注释
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            {"name": "kegg_database", "type": "string", "default": "All"},  # kegg注释库类型
            {"name": "nr_evalue", "type": "string", "default": "1e-3"},
            {"name": "kegg_evalue", "type": "string", "default": "1e-3"},
            {"name": "cog_evalue", "type": "string", "default": "1e-3"},
            {"name": "swissprot_evalue", "type": "string", "default": "1e-3"},
            {"name": "pfam_evalue", "type": "string", "default": "1e-3"},
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},  # NR比对e值
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},  # KEGG注释使用的e值
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-5},  # Swissprot比对使用的e值
            {"name": "cog_blast_evalue", "type": "float", "default": 1e-5},  # COG比对使用的e值
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},  # PFAM比对使用的e值
            {"name": "nr_blast_identity", "type": "float", "default": 0},  # NR比对identity
            {"name": "kegg_blast_identity", "type": "float", "default": 0},  # KEGG注释使用identity
            {"name": "swissprot_blast_identity", "type": "float", "default": 0},  # Swissprot比对使用identity
            {"name": "cog_blast_identity", "type": "float", "default": 0},  # COG比对使用identity
            {"name": "nr_blast_similarity", "type": "float", "default": 0},  # NR比对similarity
            {"name": "kegg_blast_similarity", "type": "float", "default": 0},  # KEGG注释使用similarity
            {"name": "swissprot_blast_similarity", "type": "float", "default": 0},  # Swissprot比对使用similarity
            {"name": "cog_blast_similarity", "type": "float", "default": 0},  # COG比对使用similarity
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg,swissprot,pfam'},
            # 表达量分析
            {"name": "express_method", "type": "string", "default": "rsem"},# 表达量分析手段: Salmon, Kallisto, RSEM
            {"name": "exp_way", "type": "string", "default": "fpkm"}, # fpkm or tpm
            # 表达差异分析
            {"name": "diff_method", "type": "string", "default": "DESeq2"},# 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "diff_fdr_ci", "type": "float", "default": 0.05},  # 显著性水平
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  #选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  #Bonferroni,Holm,BH,BY
            # SNP/INDEL分析
            {"name": "is_snp", "type": "string", "default": "True"}, #是否进行snp分析,True, False, Skip
            {"name": "snp_method", "type": "string", "default": "samtools"}, #samtools or GATK
            {"name": "qual", "type": "float", "default": 20},  #过滤低质量的SNP位点
            {"name": "dp", "type": "int", "default": 1},  #过滤低质量的SNP位点
            # 可变剪切分析
            {"name": "is_as", "type": "string", "default": "True"}, #是否进行as分析,True, False, Skip
            {"name": "as_method", "type": "string", "default": "rmats"}, #rmats or asprofile
            {"name": "filter_tpm", "type": "string", "default": None},
            {"name": "qc_soft", "type": "string", "default": None},
            {"name": "rrna_stop", "type": "string", "default": None},
            {"name": "rrna_ratio", "type": "string", "default": None},
            {"name": "rrna_sample_percent", "type": "string", "default": None},
            {"name": "mapping_ratio", "type": "string", "default": None},
            {"name": "mapping_sample_percent", "type": "string", "default": None},
            {"name": "mapping_stop", "type": "string", "default": None},
            {"name": "genome_id", "type": "string", "default": None},
            {"name": "datatype", "type": "string", "default": None},
        ]
        #获取输出目录
        self.workflow_output_tmp = self._sheet.output
        '''
        ## 暂时将上传目录设置为对象存储目录， 如果框架修改、需要进行相应的改动
        region_bucket = Config().get_project_region_bucket(project_type="ref_rna_v2")
        if re.match(r"^\w+://\S+/.+$", workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        elif workflow_output_tmp.startswith("/mnt/ilustre/tsanger-data/"):
            self.workflow_output = self.workflow_output_tmp.replace(region_bucket, '/mnt/ilustre/tsanger-data/')
        '''
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")


        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())

        #获取数据库中该物种相关文件路径
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        col = db["sg_genome_db"]
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        try:
            genome_info = col.find_one({"name" : self.option("ref_genome"), "assembly" : self.option("genome_version"), "annot_version" : self.option("genome_annot_version")})
        except:
            if not genome_info:
                self.set_error("数据库中不存在该物种注释信息，程序退出", code = "13700320")
        try:
            self.ref_annot_dir = os.path.join(db_path, genome_info["anno_path_v2"])
        except:
            if not self.ref_annot_dir:
                self.set_error("数据库中不存在该物种参考注释文件，程序退出", code = "13700321")
        try:
            self.ref_genome = os.path.join(db_path, genome_info["dna_fa"])
        except:
            if not self.ref_genome:
                self.set_error("数据库中不存在该物种参考基因组文件，程序退出", code = "13700323")
        try:
            self.ref_gtf = os.path.join(db_path, genome_info["gtf"])
        except:
            if not self.ref_gtf:
                self.set_error("数据库中不存在该物种参考基因组注释GTF文件，程序退出", code = "13700324")
        try:
            self.genome_id = genome_info["genome_id"]
        except:
            if not self.genome_id:
                self.set_error("数据库中不存在该物种参考基因组注释GTF文件，程序退出", code = "13700325")
        try:
            self.des = os.path.join(db_path, genome_info["bio_mart_annot"])
            self.des_type = genome_info["biomart_gene_annotype"]
        except:
            if not (self.des or self.des_type):
                self.set_error("数据库中不存在该物种的功能描述信息", code = "13700326")
        try:
            self.known_cds = os.path.join(db_path, genome_info["cds"])
            self.known_pep = os.path.join(db_path, genome_info["pep"])
        except:
            if not (self.known_cds or self.known_pep):
                self.set_error("数据库中不存在该物种的CDS或蛋白序列信息", code = "13700327")
        try:
            self.entrez = os.path.join(db_path, genome_info["ensemble2entrez"])
        except:
            if not self.entrez:
                self.set_error("数据库中不存在该物种entrez注释信息", code = "13700328")
        try:
            self.genome_stat = os.path.join(db_path, genome_info["gene_stat"])
        except:
            if not self.genome_stat:
                self.set_error("数据库中不存在该文件：genome_stat.xls", code = "13700329")
        try:
            self.g2t2p = os.path.join(db_path, genome_info['g2t2p'])
        except:
            if not self.g2t2p:
                self.set_error('g2t2p does not exist in the database')
        try:
            self.species_name = genome_info["name"]
            self.species = genome_info["taxon_id"]
            self.ref_anno_version = genome_info["assembly"]
            self.hyperlink = genome_info["ensemble_web"]
            self.known_ko = genome_info['kegg']
        except:
            if not (self.species_name or self.species or self.ref_anno_version or self.hyperlink):
                self.set_error("数据库中该物种注释信息不全", code = "13700330")

        '''
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.json_dict = self.get_json()
        if "anno_path_v2" not in self.json_dict[self.option("ref_genome")]:
            self.set_error("json文件中不存在参考注释文件，程序退出", code = "13700301")
        self.ref_annot_dir = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["anno_path_v2"])
        if "dna_fa" not in self.json_dict[self.option("ref_genome")]:
            self.set_error("json文件中不存在参考基因组文件，程序退出", code = "13700302")
        self.ref_genome = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["dna_fa"])
        if not os.path.exists(self.ref_genome):
            self.set_error("不存在参考基因组文件，程序退出", code = "13700303")
        if "taxon_id" not in self.json_dict[self.option("ref_genome")]:
            self.set_error("json文件中不存在taxon_id，程序退出", code = "13700304")
        self.taxon_id = self.json_dict[self.option("ref_genome")]["taxon_id"]
        if "anno_path" not in self.json_dict[self.option("ref_genome")]:
            self.set_error("json文件中不存在注释文件，程序退出", code = "13700305")
        self.anno_path = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["anno_path"])
        if not os.path.exists(self.anno_path):
            self.set_error("不存在注释文件，程序退出", code = "13700306")
        self.logger.info("注释文件路径为： " + self.anno_path)
        if "gtf" not in self.json_dict[self.option("ref_genome")]:
            self.set_error("json文件中不存在参考基因组注释GTF文件，程序退出", code = "13700307")
        self.ref_gtf = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["gtf"])
        if not os.path.exists(self.ref_gtf):
            self.set_error("不存在参考基因组注释GTF文件，程序退出", code = "13700308")
        '''


        #添加tool/module
        self.filecheck = self.add_tool("ref_rna_v2.file_check")
        if self.option('fastp'):
            self.qc = self.add_module("datasplit.fastp_rna")
        else:
            self.qc = self.add_module("ref_rna_v2.hiseq_qc")
        self.qc_stat_before = self.add_module("ref_rna_v2.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("ref_rna_v2.hiseq_reads_stat")
        self.mapping = self.add_module("ref_rna_v2.rnaseq_mapping")
        self.annot_filter_ref = self.add_module("ref_rna_v2.annot_filter")
        self.annot_class_ref = self.add_module("ref_rna_v2.annot_class")
        if self.option("is_assemble") == True:
            self.assembly = self.add_module("ref_rna_v2.refrna_assemble")
            self.annot_mapdb = self.add_module("ref_rna_v2.annot_mapdb")
            self.annot_class_new = self.add_module("ref_rna_v2.annot_class")
            self.annot_orfpfam = self.add_module("ref_rna_v2.annot_orfpfam")
            self.annot_filter_new = self.add_module("ref_rna_v2.annot_filter")
        else:
            self.transcripts = self.add_tool("ref_rna_v2.transcript_abstract")
        self.merge_annot =  self.add_tool("ref_rna_v2.annotation.merge_annot")
        if self.option("map_assess_method") != None:
            self.map_qc = self.add_module("ref_rna_v2.map_assessment")
        if self.option("is_snp").lower() == "true":
            if self.option("snp_method").lower() == "gatk":
                self.snp = self.add_module("ref_rna_v2.snp_rna")
            else:
                self.snp = self.add_module("ref_rna_v2.sam_rna")
        if self.option("is_as").lower() == "true":
            self.altersplicing = self.add_module("ref_rna_v2.rmats")
        self.express = self.add_module("ref_rna_v2.quant")
        self.exp_pca = self.add_tool("ref_rna_v2.exp_pca")
        self.exp_corr = self.add_tool("ref_rna_v2.exp_corr")
        self.diffexpress = self.add_tool("ref_rna_v2.diffexp")
        self.gene_fa = self.add_tool("ref_rna_v2.gene_fa")

        #判断流程结束tool/module list
        if self.option("sample_num") == "multiple":
            self.final_tools = [self.annot_filter_ref, self.annot_class_ref, self.express, self.diffexpress, self.merge_annot, self.gene_fa]
            if self.option("is_snp").lower() == "true":
                self.final_tools.append(self.snp)
            if self.option("group_table").prop["sample_number"] > 2:
                self.final_tools.append(self.exp_pca)
            if self.option("is_as").lower() == "true":
                self.final_tools.append(self.altersplicing)
            self.logger.info(self.final_tools)
        else:
            self.final_tools = [self.annot_filter_ref, self.annot_class_ref, self.express, self.merge_annot, self.gene_fa]
        if self.option("map_assess_method") != None:
            self.final_tools.append(self.map_qc)
        if self.option("is_assemble") == True:
            self.final_tools.append(self.assembly)
            self.final_tools.append(self.annot_mapdb)
            self.final_tools.append(self.annot_class_new)
            self.final_tools.append(self.annot_orfpfam)
            self.final_tools.append(self.annot_filter_new)

        # 添加step，显示在页面进度条
        if self.option("sample_num") == "multiple":
            if self.option("group_table").prop["sample_number"] > 2:
                self.step.add_steps("filecheck", "rna_qc", "mapping", "express", "diffexpress", "snp_rna", "assembly",
                                "map_qc", "altersplicing", "annot_mapdb", "annot_orfpfam", "merge_annot", "exp_pca",
                                "exp_corr", "qc_stat_before", "qc_stat_after", "annot_class_new", "annot_class_ref",
                                "annot_filter_new", "annot_filter_ref", "transcripts_fa", "gene_fa")
            else:
                self.step.add_steps("filecheck", "rna_qc", "mapping", "express", "diffexpress", "snp_rna", "assembly",
                                    "map_qc", "altersplicing", "annot_mapdb", "annot_orfpfam", "merge_annot","exp_corr",
                                    "qc_stat_before", "qc_stat_after", "annot_class_new", "annot_class_ref",
                                    "annot_filter_new", "annot_filter_ref", "transcripts_fa", "gene_fa")
        else:
            self.step.add_steps("filecheck", "rna_qc", "mapping", "express", "assembly",
                                "map_qc", "annot_mapdb", "annot_orfpfam", "merge_annot",
                                "qc_stat_before", "qc_stat_after", "annot_class_new", "annot_class_ref",
                                "annot_filter_new", "annot_filter_ref", "transcripts_fa", "gene_fa")

    def check_options(self):
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE", code = "13700309")
        try:
            nr_evalue = float(self.option("nr_evalue"))
            cog_evalue = float(self.option("cog_evalue"))
            kegg_evalue = float(self.option("kegg_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范", code = "13700310")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("cog_blast_evalue", cog_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围", code = "13700311")
        if not self.option("cog_blast_evalue") > 0 and not self.option("cog_blast_evalue") < 1:
            raise OptionError("Cog比对的E值超出范围", code = "13700312")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围", code = "13700313")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围", code = "13700314")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围", code = "13700315")
        if not self.option("align_method").lower() in ["tophat", "hisat"]:
            raise OptionError("比对软件应在Tophat与Hisat中选择", code = "13700316")
        for i in self.option('map_assess_method').split(','):
            if i.lower() not in ["saturation", "distribution", "coverage", "chr_stat"]:
                raise OptionError("比对质量评估分析没有%s，请检查", variables = (i), code = "13700317")
        if self.option("is_assemble") == True:
            if self.option("assemble_method").lower() not in ["cufflinks", "stringtie"]:
                raise OptionError("拼接软件应在cufflinks和stringtie中选择", code = "13700318")
        if self.option("kegg_database") == "All":
            self.option("kegg_database","All")
        elif self.option("kegg_database") == "Animal":
            self.option("kegg_database","Animals")
        elif self.option("kegg_database") == "Plant":
            self.option("kegg_database","Plants")
        elif self.option("kegg_database") == "Protist":
            self.option("kegg_database","Protists")
        if self.option("nr_database") == "All":
            self.option("nr_database", "nr")
        elif self.option("nr_database") == "Animal":
            self.option("nr_database", "metazoa")
        elif self.option("nr_database") == "Plant":
            self.option("nr_database", "viridiplantae")
        else:
            nr = self.option("nr_database").lower()
            self.option("nr_database", nr)
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_json(self):
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    def run(self):
        """
        ref-rna workflow run方法
        :return:
        """
        if self.option("sample_num") == "multiple":
            self.filecheck.on('end', self.run_qc)
            self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
            self.filecheck.on('end', self.run_annot_filter_ref)
            self.annot_filter_ref.on('end', self.run_annot_class_ref)
            self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
            self.qc.on('end', self.run_mapping)
            if self.option("map_assess_method") != None:
                self.mapping.on('end', self.run_map_assess)
            if self.option("is_assemble") == True:
                self.mapping.on('end', self.run_assembly)
            if self.option("is_snp").lower() == "true":
                self.mapping.on("end", self.run_snp)
            if self.option("is_assemble") == True:
                self.assembly.on('end', self.run_annot_mapdb)
                self.assembly.on('end', self.run_annot_orfpfam)
                self.on_rely([self.annot_mapdb,self.annot_orfpfam], self.run_annot_filter_new)
                self.annot_filter_new.on('end', self.run_annot_class_new)
                self.on_rely([self.annot_class_new,self.annot_class_ref], self.run_merge_annot)
                self.assembly.on('end', self.run_express)
                self.assembly.on('end', self.run_gene_fa)
            else:
                self.mapping.on('end', self.run_transcripts)
                self.transcripts.on('end', self.run_express)
                self.annot_class_ref.on('end', self.run_merge_annot)
                self.filecheck.on('end', self.run_gene_fa)
            if self.option("is_as").lower() == "true":
                self.mapping.on('end', self.run_altersplicing)
            self.express.on('end', self.run_diffexpress)
            self.diffexpress.on('end', self.run_exp_corr)
            if self.option("group_table").prop["sample_number"] > 2:
                self.diffexpress.on('end', self.run_exp_pca)
        else:
            self.filecheck.on('end', self.run_qc)
            self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
            self.filecheck.on('end', self.run_annot_filter_ref)
            self.annot_filter_ref.on('end', self.run_annot_class_ref)
            self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
            self.qc.on('end', self.run_mapping)
            if self.option("map_assess_method") != None:
                self.mapping.on('end', self.run_map_assess)
            if self.option("is_assemble") == True:
                self.mapping.on('end', self.run_assembly)
            if self.option("is_assemble") == True:
                self.assembly.on('end', self.run_annot_mapdb)
                self.assembly.on('end', self.run_annot_orfpfam)
                self.on_rely([self.annot_mapdb,self.annot_orfpfam], self.run_annot_filter_new)
                self.annot_filter_new.on('end', self.run_annot_class_new)
                self.on_rely([self.annot_class_new,self.annot_class_ref], self.run_merge_annot)
                self.assembly.on('end', self.run_express)
                self.assembly.on('end', self.run_gene_fa)
            else:
                self.mapping.on('end', self.run_transcripts)
                self.transcripts.on('end', self.run_express)
                self.annot_class_ref.on('end', self.run_merge_annot)
                self.filecheck.on('end', self.run_gene_fa)
        self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(RefrnaWorkflow, self).run()

    def run_gene_fa(self):
        if self.option("is_assemble") == True:
            opts = {
                "ref_new_gtf": self.assembly.option("ref_and_new_gtf").prop['path'],  # ref和new合并后的gtf
                "ref_genome_custom": self.ref_genome
            }
        else:
            opts = {
                "ref_new_gtf": self.ref_gtf,
                "ref_genome_custom": self.ref_genome
            }
        self.gene_fa.set_options(opts)
        self.gene_fa.on('end', self.set_output, 'gene_fa')
        self.gene_fa.on('start', self.set_step, {'start': self.step.gene_fa})
        self.gene_fa.on('end', self.set_step, {'end': self.step.gene_fa})
        self.gene_fa.run()
        pass

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            "in_gtf": self.ref_gtf,
            'sample_num': self.option('sample_num'),
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
        self.logger.info("开始运行质控")
        if self.option('fastp'):
            with open(self.option('fastq_dir').prop['path'] + '/list.txt', 'r') as list_r, open(self.option('fastq_dir').prop['path'] + '/tmp.txt', 'w') as tmp_w:
                for line in list_r:
                    tmp_w.write(self.option('fastq_dir').prop['path'] + '/' + line)
            self.qc.set_options({
                'sample_path': self.option('fastq_dir').prop['path'] + '/tmp.txt',
                'length_required': '30'
            })
        else:
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
            "ref_genome": self.option("ref_genome"),
            "genome_version": self.option("genome_version"), # 参考基因组版本
            "genome_annot_version": self.option("genome_annot_version"),  # 参考基因组注释版本
            "mapping_method": self.option("align_method").lower(),  # 比对软件
            "seq_method": self.option("fq_type"),   # PE or SE
            "fastq_dir": self.qc.option("sickle_dir"),
            "assemble_method": self.option("assemble_method").lower(),
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

    def run_assembly(self):
        self.logger.info("开始运行拼接步骤")
        opts = {
            "sample_bam_dir": self.mapping.option("bam_output"),
            "assemble_method": self.option("assemble_method").lower(),
            "ref_gtf": self.filecheck.option("gtf"),
            "ref_fa": self.ref_genome,
        }
        # 如果具有链特异性
        if self.option("strand_specific"):
            if self.option("strand_dir") == "forward":
                strand_dir = "firststrand"
            else:
                strand_dir = "secondstrand"
            opts.update({
                "strand_direct": strand_dir,
                "fr_stranded": "fr-stranded"
                })
        else:
            opts.update({
                "fr_stranded": "fr-unstranded"
                })
        self.assembly.set_options(opts)
        self.assembly.on("end", self.set_output, "assembly")
        self.assembly.on('start', self.set_step, {'start': self.step.assembly})
        self.assembly.on('end', self.set_step, {'end': self.step.assembly})
        self.assembly.run()

    def run_annot_filter_ref(self):
        options = {
            "blast_nr_xml" : self.ref_annot_dir + "/annot_mapdb/nr/blast.xml",
            "blast_eggnog_xml" : self.ref_annot_dir + "/annot_mapdb/eggnog/blast.xml",
            "blast_kegg_xml": self.ref_annot_dir + "/annot_mapdb/kegg/blast.xml",
            "blast_swissprot_xml" : self.ref_annot_dir + "/annot_mapdb/swissprot/blast.xml",
            "pfam_domain" : self.ref_annot_dir + "/annot_orfpfam/pfam_domain",
            "blast2go_annot" : self.ref_annot_dir + "/annot_mapdb/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_blast_evalue'),
            'nr_identity': self.option('nr_blast_identity'),
            'nr_similarity': self.option('nr_blast_similarity'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'swissprot_identity': self.option('swissprot_blast_identity'),
            'swissprot_similarity': self.option('swissprot_blast_similarity'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'eggnog_identity': self.option('cog_blast_identity'),
            'eggnog_similarity': self.option('cog_blast_similarity'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'kegg_identity': self.option('kegg_blast_identity'),
            'kegg_similarity': self.option('kegg_blast_similarity'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        }
        self.annot_filter_ref.set_options(options)
        self.annot_filter_ref.on('start', self.set_step, {'start': self.step.annot_filter_ref})
        self.annot_filter_ref.on('end', self.set_step, {'end': self.step.annot_filter_ref})
        self.annot_filter_ref.on('end', self.set_output, "annot_filter_ref")
        self.annot_filter_ref.run()

    def run_annot_filter_new(self):
        options = {
            "blast_nr_xml" : self.annot_mapdb.output_dir + "/nr/blast.xml",
            "blast_eggnog_xml" : self.annot_mapdb.output_dir + "/eggnog/blast.xml",
            "blast_kegg_xml": self.annot_mapdb.output_dir + "/kegg/blast.xml",
            "blast_swissprot_xml" : self.annot_mapdb.output_dir + "/swissprot/blast.xml",
            "pfam_domain" : self.annot_orfpfam.output_dir + "/pfam_domain",
            "blast2go_annot" : self.annot_mapdb.output_dir + "/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_blast_evalue'),
            'nr_identity': self.option('nr_blast_identity'),
            'nr_similarity': self.option('nr_blast_similarity'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'swissprot_identity': self.option('swissprot_blast_identity'),
            'swissprot_similarity': self.option('swissprot_blast_similarity'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'eggnog_identity': self.option('cog_blast_identity'),
            'eggnog_similarity': self.option('cog_blast_similarity'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'kegg_identity': self.option('kegg_blast_identity'),
            'kegg_similarity': self.option('kegg_blast_similarity'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        }
        self.annot_filter_new.set_options(options)
        self.annot_filter_new.on('start', self.set_step, {'start': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_step, {'end': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_output, "annot_filter_new")
        self.annot_filter_new.run()

    def run_annot_class_ref(self):
        ref_filter_dir = self.annot_filter_ref.output_dir
        options = {
            'taxonomy': self.option("kegg_database"),
            'blast_nr_xml': ref_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml': ref_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': ref_filter_dir + "/pfam/pfam_domain.filter.xls",
            'blast2go_annot': ref_filter_dir + "/go/blast2go_merge.xls.filter.xls",
            'g2t2p': self.g2t2p,
            'gtf': self.ref_gtf,
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'type': "ref",
        }

        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        if self.known_ko != "" and os.path.exists(os.path.join(db_path, self.known_ko)):
            options.update({
                "known_ko": os.path.join(db_path, self.known_ko)
            })

        self.annot_class_ref.set_options(options)
        self.annot_class_ref.on('start', self.set_step, {'start': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_step, {'end': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_output, "annot_class_ref")
        self.annot_class_ref.run()

    def run_annot_class_new(self):
        new_filter_dir = self.annot_filter_new.output_dir
        options = {
            'taxonomy': self.option("kegg_database"),
            'blast_nr_xml': new_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml': new_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': new_filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': new_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': new_filter_dir + "/pfam/pfam_domain.filter.xls",
            'blast2go_annot': new_filter_dir + "/go/blast2go_merge.xls.filter.xls",
            'gene2trans': self.annot_orfpfam.output_dir + "/all_tran2gen.txt",
            'gtf': self.assembly.option("new_transcripts_gtf"),
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
        }
        self.annot_class_new.set_options(options)
        self.annot_class_new.on('start', self.set_step, {'start': self.step.annot_class_new})
        self.annot_class_new.on('end', self.set_step, {'end': self.step.annot_class_new})
        self.annot_class_new.on('end', self.set_output, "annot_class_new")
        self.annot_class_new.run()

    def run_transcripts(self):
        if self.option("is_assemble") == True:
            opts = {
                "ref_genome": self.ref_genome,
                "ref_genome_gtf": self.assembly.option("ref_and_new_gtf")
            }
            self.transcripts.set_options(opts)
            self.transcripts.run()
        else:
            opts = {
                "ref_genome": self.ref_genome,
                "ref_genome_gtf": self.ref_gtf
            }
            self.transcripts.set_options(opts)
            self.transcripts.on("end", self.set_output, "transcripts")
            self.transcripts.on('start', self.set_step, {'start': self.step.transcripts_fa})
            self.transcripts.on('end', self.set_step, {'end': self.step.transcripts_fa})
            self.transcripts.run()

    def run_express(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option('strand_dir') == 'forward':
                if self.option("fq_type") == "PE":
                    self.libtype = "rf"
                else:
                    self.libtype = "r"
            else:
                if self.option("fq_type") == "PE":
                    self.libtype = "fr"
                else:
                    self.libtype = "f"
        else:
            self.libtype = None
        if self.option("is_assemble") == True:
            opts = {
                "fastq" : self.qc.option("fq_list"),
                "method" : self.option("express_method"),
                "libtype" : self.libtype,
                "transcriptome" : self.assembly.option("all_transcripts_fa"),
                "t2g" : self.assembly.option("trans2gene").prop["path"],
            }
        else:
            opts = {
                "fastq" : self.qc.option("fq_list"),
                "method" : self.option("express_method"),
                "libtype" : self.libtype,
                "transcriptome" : self.transcripts.option("trans_fa"),
                "t2g" : self.transcripts.option("trans2gene").prop["path"],
            }
        self.express.set_options(opts)
        self.express.on("end", self.set_output, "express")
        self.express.on('start', self.set_step, {'start': self.step.express})
        self.express.on('end', self.set_step, {'end': self.step.express})
        self.express.run()

    def run_exp_pca(self):
        self.logger.info("开始运行pca")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : self.work_dir + "/gene.fpkm.matrix"
            }
        else:
            opts = {
                "exp" : self.work_dir + "/gene.tpm.matrix"
            }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.run()

    def run_exp_corr(self):
        self.logger.info("开始运行聚类分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : self.work_dir + "/gene.fpkm.matrix"
            }
        else:
            opts = {
                "exp" : self.work_dir + "/gene.tpm.matrix"
            }
        self.exp_corr.set_options(opts)
        self.exp_corr.on("end", self.set_output, "exp_corr")
        self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.run()

    def run_diffexpress(self):
        self.logger.info("开始运行基因差异表达分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            count_file = self.work_dir + "/gene.count.matrix"
            fpkm_file = self.work_dir + "/gene.fpkm.matrix"
            trans_fpkm_file = self.work_dir + "/transcript.fpkm.matrix"
            with open(self.express.option("gene_count").prop['path'], "r") as f1, open(count_file, "w") as w1:
                for line in f1:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w1.write(line)
            with open(self.express.output_dir + "/gene.fpkm.matrix", "r") as f2, open(fpkm_file, "w") as w2:
                for line in f2:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w2.write(line)
            with open(self.express.output_dir + "/transcript.fpkm.matrix", "r") as f2, open(trans_fpkm_file, "w") as w2:
                for line in f2:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w2.write(line)
            opts = {
                "count" : count_file,
                "exp" : fpkm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                "pvalue_padjust" : self.option("pvalue_padjust"),
                "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
        else:
            count_file = self.work_dir + "/gene.count.matrix"
            tpm_file = self.work_dir + "/gene.tpm.matrix"
            trans_tpm_file = self.work_dir + "/transcript.tpm.matrix"
            with open(self.express.option("gene_count").prop['path'], "r") as f1, open(count_file, "w") as w1:
                for line in f1:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w1.write(line)
            with open(self.express.option("gene_tpm").prop['path'], "r") as f2, open(tpm_file, "w") as w2:
                for line in f2:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w2.write(line)
            with open(self.express.output_dir + "/transcript.tpm.matrix", "r") as f2, open(trans_tpm_file, "w") as w2:
                for line in f2:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        w2.write(line)
            opts = {
                "count" : count_file,
                "exp" : tpm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                "pvalue_padjust" : self.option("pvalue_padjust"),
                "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

    def run_annot_mapdb(self, event):
        self.logger.info("开始运行diamond注释")
        opts = {
            "query" : self.assembly.option("new_transcripts_fa"),
            "nr_db" : self.option("nr_database")
        }
        self.annot_mapdb.set_options(opts)
        self.annot_mapdb.on("end", self.set_output, "annot_mapdb")
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.logger.info("开始运行pfam注释")
        opts = {
            "fasta" : self.assembly.option("new_transcripts_fa"),
            "gtf" : self.assembly.option("new_transcripts_gtf")
        }
        self.annot_orfpfam.set_options(opts)
        self.annot_orfpfam.on("end", self.set_output, "annot_orfpfam")
        self.annot_orfpfam.on('start', self.set_step, {'start': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_step, {'end': self.step.annot_orfpfam})
        self.annot_orfpfam.run()

    def run_merge_annot(self):
        options = {
            "annot_class_ref":self.annot_class_ref.output_dir,
            "is_assemble": self.option("is_assemble")
        }
        if self.option("is_assemble") == True:
            options.update({
                "annot_class": self.annot_class_new.output_dir,
                "annot_db": self.annot_mapdb.output_dir
            })
        self.merge_annot.set_options(options)
        self.merge_annot.on('start', self.set_step, {'start': self.step.merge_annot})
        self.merge_annot.on('end', self.set_step, {'end': self.step.merge_annot})
        self.merge_annot.on('end', self.set_output, "merge_annot")
        self.merge_annot.run()

    def run_snp(self):
        self.logger.info("开始运行snp步骤")
        if self.option("snp_method").lower() == "gatk":
            opts = {
                "ref_genome_custom": self.ref_genome,
                "ref_genome":  "customer_mode",
                "ref_gtf": self.ref_gtf,
                "in_bam": self.mapping.option("bam_output"),
                'des': self.des,
                'des_type': self.des_type,
            }
        else:
            opts = {
                "ref_genome_custom": self.ref_genome,
                "ref_genome":  "customer_mode",
                "ref_gtf": self.ref_gtf,
                "fq_list": self.qc.option("fq_list").prop["path"],
                "bamlist": self.mapping.option("bamlist"),
                'des': self.des,
                'des_type': self.des_type,
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
            "analysis": self.option("map_assess_method")
        }
        self.map_qc.set_options(opts)
        self.map_qc.on("start", self.set_step, {"start": self.step.map_qc})
        self.map_qc.on("end", self.set_step, {"end": self.step.map_qc})
        self.map_qc.on("end", self.set_output, "map_qc")
        self.map_qc.run()

    def run_altersplicing(self):
        '''
        forward ~ RF ~ fr-firststrand
        reverse ~ FR ~ fr-secondstrand
        '''
        if self.option('strand_specific'):
            if option('strand_dir') == 'forward':
                lib_type = 'fr-firststrand'
            elif option('strand_dir') == 'reverse':
                lib_type = 'fr-secondstrand'
        else:
            lib_type = 'fr-unstranded'
        if self.option('fq_type') == 'PE':
            seq_type = 'paired'
        elif self.option('fq_type') == 'SE':
            seq_type = 'single'
        opts = {
            'sample_bam_dir': self.mapping.option('bam_output'),
            'lib_type': lib_type,
            'seq_type': seq_type,
            'ref_gtf': self.filecheck.option('gtf'),
            'group_table': self.option('group_table'),
            'rmats_control': self.option('control_file')
        }
        self.altersplicing.set_options(opts)
        self.altersplicing.on('start', self.set_step, {'start': self.step.altersplicing})
        self.altersplicing.on('end', self.set_step, {'end': self.step.altersplicing})
        self.altersplicing.on('end', self.set_output, 'altersplicing')
        self.altersplicing.run()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "13700319")
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
        if event['data'] == 'map_qc':
            self.move2outputdir(obj.output_dir, 'map_qc')
        if event['data'] == 'assembly':
            self.move2outputdir(obj.output_dir, 'assembly')
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'snp')
        if event['data'] == 'express':
            self.move2outputdir(obj.output_dir, 'express')
        if event['data'] == 'exp_pca':
            self.move2outputdir(obj.output_dir, 'exp_pca')
        if event['data'] == 'exp_corr':
            self.move2outputdir(obj.output_dir, 'exp_corr')
        if event['data'] == 'diffexpress':
            self.move2outputdir(obj.output_dir, 'diffexpress')
        if event['data'] == 'altersplicing':
            self.move2outputdir(obj.output_dir, 'altersplicing')
        if event['data'] == 'annot_mapdb':
            self.move2outputdir(obj.output_dir, 'annot_mapdb')
        if event['data'] == 'annot_orfpfam':
            self.move2outputdir(obj.output_dir, 'annot_orfpfam')
        if event['data'] == 'merge_annot':
            self.move2outputdir(obj.output_dir, 'merge_annot')
        if event['data'] == 'transcripts':
            self.move2outputdir(obj.output_dir, 'abstract_transcripts')
        if event['data'] == 'annot_class_new':
            self.move2outputdir(obj.output_dir, 'annot_class_new')
        if event['data'] == 'annot_class_new':
            self.move2outputdir(obj.output_dir, 'annot_class_new')
        if event['data'] == 'annot_class_ref':
            self.move2outputdir(obj.output_dir, 'annot_class_ref')
        if event['data'] == 'annot_filter_new':
            self.move2outputdir(obj.output_dir, 'annot_filter_new')
        if event['data'] == 'annot_filter_ref':
            self.move2outputdir(obj.output_dir, 'annot_filter_ref')

    def end(self):
        self.run_api()
        ## 更新一系列主表的字段，用于页面交互分析
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        col = db["sg_task"]
        col.update({"task_id" : self.task_id}, {"$set": {"refrna_seqdb": self.workflow_output + "/Sequence_database/refrna_seqs.db"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"ref_gtf": self.ref_gtf}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"ref_genome": self.ref_genome}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"genome_id": self.genome_id}}, upsert=True)
        if self.option("is_assemble") == True:
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_t2g": self.workflow_output + "/Assemble/trans2gene.txt"}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.workflow_output + "/Assemble/all_transcripts.fa"}}, upsert=True)
        else:
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_t2g": self.workflow_output + "/Transcripts/trans2gene.txt"}}, upsert=True)
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.workflow_output + "/Transcripts/all_transcripts.fa"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"fastq": self.workflow_output + "/QC/fq_list.txt"}}, upsert=True)
        # begin for rmats
        col.update({"task_id" : self.task_id}, {"$set": {'as_gtf': self.assembly.option('ref_and_new_gtf').path}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {'fq_type': self.option('fq_type')}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {'strand_specific': self.option('strand_specific')}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {'strand_dir': self.option('strand_dir')}}, upsert=True)
        # final for rmats
        col1 = db["sg_annotation_stat"]
        col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_exp"]
        col2.update({"task_id" : self.task_id, "exp_level" : "G"}, {"$set": {"count_file": self.workflow_output + "/Express/ExpAnnalysis/gene.count.matrix.xls"}}, upsert=True)
        col2.update({"task_id" : self.task_id, "exp_level" : "T"}, {"$set": {"count_file": self.workflow_output + "/Express/ExpAnnalysis/transcript.count.matrix.xls"}}, upsert=True)
        if self.option("sample_num") == "multiple":
            if self.option("is_as").lower() == "true":
                col3 = db["sg_splicing_rmats"]
                for file in os.listdir(self.altersplicing.output_dir):
                    tmp_group_list = file.split("_vs_")
                    group_a = tmp_group_list[0]
                    group_b = tmp_group_list[1]
                    compare_plan = group_a + "|" + group_b
                    col3.update({"task_id" : self.task_id, "compare_plan": compare_plan}, {"$set": {"result_dir": self.workflow_output + "/AS/" + file}}, upsert=True)
        self.merge_annotation_exp_matrix() # 表达量表增加注释信息
        if self.option("sample_num") == "multiple":
            self.merge_annotation_diffexp_matrix() # 差异表达量表增加注释信息
        self.modify_output()

        super(RefrnaWorkflow, self).end()

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
        seq_db = origin_dir + "/Sequence_database/refrna_seqs.db"
        os.link(seq_db, target_dir + "/Sequence_database/refrna_seqs.db")
        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        fq_list = origin_dir + "/QC_stat/sickle_dir/fq_list.txt"
        os.mkdir(target_dir + "/QC")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        a_b = pd.read_table(target_dir + "/QC/rawdata_statistics.xls",index_col="#Sample_ID")
        b_a = pd.read_table(target_dir + "/QC/cleandata_statistics.xls",index_col="#Sample_ID")
        a_bf =pd.DataFrame(a_b,columns=["Total_Reads","Total_Bases"])
        a_bf.rename(columns={"Total_Reads":"Raw_Reads","Total_Bases":"Raw_Bases"},inplace = True)
        b_a.rename(columns={"Total_Reads":"Clean_Reads","Total_Bases":"Clean_Bases"},inplace = True)
        fq_stat_final = pd.concat([a_bf, b_a], axis=1, join_axes=[b_a.index])
        fq_stat_final.to_csv(target_dir + "/QC/qc_statistics.xls", sep='\t')
        os.link(fq_list, target_dir + "/QC/fq_list.txt")
        ## 质控数据量太大，取消上传至结果目录，会影响页面上的表达定量分析
        # os.mkdir(target_dir + "/QC/cleandata")
        # for file in os.listdir(origin_dir + "/QC_stat/sickle_dir"):
        #     file_path = os.path.join(origin_dir + "/QC_stat/sickle_dir", file)
        #     os.link(file_path, target_dir + "/QC/cleandata/" + file)
        ## modify fq_list.txt
        # with open(target_dir + "/QC/cleandata/fq_list.txt", "r") as f, open(target_dir + "/QC/cleandata/fq_list_tmp.txt", "w") as w:
        #     for line in f:
        #         items = line.strip().split("\t")
        #         if len(items) == 3:
        #             w.write(items[0] + "\t" + self.workflow_output + "/QC/cleandata/" + os.path.basename(items[1]) + "\t" + self.workflow_output + "/QC/cleandata/" + os.path.basename(items[2]) + "\n")
        #         else:
        #             w.write(items[0] + "\t" + self.workflow_output + "/QC/cleandata/" + os.path.basename(items[1]) + "\n")
        # os.rename(target_dir + "/QC/cleandata/fq_list_tmp.txt", target_dir + "/QC/cleandata/fq_list.txt")
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
        #total_stat
        os.link(os.path.join(origin_dir, 'mapping/Comparison_results'), os.path.join(target_dir, 'Align/AlignStat/Comparison_results'))
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
        # Assemble
        if self.option("is_assemble") == True:
            os.makedirs(target_dir + "/Assemble")
            file_path = os.path.join(origin_dir + "/assembly/NewTranscripts/ref_and_new.gtf")
            os.link(file_path, target_dir + "/Assemble/ref_and_new.gtf")
            file_path = os.path.join(origin_dir + "/assembly/NewTranscripts/new_transcripts.fa")
            os.link(file_path, target_dir + "/Assemble/new_transcripts.fa")
            file_path = os.path.join(origin_dir + "/assembly/NewTranscripts/new_transcripts.gtf")
            os.link(file_path, target_dir + "/Assemble/new_transcripts.gtf")
            file_path = os.path.join(origin_dir + "/assembly/NewTranscripts/all_transcripts.fa")
            os.link(file_path, target_dir + "/Assemble/all_transcripts.fa")
            file_path = os.path.join(origin_dir + "/assembly/NewTranscripts/trans2gene")
            os.link(file_path, target_dir + "/Assemble/trans2gene.txt")
        else:
            os.makedirs(target_dir + "/Transcripts")
            file_path = os.path.join(self.transcripts.work_dir + "/trans2gene")
            os.link(file_path, target_dir + "/Transcripts/trans2gene.txt")
            file_path = os.path.join(self.transcripts.work_dir + "/exons.fa")
            os.link(file_path, target_dir + "/Transcripts/all_transcripts.fa")
        #Annotation
        os.makedirs(target_dir + "/Annotation")
        CopyFile().linkdir(origin_dir + "/merge_annot", target_dir + "/Annotation")
        rm_file_dir1 = glob.glob(os.path.join(self.work_dir, "upload_results/Annotation/*/*/*/*.html.mark"))
        rm_file_dir2 = glob.glob(os.path.join(self.work_dir, "upload_results/Annotation/*/*/*/*.KOs.txt"))
        rm_file_dir3 = glob.glob(os.path.join(self.work_dir, "upload_results/Annotation/*/*/*/*.pdf"))
        for file in rm_file_dir1:
            os.remove(file)
        for file in rm_file_dir2:
            os.remove(file)
        for file in rm_file_dir3:
            os.remove(file)
        #Expression
        os.mkdir(target_dir + "/Express")
        os.mkdir(target_dir + "/Express/ExpAnnalysis")
        if self.option("express_method") == "RSEM":
            gene_count = os.path.join(self.express.output_dir + "/gene.count.matrix")
            os.link(gene_count, target_dir + "/Express/ExpAnnalysis/gene.count.matrix.xls")
            gene_tpm = os.path.join(self.express.output_dir + "/gene.tpm.matrix")
            os.link(gene_tpm, target_dir + "/Express/ExpAnnalysis/gene.tpm.matrix.xls")
            gene_fpkm = os.path.join(self.express.output_dir + "/gene.fpkm.matrix")
            os.link(gene_fpkm, target_dir + "/Express/ExpAnnalysis/gene.fpkm.matrix.xls")
            gene_tpm_anno = os.path.join(self.express.output_dir + "/gene.tpm.matrix.annot.xls")
            os.link(gene_tpm_anno, target_dir + "/Express/ExpAnnalysis/gene.tpm.matrix.annot.xls")
            gene_fpkm_anno = os.path.join(self.express.output_dir + "/gene.fpkm.matrix.annot.xls")
            os.link(gene_fpkm_anno, target_dir + "/Express/ExpAnnalysis/gene.fpkm.matrix.annot.xls")
            if self.option("level").lower() == "transcript":
                trans_tpm = os.path.join(self.express.output_dir + "/transcript.tpm.matrix")
                os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
                trans_count = os.path.join(self.express.output_dir + "/transcript.count.matrix")
                os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
                trans_fpkm = os.path.join(self.express.output_dir + "/transcript.fpkm.matrix")
                os.link(trans_fpkm, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.xls")
                trans_tpm_anno = os.path.join(self.express.output_dir + "/transcript.tpm.matrix.annot.xls")
                os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
                trans_fpkm_anno = os.path.join(self.express.output_dir + "/transcript.fpkm.matrix.annot.xls")
                os.link(trans_fpkm_anno, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls")
        else:
            gene_count = os.path.join(self.express.output_dir + "/gene.count.matrix")
            os.link(gene_count, target_dir + "/Express/ExpAnnalysis/gene.count.matrix.xls")
            gene_tpm = os.path.join(self.express.output_dir + "/gene.tpm.matrix")
            os.link(gene_tpm, target_dir + "/Express/ExpAnnalysis/gene.tpm.matrix.xls")
            gene_tpm_anno = os.path.join(self.express.output_dir + "/gene.tpm.matrix.annot.xls")
            os.link(gene_tpm_anno, target_dir + "/Express/ExpAnnalysis/gene.tpm.matrix.annot.xls")
            if self.option("level").lower() == "transcript":
                trans_tpm = os.path.join(self.express.output_dir + "/transcript.tpm.matrix")
                os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
                trans_count = os.path.join(self.express.output_dir + "/transcript.count.matrix")
                os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
                trans_tpm_anno = os.path.join(self.express.output_dir + "/transcript.tpm.matrix.annot.xls")
                os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
        if self.option("sample_num") == "multiple":
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
            if self.option("is_snp").lower() == "true":
                os.makedirs(target_dir + "/SNP")
                CopyFile().linkdir(self.work_dir + "/snp_tmp", target_dir + "/SNP")
                # self.merge1(target_dir + "/SNP/data_anno_pre.xls", target_dir + "/SNP/snp_annotation.xls")
                if os.path.exists(target_dir + "/SNP/snp_annotation.xls"):
                    os.remove(target_dir + "/SNP/snp_annotation.xls")
                if os.path.exists(target_dir + "/SNP/data_anno_pre.xls"):
                    os.remove(target_dir + "/SNP/data_anno_pre.xls")
            # Alternative Splicing
            if self.option("is_as").lower() == "true":
                os.makedirs(target_dir + "/AS")
                CopyFile().linkdir(origin_dir + "/altersplicing", target_dir + "/AS")
                ## 修改A_group_bam.txt和B_group_bam.txt文件中bam的绝对路径
                groups = glob.glob(target_dir + "/AS/*_vs_*")
                for group in groups:
                    A_group_bam = group + "/A_group_bam.txt"
                    B_group_bam = group + "/B_group_bam.txt"
                    with open(A_group_bam, "r") as f1, open(A_group_bam + ".tmp", "w") as w1:
                        bam_str = f1.readline().strip()
                        bam_str_new = ""
                        for each in bam_str.split(","):
                            bam_basename = os.path.basename(each)
                            if bam_str_new == "":
                                bam_str_new = self.workflow_output + "/Align/AlignBam/" + bam_basename
                            else:
                                bam_str_new += "," + self.workflow_output + "/Align/AlignBam/" + bam_basename
                        w1.write(bam_str_new + "\n")
                        os.rename(A_group_bam + ".tmp", A_group_bam)
                    with open(B_group_bam, "r") as f2, open(B_group_bam + ".tmp", "w") as w2:
                        bam_str = f2.readline().strip()
                        bam_str_new = ""
                        for each in bam_str.split(","):
                            bam_basename = os.path.basename(each)
                            if bam_str_new == "":
                                bam_str_new = self.workflow_output + "/Align/AlignBam/" + bam_basename
                            else:
                                bam_str_new += "," + self.workflow_output + "/Align/AlignBam/" + bam_basename
                        w2.write(bam_str_new + "\n")
                        os.rename(B_group_bam + ".tmp", B_group_bam)

        sdir = self.add_upload_dir(target_dir)
        sdir.add_regexp_rules([
            [r"Background/.*\.gtf\.genome_stat\.xls", "xls", " 参考基因组注释信息表", 0, "211001"],
            [r"Align/AlignBam/.*\.bam", "", "样本比对bam文件", 1, "211002"],
            [r"Align/AlignStat/.*_align_stat\.txt", "", "样本比对统计结果表", 0, "211003"],
            [r"Align/QualityAssessment/.*\.chr_distribution\.xls", "", "Reads在不同染色体的分布统计表", 0, "211004"],
            [r"Align/QualityAssessment/.*\.geneBodyCoverage\.txt", "", "基因覆盖度分布结果文件", 0, "211005"],
            [r"Align/QualityAssessment/.*\.region_distribution\.xls", "", "Reads在不同区域的分布统计表", 0, "211006"],
            [r"Align/QualityAssessment/.*\.saturation\.xls", "", "测序饱和度分析结果文件", 0, "211007"],
            [r"DiffExpress/.*_vs_.*\.xls", "xls", "差异分析结果表", 0, "211008"],
            [r"DiffExpress/.*\.DE\.list", "xls", "差异基因列表", 0, "211009"],
            [r"DiffExpress/.*summary\.xls", "xls", "差异统计结果表", 0, "211010"],
            [r"DiffExpress/.*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表", 0, "211011"],
            [r"DiffExpress/.*_vs_.*\..*annot\.xls", "xls", "差异统计注释结果表", 0, "211012"],
            [r"AS/.*/all_events_detail_big_table\.txt", "", "结果详情表", 0, "211013"],
            [r"AS/.*/event_stats\.file\.txt", "", "差异可变剪切事件统计表", 0, "211014"],
            [r"AS/.*/event_type\.file\.txt", "", "可变剪切事件类型统计表", 0, "211015"],
            [r"AS/.*/psi_stats\.file\.txt", "", "差异可变剪切模式变化统计表", 0, "211016"],
            [r"AS/.*/fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '全部可变剪接事件表', 0, "211017"],
            [r"AS/.*/fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '新发现可变剪接事件表', 0, "211018"],
            [r"AS/.*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt", 'txt', '差异事件详情表（JCEC）', 0, "211019"],
            [r"AS/.*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt", 'txt', '差异事件详情表（JC）', 0, "211020"],
            [r"AS/.*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt", 'txt', '事件详情表（JCEC）', 0, "211021"],
            [r"AS/.*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt", 'txt', '事件详情表（JC）', 0, "211022"],
            ])
        sdir.add_relpath_rules([
            [".", "", "流程分析结果目录", 0, "211023"],
            ["Background", "", "项目背景目录", 0, "211024"],
            ["Sequence_database", "", "序列文件数据库", 1, "211025"],
            ["QC", "", "测序数据统计与质控结果", 0, "211026"],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表", 0, "211027"],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表", 0, "211028"],
            ["QC/cleandata", "", "质控结果目录", 1, "211029"],
            ["Align", "", "比对结果目录", 0, "211030"],
            ["Align/AlignBam", "", "比对结果bam目录", 1, "211031"],
            ["Align/AlignStat", "", "比对结果统计目录", 0, "211032"],
            ["Align/QualityAssessment", "", "比对结果整体评估目录", 0, "211033"],
            ["Transcripts", "", "转录本提取结果目录", 0, "211034"],
            ["Transcripts/all_transcripts.fa", "", "所有转录本序列文件", 0, "211035"],
            ["Transcripts/trans2gene.txt", "", "转录本与基因的对应关系文件", 0, "211036"],
            ["Assemble", "", "组装结果目录", 0, "211037"],
            ["Assemble/ref_and_new.gtf", "", "组装结果GTF文件", 0, "211038"],
            ["Assemble/new_transcripts.fa", "", "新预测转录本序列", 0, "211039"],
            ["Assemble/new_transcripts.gtf", "", "新预测转录本GTF文件", 0, "211040"],
            ["Assemble/all_transcripts.fa", "", "所有转录本序列文件", 0, "211041"],
            ["Assemble/trans2gene.txt", "", "转录本与基因的对应关系文件", 0, "211042"],
            ["Annotation", "", "注释结果目录", 0, "211043"],
            ["Annotation/allannot_class", "", "已知基因与新预测合并注释结果目录", 0, "211044"],
            ["Annotation/allannot_class/all_annot.xls", "", "注释结果详情表", 0, "211045"],
            ["Annotation/allannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0, "211046"],
            ["Annotation/allannot_class/all_stat.xls", "", "注释结果统计表", 0, "211047"],
            ["Annotation/allannot_class/cog", "", "cog注释结果目录", 0, "211048"],
            ["Annotation/allannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0, "211049"],
            ["Annotation/allannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0, "211050"],
            ["Annotation/allannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0, "211051"],
            ["Annotation/allannot_class/go/", "", "go注释结果目录", 0, "211052"],
            ["Annotation/allannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0, "211053"],
            ["Annotation/allannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0, "211054"],
            ["Annotation/allannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0, "211055"],
            ["Annotation/allannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0, "211056"],
            ["Annotation/allannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0, "211057"],
            ["Annotation/allannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0, "211058"],
            ["Annotation/allannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0, "211059"],
            ["Annotation/allannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0, "211060"],
            ["Annotation/allannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0, "211061"],
            ["Annotation/allannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0, "211062"],
            ["Annotation/allannot_class/kegg", "", "kegg注释结果目录", 0, "211270"],
            ["Annotation/allannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0, "211271"],
            ["Annotation/allannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0, "211272"],
            ["Annotation/allannot_class/kegg/kegg_pathway_gene_dir:", "", "基因kegg注释通路图", 0, "211273"],
            ["Annotation/allannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0, "211274"],
            ["Annotation/allannot_class/kegg/kegg_pathway_tran_dir:", "", "转录本kegg通路图", 0, "211275"],
            ["Annotation/allannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0, "211276"],
            ["Annotation/allannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0, "211277"],
            ["Annotation/allannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0, "211278"],
            ["Annotation/allannot_class/nr", "", "nr注释结果目录", 0, "211279"],
            ["Annotation/allannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0, "211280"],
            ["Annotation/allannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0, "211281"],
            ["Annotation/allannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0, "211282"],
            ["Annotation/allannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0, "211283"],
            ["Annotation/allannot_class/pfam", "", "pfam注释结果目录", 0, "211284"],
            ["Annotation/allannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0, "211285"],
            ["Annotation/allannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0, "211286"],
            ["Annotation/allannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0, "211287"],
            ["Annotation/allannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0, "211288"],
            ["Annotation/allannot_class/swissprot", "", "swissprot注释结果目录", 0, "211289"],
            ["Annotation/allannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0, "211290"],
            ["Annotation/allannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0, "211291"],
            ["Annotation/allannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0, "211292"],
            ["Annotation/allannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0, "211293"],
            ["Annotation/newannot_class", "", "新预测转录本注释结果目录", 0, "211294"],
            ["Annotation/newannot_class/all_annot.xls", "", "注释结果详情表", 0, "211295"],
            ["Annotation/newannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0, "211296"],
            ["Annotation/newannot_class/all_stat.xls", "", "注释结果统计表", 0, "211297"],
            ["Annotation/newannot_class/cog", "", "cog注释结果目录", 0, "211298"],
            ["Annotation/newannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0, "211299"],
            ["Annotation/newannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0, "211300"],
            ["Annotation/newannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0, "211301"],
            ["Annotation/newannot_class/go/", "", "go注释结果目录", 0, "211302"],
            ["Annotation/newannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0, "211303"],
            ["Annotation/newannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0, "211304"],
            ["Annotation/newannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0, "211305"],
            ["Annotation/newannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0, "211306"],
            ["Annotation/newannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0, "211307"],
            ["Annotation/newannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0, "211308"],
            ["Annotation/newannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0, "211309"],
            ["Annotation/newannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0, "211310"],
            ["Annotation/newannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0, "211311"],
            ["Annotation/newannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0, "211312"],
            ["Annotation/newannot_class/kegg", "", "kegg注释结果目录", 0, "211313"],
            ["Annotation/newannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0, "211314"],
            ["Annotation/newannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0, "211315"],
            ["Annotation/newannot_class/kegg/kegg_pathway_gene_dir:", "", "基因kegg注释通路图", 0, "211316"],
            ["Annotation/newannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0, "211317"],
            ["Annotation/newannot_class/kegg/kegg_pathway_tran_dir:", "", "转录本kegg通路图", 0, "211318"],
            ["Annotation/newannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0, "211319"],
            ["Annotation/newannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0, "211320"],
            ["Annotation/newannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0, "211321"],
            ["Annotation/newannot_class/nr", "", "nr注释结果目录", 0, "211322"],
            ["Annotation/newannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0, "211323"],
            ["Annotation/newannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0, "211324"],
            ["Annotation/newannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0, "211325"],
            ["Annotation/newannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0, "211326"],
            ["Annotation/newannot_class/pfam", "", "pfam注释结果目录", 0, "211327"],
            ["Annotation/newannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0,  "211328"],
            ["Annotation/newannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0,  "211329"],
            ["Annotation/newannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0,  "211330"],
            ["Annotation/newannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0,  "211331"],
            ["Annotation/newannot_class/swissprot", "", "swissprot注释结果目录", 0,  "211332"],
            ["Annotation/newannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0,  "211333"],
            ["Annotation/newannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0,  "211334"],
            ["Annotation/newannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0,  "211335"],
            ["Annotation/newannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0,  "211336"],
            ["Annotation/refannot_class", "", "已知基因/转录本注释结果目录", 0,  "211337"],
            ["Annotation/refannot_class/all_annot.xls", "", "注释结果详情表", 0,  "211338"],
            ["Annotation/refannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0,  "211339"],
            ["Annotation/refannot_class/all_stat.xls", "", "注释结果统计表", 0,  "211340"],
            ["Annotation/refannot_class/cog", "", "cog注释结果目录", 0,  "211341"],
            ["Annotation/refannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0,  "211342"],
            ["Annotation/refannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0,  "211343"],
            ["Annotation/refannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0,  "211344"],
            ["Annotation/refannot_class/go/", "", "go注释结果目录", 0,  "211345"],
            ["Annotation/refannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0,  "211346"],
            ["Annotation/refannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0,  "211347"],
            ["Annotation/refannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0,  "211348"],
            ["Annotation/refannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0,  "211349"],
            ["Annotation/refannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0,  "211350"],
            ["Annotation/refannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0,  "211351"],
            ["Annotation/refannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0,  "211352"],
            ["Annotation/refannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0,  "211353"],
            ["Annotation/refannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0,  "211354"],
            ["Annotation/refannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0,  "211355"],
            ["Annotation/refannot_class/kegg", "", "kegg注释结果目录", 0,  "211356"],
            ["Annotation/refannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0,  "211357"],
            ["Annotation/refannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0,  "211358"],
            ["Annotation/refannot_class/kegg/kegg_pathway_gene_dir:", "", "基因kegg注释通路图", 0,  "211359"],
            ["Annotation/refannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0,  "211360"],
            ["Annotation/refannot_class/kegg/kegg_pathway_tran_dir:", "", "转录本kegg通路图", 0,  "211361"],
            ["Annotation/refannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0,  "211362"],
            ["Annotation/refannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0,  "211363"],
            ["Annotation/refannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0,  "211364"],
            ["Annotation/refannot_class/nr", "", "nr注释结果目录", 0,  "211365"],
            ["Annotation/refannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0,  "211366"],
            ["Annotation/refannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0,  "211367"],
            ["Annotation/refannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0,  "211368"],
            ["Annotation/refannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0,  "211369"],
            ["Annotation/refannot_class/pfam", "", "pfam注释结果目录", 0,  "211370"],
            ["Annotation/refannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0,  "211371"],
            ["Annotation/refannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0,  "211372"],
            ["Annotation/refannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0,  "211373"],
            ["Annotation/refannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0,  "211374"],
            ["Annotation/refannot_class/swissprot", "", "swissprot注释结果目录", 0,  "211375"],
            ["Annotation/refannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0,  "211376"],
            ["Annotation/refannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0,  "211377"],
            ["Annotation/refannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0,  "211378"],
            ["Annotation/refannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0,  "211379"],
            ["Annotation/newannot_mapdb", "", "新转录本与数据库比对结果目录", 0,  "211380"],
            ["Annotation/newannot_mapdb/eggnog", "", "新转录本eggnog数据库比对结果目录", 0,  "211381"],
            ["Annotation/newannot_mapdb/go", "", "新转录本GO数据库比对结果目录", 0,  "211382"],
            ["Annotation/newannot_mapdb/kegg", "", "新转录本kegg数据库比对结果目录", 0,  "211383"],
            ["Annotation/newannot_mapdb/nr", "", "新转录本nr数据库比对结果目录", 0,  "211384"],
            ["Annotation/newannot_mapdb/swissprot", "", "新转录本swissprot数据库比对结果目录", 0,  "211385"],
            ["Express", "", "表达量分析结果目录", 0,  "211386"],
            ["Express/ExpAnnalysis", "", "表达定量结果目录", 0,  "211387"],
            ["Express/ExpAnnalysis/gene.count.matrix.xls", "", "gene count表达定量结果表", 0,  "211388"],
            ["Express/ExpAnnalysis/gene.tpm.matrix.xls", "", "gene tpm表达定量结果表", 0,  "211389"],
            ["Express/ExpAnnalysis/gene.fpkm.matrix.xls", "", "gene fpkm表达定量结果表", 0,  "211390"],
            ["Express/ExpAnnalysis/transcript.count.matrix.xls", "", "transcript count表达定量结果表", 0,  "211391"],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.xls", "", "transcript tpm表达定量结果表", 0,  "211392"],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.xls", "", "transcript fpkm表达定量结果表", 0,  "211393"],
            ["Express/ExpAnnalysis/gene.tpm.matrix.annot.xls", "", "gene tpm表达定量注释结果表", 0,  "211394"],
            ["Express/ExpAnnalysis/gene.fpkm.matrix.annot.xls", "", "gene fpkm表达定量注释结果表", 0,  "211395"],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls", "", "transcript tpm表达定量注释结果表", 0,  "211396"],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量注释结果表", 0,  "211397"],
            ["Express/ExpCorr", "", "样本间相关性分析", 0,  "211398"],
            ["Express/ExpCorr/sample_correlation.xls ", "", "样本间相关性分析矩阵表", 0,  "211399"],
            ["Express/ExpPCA", "", "样本间PCA分析", 0,  "211400"],
            ["Express/ExpPCA/PCA.xls", "", "样本间PCA分析结果表", 0,  "211401"],
            ["Express/ExpPCA/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表", 0,  "211402"],
            ["DiffExpress", "", "表达量差异分析", 0,  "211403"],
            ["SNP", "", "SNP/InDel分析结果文件", 0,  "211404"],
            ["SNP/snp_annotation_statistics.xls", "", "snp分析结果注释详情表格", 0,  "211405"],
            ["SNP/snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格", 0,  "211406"],
            ["SNP/snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0,  "211407"],
            ["SNP/snp_depth_statistics.xls", "", "SNP深度统计结果表格", 0,  "211408"],
            ["SNP/snp_position_distribution.xls", "", "SNP不同区域布析结果表格", 0,  "211409"],
            ["SNP/indel_position_distribution.xls", "", "InDel不同区域布析结果表格", 0,  "211410"],
            ["AS", "", "可变剪切分析结果文件", 0,  "211411"],
            ])

    def run_api(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.ref_rna_v2')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.stop_timeout_check()
        self.export_genome_info()
        self.export_qc()
        if self.option("map_assess_method") != None:
            self.export_map_assess()
        if self.option("is_assemble") == True:
            self.export_assembly()
        if self.option("sample_num") == "multiple":
            if self.option("is_snp").lower() == "true":
                self.export_snp()
            if self.option("is_as").lower() == "true":
                self.export_as()
        self.export_annotation()
        self.export_expression()
        self.export_gene_detail()
        self.logger.info("导表完成")

    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("ref_rna_v2.genome_info")
        file_path = self.genome_stat
        species_name = self.species_name
        species = self.species
        ref_anno_version = self.ref_anno_version
        hyperlink = self.hyperlink
        self.api_geno.add_genome_info(file_path=file_path,species_name=species_name, species=species, ref_anno_version=ref_anno_version,hyperlink=hyperlink)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("ref_rna_v2.ref_rna_qc")
        qc_stat = self.qc_stat_before.output_dir
        fq_type = self.option("fq_type").lower()
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before")
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat"
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat"  # 将qc前导表加于该处
        self.api_qc.add_gragh_info(quality_stat_before, "before")
        self.api_qc.add_bam_path(self.workflow_output)
        qc_stat = self.qc_stat_after.output_dir
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="after")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        if self.option("group_table").is_set:
            self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(self.option("group_table").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
        if self.option("control_file").is_set:
            self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control_file").prop["path"], self.group_id)
            self.compare_detail = compare_detail
        self.api_qc.add_bam_path(self.workflow_output)

    @time_count
    def export_assembly(self):
        gevent.sleep()
        self.api_assembly = self.api.api("ref_rna_v2.ref_assembly")
        if self.option("assemble_method").lower() == "cufflinks":
            all_gtf_path = self.assembly.output_dir + "/Cufflinks"
            merged_path = self.assembly.output_dir + "/Cuffmerge"
        else:
            all_gtf_path = self.assembly.output_dir + "/Stringtie"
            merged_path = self.assembly.output_dir + "/StringtieMerge"
        params = self.option("assemble_method").lower()
        self.api_assembly.add_assembly_result(all_gtf_path=all_gtf_path, params=params, merged_path=merged_path, statistics_path=self.assembly.output_dir + "/Statistics")

    @time_count
    def export_map_assess(self):
        gevent.sleep()
        self.api_map = self.api.api("ref_rna_v2.ref_rna_qc")
        stat_dir = self.mapping.output_dir + "/stat"
        if self.option("align_method").lower() == "tophat":
            self.api_map.add_tophat_mapping_stat(stat_dir)
        else:
            self.api_map.add_hisat_mapping_stat(stat_dir)
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
        self.api_annotation = self.api.api("ref_rna_v2.ref_annotation")
        annot_dir = self.merge_annot.output_dir
        if self.option("is_assemble") == False:
            self.api_annotation.has_new = False
            trans2gene  = None
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        else:
            trans2gene = annot_dir + "/newannot_class/all_tran2gene.txt"
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        self.api_annotation.species_name = self.option("ref_genome")
        params_dict = {
            "nr_evalue": str(self.option("nr_blast_evalue")),
            "nr_similarity": self.option("nr_blast_similarity"),
            "nr_identity": self.option("nr_blast_identity"),
            "swissprot_evalue":str(self.option("swissprot_blast_evalue")),
            "swissprot_similarity": self.option("swissprot_blast_similarity"),
            "swissprot_identity": self.option("swissprot_blast_identity"),
            "cog_evalue": str(self.option("cog_blast_evalue")),
            "cog_similarity": self.option("cog_blast_similarity"),
            "cog_identity": self.option("cog_blast_identity"),
            "kegg_evalue": str(self.option("kegg_blast_evalue")),
            "kegg_similarity": self.option("kegg_blast_similarity"),
            "kegg_identity": self.option("kegg_blast_identity"),
            "pfam_evalue": str(self.option("pfam_blast_evalue")),
        }
        self.api_annotation.run(
            annot_dir,
            trans2gene,
            trans2gene_ref,
            params_dict=params_dict,
            taxonomy=self.option("kegg_database"),
            exp_level=self.option("level").lower(),
            version = "v3",
            gene_exp = self.work_dir + "/gene.tpm.matrix",
            trans_exp = self.work_dir + "/transcript.tpm.matrix",
        )

    @time_count
    def export_as(self):
        gevent.sleep()
        self.api_as = self.api.api("ref_rna_v2.splicing_rmats")
        for file in os.listdir(self.altersplicing.output_dir):
            tmp_group_list = file.split("_vs_")
            group_a = tmp_group_list[0]
            group_b = tmp_group_list[1]
            compare_plan = group_a + "|" + group_b
            group_a_member = self.option('group_table').prop['group_dict'][group_a]
            group_b_member = self.option('group_table').prop['group_dict'][group_b]
            group_dict = dict()
            group_dict[group_a] = group_a_member
            group_dict[group_b] = group_b_member
            params = {
            "submit_location": "splicingrmats",
            "task_type": 2,
            "group_id": str(self.group_id),
            "compare_plan": compare_plan,
            "control_id": str(self.control_id),
            "group_dict": group_dict,
            "task_id": self.task_id
        }
            self.logger.info(params)
            outpath = os.path.join(self.altersplicing.output_dir, file)
            self.api_as.add_sg_splicing_rmats(params=params, major=True, group_dict=group_dict, compare_plan=compare_plan,name=None, outpath=outpath)

    @time_count
    def export_snp(self):
        gevent.sleep()
        self.logger.info("开始进行Snpfinal的导表")
        task_id = self.task_id
        project_sn = self.project_sn
        if self.option("snp_method").lower() == "gatk":
            params = dict(
                task_id=task_id,
                submit_location="snp",
                task_type=2,
                method_type="gatk"
            )
            self.api_snp = self.api.api("ref_rna_v2.ref_snp")
            snp_anno = self.snp.output_dir
            if os.path.exists(self.work_dir + "/snp_tmp"):
                shutil.rmtree(self.work_dir + "/snp_tmp")
            os.mkdir(self.work_dir + "/snp_tmp")

            self.api_snp.add_snp_main(snp_anno=snp_anno, params=params,
                                      task_id=task_id, method_type="gatk", project_sn=project_sn, new_output=self.work_dir + "/snp_tmp")
        if self.option("snp_method").lower() == "samtools":
            params = dict(
                task_id=task_id,
                submit_location="snp",
                task_type=2,
                method_type="samtools"
            )
            self.api_snp = self.api.api("ref_rna_v2.ref_snp")
            snp_anno = self.snp.output_dir
            if os.path.exists(self.work_dir + "/snp_tmp"):
                shutil.rmtree(self.work_dir + "/snp_tmp")
            os.mkdir(self.work_dir + "/snp_tmp")
            # os.link(os.path.join(snp_anno, 'snp_anno.xls'), os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))
            self.api_snp.add_snp_main(snp_anno=snp_anno, params=params,
                                      task_id=task_id, method_type="samtools", project_sn=project_sn, new_output=self.work_dir + "/snp_tmp")
        if os.path.exists(os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls')):
            os.remove(os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))
        os.link(os.path.join(snp_anno, 'data_anno_pre.xls'), os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))

    @time_count
    def export_expression(self):
        gevent.sleep()
        all_exp = self.api.api("ref_rna_v2.all_exp")
        quant_method = self.option('express_method')
        task_id = self.task_id
        project_sn = self.project_sn
        if self.option("group_table").is_set:
            group_dict = self.option('group_table').prop['group_dict']
            group_id = self.group_id
        else:
            group_dict = None
            group_id = None
        if self.option("control_file").is_set:
            control_id = self.control_id
        else:
            control_id = None

        # add exp matrix
        ## 还需确认路径信息，因为后续交互需要用到express.output_dir
        exp_output = self.express.output_dir
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
        if self.option("strand_specific") == True:
            if self.option('strand_dir') == 'forward':
                if self.option("fq_type") == "PE":
                    self.libtype = "rf"
                else:
                    self.libtype = "r"
            else:
                if self.option("fq_type") == "PE":
                    self.libtype = "fr"
                else:
                    self.libtype = "f"
        else:
            self.libtype = None
        if self.option("exp_way").lower() == "tpm":
            if self.option("level").lower() == "transcript":
                exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
                ref_exp_matrix = self.work_dir + "/transcript.tpm.matrix"
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T', lib_type=self.libtype,
                                               group_dict=group_dict,  group_id=group_id, add_distribution=False,
                                               exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
                if self.option("sample_num") == "multiple":
                    params_distribution = dict(
                        task_id=task_id,
                        exp_id=str(trans_exp_id),
                        group_dict=group_dict,
                        group_id=str(group_id),
                        submit_location="expgraph",
                        task_type=2,
                        exp_level='T',
                        type='ref',
                    )
                    all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='T',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)

            exp_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
            ref_exp_matrix = self.work_dir + "/gene.tpm.matrix"
            gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G', lib_type=self.libtype,
                                          group_dict=group_dict,  group_id=group_id, add_distribution=False,
                                          exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
            if self.option("sample_num") == "multiple":
                params_distribution = dict(
                    task_id=task_id,
                    exp_id=str(gene_exp_id),
                    group_dict=group_dict,
                    group_id=str(group_id),
                    submit_location="expgraph",
                    task_type=2,
                    exp_level='G',
                    type='ref',
                )
                all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='G',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)
        else:
            if self.option("level").lower() == "transcript":
                exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
                ref_exp_matrix = self.work_dir + "/transcript.fpkm.matrix"
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T', lib_type=self.libtype,
                                               group_dict=group_dict,  group_id=group_id, add_distribution=False,
                                               exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
                if self.option("sample_num") == "multiple":
                    params_distribution = dict(
                        task_id=task_id,
                        exp_id=str(trans_exp_id),
                        group_dict=group_dict,
                        group_id=str(group_id),
                        submit_location="expgraph",
                        task_type=2,
                        exp_level='T',
                        type='ref',
                    )
                    all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='T',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)

            exp_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            ref_exp_matrix = self.work_dir + "/gene.fpkm.matrix"
            gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G', lib_type=self.libtype,
                                          group_dict=group_dict,  group_id=group_id, add_distribution=False,
                                          exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
            if self.option("sample_num") == "multiple":
                params_distribution = dict(
                    task_id=task_id,
                    exp_id=str(gene_exp_id),
                    group_dict=group_dict,
                    group_id=str(group_id),
                    submit_location="expgraph",
                    task_type=2,
                    exp_level='G',
                    type='ref',
                )
                all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='G',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)

        if self.option("sample_num") == "multiple":
            # add gene corr
            corr_output = self.exp_corr.work_dir
            params = dict(
                task_id=task_id,
                submit_location='expcorr',
                task_type=2,
                exp_id=str(gene_exp_id),
                group_id=str(group_id),
                exp_level="G",
                group_dict=group_dict,
                scm="complete",
                scd="euclidean",
                # quant_method=quant_method,
                corr_method="pearson",
                type="ref"
            )
            all_exp.add_exp_corr2(corr_output, exp_level='G', quant_method=quant_method, params=params,
                                  project_sn=project_sn, task_id=task_id)

            # add gene pca
            if self.option("group_table").prop["sample_number"] > 2:
                pca_output = self.exp_pca.output_dir
                params = dict(
                    task_id=task_id,
                    submit_location="exppca",
                    task_type=2,
                    exp_id=str(gene_exp_id),
                    group_id=str(group_id),
                    exp_level="G",
                    group_dict=group_dict,
                    type="ref"
                    # quant_method=quant_method,
                )
                all_exp.add_exp_pca2(pca_output, quant_method=quant_method, exp_id=gene_exp_id, exp_level="G",
                                     params=params, project_sn=project_sn, task_id=task_id)

            # add gene diff
            diff_output = self.diffexpress.output_dir
            exp_id, exp_level = gene_exp_id, 'G'
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
                type='ref'
            )
            all_exp.add_diffexp(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
                                exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
                                project_sn=project_sn, task_id=task_id, params=params,
                                pvalue_padjust=stat_type)

    @time_count
    def export_gene_detail(self):
        """
        导入基因详情表
        :return:
        """
        gevent.sleep()
        self.api_gene_detail = self.api.api('ref_rna_v2.add_gene_detail')
        db_path = self.output_dir + "/Sequence_database/refrna_seqs.db"
        if os.path.exists(self.output_dir + "/Sequence_database/"):
            shutil.rmtree(self.output_dir + "/Sequence_database/")
        os.mkdir(self.output_dir + "/Sequence_database/")
        gene_bed = self.gene_fa.option("gene_bed").prop["path"]
        transcript_bed = self.gene_fa.option("transcript_bed").prop["path"]
        if self.option("is_assemble") == True:
            transcript_fasta = self.assembly.option("all_transcripts_fa").prop["path"]
            t2g = self.assembly.option("trans2gene").prop["path"]
        else:
            transcript_fasta = self.transcripts.option("trans_fa").prop["path"]
            t2g = self.transcripts.option("trans2gene").prop["path"]
        gene_fasta = self.gene_fa.option("gene_fa").prop["path"]
        species_urls = self.hyperlink
        biomart_file = self.des
        biomart_type = self.des_type
        known_cds = self.known_cds
        known_pep = self.known_pep
        if self.option("is_assemble") == True:
            new_cds = self.annot_orfpfam.output_dir + "/new_transcripts.fa.transdecoder.cds"
            new_pep = self.annot_orfpfam.output_dir + "/new_transcripts.fa.transdecoder.pep"
        else:
            new_cds = None
            new_pep = None

        self.api_gene_detail.add_gene_detail(db_path, gene_bed, transcript_bed, species_urls, biomart_file, biomart_type,
                        known_cds, known_pep, new_cds, new_pep, transcript_fasta, gene_fasta, t2g)

    # 添加注释信息
    def merge_annotation_exp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        exp_output = self.express.output_dir
        if self.option("is_assemble") == True:
            annot = os.path.join(self.output_dir, 'merge_annot/allannot_class/all_annot.xls')
        else:
            annot = os.path.join(self.output_dir, 'merge_annot/refannot_class/all_annot.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript_id', 'is_gene'])
        trans_annot_pd = all_annot.reset_index().drop(columns=['gene_id', 'is_gene']).set_index('transcript_id')

        # for gene
        # gene tpm
        gene_tpm_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
        gene_pd = pd.read_table(gene_tpm_matrix, header=0, index_col=0)
        gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1)
        gene_out = os.path.join(exp_output, 'gene.tpm.matrix.annot.xls')
        gene_result.to_csv(gene_out, header=True, index=True, sep='\t')
        # gene fpkm
        if self.option("express_method") == "RSEM":
            gene_fpkm_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            gene_pd = pd.read_table(gene_fpkm_matrix, header=0, index_col=0)
            gene_result = pd.concat([gene_pd, gene_annot_pd], join='inner', axis=1)
            gene_out = os.path.join(exp_output, 'gene.fpkm.matrix.annot.xls')
            gene_result.to_csv(gene_out, header=True, index=True, sep='\t')

        # for transcript
        if self.option("level").lower() == "transcript":
            # trans tpm
            trans_tpm_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            trans_pd = pd.read_table(trans_tpm_matrix, header=0, index_col=0)
            trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
            trans_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            trans_result.to_csv(trans_out, header=True, index=True, sep='\t')
            # trans fpkm
            if self.option("express_method") == "RSEM":
                trans_fpkm_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
                trans_pd = pd.read_table(trans_fpkm_matrix, header=0, index_col=0)
                trans_result = pd.concat([trans_pd, trans_annot_pd], join='inner', axis=1)
                trans_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
                trans_result.to_csv(trans_out, header=True, index=True, sep='\t')

    def merge_annotation_diffexp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        if self.option("is_assemble") == True:
            annot = os.path.join(self.output_dir, 'merge_annot/allannot_class/all_annot.xls')
        else:
            annot = os.path.join(self.output_dir, 'merge_annot/refannot_class/all_annot.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript_id', 'is_gene'])
        diff_output = self.diffexpress.output_dir

        duplicate_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls') ## 为了防止流程重跑的时候反复增加注释结果表
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + '/' + '*_vs_*.*.xls')
        for each in target_files:
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_pd, gene_annot_pd], join='inner', axis=1)
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
        df_join_pre.to_csv(self.work_dir + "/upload_results/SNP/snp_annotation_statistics.xls", sep="\t", index=False)
