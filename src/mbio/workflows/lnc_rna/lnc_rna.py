# -*- coding:utf-8 -*-
# __author__ = 'liubinxu'
"""lncrna一键化工作流"""

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
from mbio.packages.lnc_rna.copy_file import CopyFile
from mbio.packages.lnc_rna.lnc_file_des import LncFileDes
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo

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


class LncRnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(LncRnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 单样本 OR 多样本
            {"name": "sample_num", "type": "string", "default": "multiple"},  # 测序类型，single OR multiple
            ## 基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE
            {"name": "datatype", "type": "string", "default": "rawdata"},
            {"name": "qc_soft", "type" : "string", "default": "seqprep"},
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "strand_specific", "type": "bool", "default": True}, # 当为PE测序时，是否为链特异性
            {"name": "strand_dir", "type": "string", "default": "RF"}, # 当链特异性时为True时，RF R or FR F
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "taxonmy", "type":"string", "default": "Animal"}, # 物种类别
            {"name": "ref_genome", "type": "string", "default": "Custom"},  # 参考基因组，具体物种名称
            {"name": "genome_version", "type": "string", "default": "Custom"},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": "Custom"},  # 参考基因组注释版本
            {"name": "genome_id", "type": "string", "default": None},  # 参考基因组ID 参考基因组数据库中唯一ID
            {"name": "knownlnc_fasta", "type": "infile", "format": "lnc_rna.fasta"},  # 已知lncRNA 序列文件
            {"name": "knownlnc_evalue", "type": "string", "default": "1e-3"},

            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {"name": "knownlnc_qcov", "type": "float", "default": 80},
            {"name": "knownlnc_scov", "type": "float", "default": 80},
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "qc_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            # 上机名称
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'},
            {"name": "control_file", "type": "infile", "format": "sample.control_table"}, # 对照表
            {"name": "mirbase_category", "type":"string", "default": "Animal"}, # mirbase物种类别，Animal OR Plant
            {"name": "mirbase_specie", "type": "string", "default": ""},  # mirbase参考物种名，多选时分号分隔
            {"name": "mapping_stop", "type": "string", "default": 'False'}, # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_sample_percent", "type": "float", "default": 50}, # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_ratio", "type": "float", "default": 60}, # 设定比对率小于多少时判定样本不达标
            {"name": "rrna_stop", "type": "string", "default": 'False'}, # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_sample_percent", "type": "float", "default": 50}, # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_ratio", "type": "float", "default": 15}, # 设定rRNA比率大于多少时判定样本不达标

            ## 高级参数设置
            # 分析水平
            {"name": "level", "type": "string", "default": "gene"}, # 分析水平，gene or genetrans
            # 序列比对分析
            {"name": "align_method", "type": "string", "default": "Hisat"},  # 比对方法，Tophat or Hisat
            {"name": "map_assess_method", "type": "string", "default": "rrna,saturation,distribution,coverage,chr_stat"},#质量评估分析项
            # 转录本组装
            {"name": "is_assemble", "type": "string", "default": True}, # 是否进行有参拼接
            {"name": "is_filter", "type": "string", "default": "True"}, # 是否进行过滤
            {"name": "assemble_method", "type": "string", "default": "stringtie"},# 拼接方法，Cufflinks or Stringtie
            # lncRNA鉴定与预测
            {"name": "lnc_ref_db", "type": "string", "default": ""},
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
            {"name": "cpc", "type": "string", "default": "True"},
            {"name": "cnci", "type": "string", "default": "True"},
            {"name": "cpat", "type": "string", "default": "False"},
            {"name": "pfamscan", "type": "string", "default": "True"},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'lnc_identify_num', 'type': 'int', 'default': 2},
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
            {"name": "filter_tpm", "type": "float", "default": -1},  # 表达量过滤标准
            # 表达差异分析
            {"name": "filter_method", "type": "string", "default": "None"},  # 表达量过滤标准 差异分析时使用
            {"name": "tpm_filter_threshold", "type": "string", "default": "NA"},  # 表达量过滤标准 差异分析时使用
            {"name": "diff_method", "type": "string", "default": "DESeq2"},# 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "diff_fdr_ci", "type": "float", "default": 0.05},  # 显著性水平
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  #选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  #Bonferroni,Holm,BH,BY
            # lncRNA靶基因分析
            {"name": "is_target_cistrans", "type": "string", "default": "True"}, #是否进行cis把基因预测
            {"name": "is_target_cis", "type": "string", "default": "True"}, #是否进行cis把基因预测
            {"name": "is_target_trans", "type": "string", "default": "True"}, #是否进行trans靶基因预测
            {"name": "target_trans_method", "type": "string", "default": "rnaplex"}, #trans靶基因预测方法
            {"name": "is_mirna_pre", "type": "string", "default": "True"}, #是否进行mirna前体预测
            {'name': 'target_pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'target_pvalue_type', 'type': 'string', 'default': "padjust"},
            {'name': 'target_cor_cutoff', 'type': 'float', 'default': 0.9},
            {'name': 'target_corr_way', 'type': 'string', 'default': "spearmanr"},
            {'name': 'target_padjust_way', 'type': 'string', 'default': "BH"},
            {'name': 'up_dis', 'type': 'int', 'default': 10},
            {'name': 'down_dis', 'type': 'int', 'default': 10},

            {"name": "is_loc_cons", "type": "string", "default": "True"}, #是否进行位置保守性分析
            {"name": "is_seq_cons", "type": "string", "default": "True"}, #是否进行序列保守性分析
            {"name": "is_rfam", "type": "string", "default": "True"}, # 是否比对rfam数据库

            # SNP/INDEL分析
            {"name": "is_snp", "type": "string", "default": "True"}, #是否进行snp分析,True, False, Skip
            {"name": "snp_method", "type": "string", "default": "samtools"}, #samtools or GATK
            {"name": "qual", "type": "float", "default": 20},  #过滤低质量的SNP位点
            {"name": "dp", "type": "int", "default": 1},  #过滤低质量的SNP位点
            # 可变剪切分析
            {"name": "is_as", "type": "string", "default": "True"}, #是否进行as分析,True, False, Skip
            {"name": "as_method", "type": "string", "default": "rmats"}, #rmats or asprofile
            {"name": "version", "type": "string", "default": "v1.1"},
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_2019"},
        ]
        #获取输出目录

        # taxon dict 页面上传物种分类与module对应关系
        self.taxon_dict = {
            "All": ["nr", "All"],
            "Animal": ["metazoa", "Animals"],
            "Plant": ["viridiplantae", "Plants"],
            "Fungi": ["fungi", "Fungi"],
            "Protist": ["protist", "Protists"],
        }

        self.workflow_output_tmp = self._sheet.output
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

        self.logger.info("fc is {}".format(self.option("fc")))
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))

        #获取数据库中该物种相关文件路径
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        col = db["sg_genome_db"]
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"

        db_small = Config().get_mongo_client(mtype="small_rna", dydb_forbid=True)[Config().get_mongo_dbname("small_rna", dydb_forbid=True)]
        small_col = db_small["mirbase"]
        try:
            if self.option("genome_id"):
                self.genome_info = col.find_one({"genome_id" : self.option("genome_id")})
            else:
                self.genome_info = col.find_one({"organism_name" : self.option("ref_genome"), "assembly" : self.option("genome_version"), "annot_version" : self.option("genome_annot_version")})
        except:
            self.set_error("数据库中不存在该物种注释信息，程序退出")

        if not self.genome_info:
            self.set_error("数据库中不存在该物种注释信息，程序退出")
        try:
            # print self.genome_info
            self.ref_annot_dir = os.path.join(db_path, self.genome_info["anno_path_v2"])
            self.ref_genome = os.path.join(db_path, self.genome_info["dna_fa"])
            self.ref_gtf = os.path.join(db_path, self.genome_info["gtf"])
            self.genome_id = self.genome_info["genome_id"]
            self.des = os.path.join(db_path, self.genome_info["bio_mart_annot"])
            self.des_type = self.genome_info["biomart_gene_annotype"]
            self.known_cds = os.path.join(db_path, self.genome_info["cds"])
            self.known_pep = os.path.join(db_path, self.genome_info["pep"])
            '''
            self.lnc_dir = os.path.join(db_path, self.genome_info["lnc_dir"])
            self.lnc_gtf = os.path.join(db_path, self.genome_info["lnc_dir"], "lncrna.gtf")
            self.lnc_ids = os.path.join(db_path, self.genome_info["lnc_dir"], "ids_matrix.xls")
            '''
            if "lnc_dir" in self.genome_info and \
               os.path.exists(os.path.join(db_path, self.genome_info["lnc_dir"], "mrna.gtf")) and \
               os.path.exists(os.path.join(db_path, self.genome_info["lnc_dir"], "mrna.fa")):
                self.mrna_gtf = os.path.join(db_path, self.genome_info["lnc_dir"], "mrna.gtf")
                self.mrna_fa = os.path.join(db_path, self.genome_info["lnc_dir"], "mrna.fa")
            else:
                self.mrna_gtf = self.ref_gtf
                self.mrna_fa =  os.path.join(db_path, self.genome_info["transcript"])
            if "lnc_dir" in self.genome_info and \
               os.path.exists(os.path.join(db_path, self.genome_info["lnc_dir"], "lncrna.gtf")) and \
               os.path.exists(os.path.join(db_path, self.genome_info["lnc_dir"], "lncrna.fa")):
                self.lnc_dir = os.path.join(db_path, self.genome_info["lnc_dir"])
                self.lnc_fa = os.path.join(db_path, self.genome_info["lnc_dir"], "lncrna.fa")
                self.lnc_gtf = os.path.join(db_path, self.genome_info["lnc_dir"], "lncrna.gtf")
                self.lnc_ids = os.path.join(db_path, self.genome_info["lnc_dir"], "ids_matrix.xls")
                self.has_known_lnc = True
            else:
                self.has_known_lnc = False
                os.system("touch {}".format(self.work_dir + '/lncrna.fa'))
                os.system("touch {}".format(self.work_dir + '/lncrna.gtf'))
                os.system("touch {}".format(self.work_dir + '/ids_matrix.xls'))
                self.lnc_fa = self.work_dir + '/lncrna.fa'
                self.lnc_gtf = self.work_dir + '/lncrna.gtf'
                self.lnc_ids = os.path.join(db_path, self.work_dir + '/ids_matrix.xls')

            self.entrez = os.path.join(db_path, self.genome_info["ensemble2entrez"])
            self.genome_stat = os.path.join(db_path, self.genome_info["gene_stat"])
            self.g2t2p = os.path.join(db_path, self.genome_info['g2t2p'])
            self.species_name = self.genome_info["name"]
            self.species = self.genome_info["taxon_id"]
            self.ref_anno_version = self.genome_info["assembly"]
            self.hyperlink = self.genome_info["ensemble_web"]
            self.known_ko = self.genome_info['kegg']
        except Exception as e:
            self.set_error("数据库中该物种注释信息不全 {}".format(e))

        self.mirpre =  True

        if self.option("is_mirna_pre").lower() == "true":
            if self.option("mirbase_specie") == "" or self.option("mirbase_specie").lower() == "none":
                self.mirpre = False
            else:
                self.mirna_info = small_col.find_one({"Name" : self.option("mirbase_specie")})

                if not self.mirna_info:
                    self.mirna_info = small_col.find_one({"organism" : self.option("mirbase_specie")})
                if not self.mirna_info:
                    self.set_error("small RNA物种数据库选择错误")
        else:
            self.mirpre  = False
            pass



        for ref_file in ["anno_path_v2", "dna_fa", "gtf", "bio_mart_annot", "cds", "pep", "ensemble2entrez", "gene_stat", 'g2t2p']:
            if os.path.exists(os.path.join(db_path, self.genome_info[ref_file])):
                pass
            else:
                self.set_error("数据库中该物种注释信息不全,文件{}不存在".format(os.path.join(db_path, self.genome_info[ref_file])))


        #添加tool/module
        self.filecheck = self.add_tool("lnc_rna.file_check")
        if self.option('qc_soft') == "fastp":
            self.qc = self.add_module("datasplit.fastp_rna")
        else:
            self.qc = self.add_module("lnc_rna.lncrna_qc")

        self.qc_stat_before = self.add_module("lnc_rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("lnc_rna.hiseq_reads_stat")
        self.mapping = self.add_module("lnc_rna.rnaseq_mapping")
        self.assembly = self.add_module("lnc_rna.assemble")

        self.logger.info("before import  lnc_indent")
        self.lncrna_identify = self.add_module("lnc_rna.known_lnc_identify")
        self.logger.info("before import  lnc_stat")
        self.lncrna_stat = self.add_tool("lnc_rna.lncrna_identification.lncrna_stat")
        self.logger.info("before import  lnc_annot")
        self.lncrna_new_annot = self.add_tool("lnc_rna.lncrna_annot")
        self.logger.info("before import  predict")
        self.lncrna_new = self.add_module("lnc_rna.new_lncrna_predict")

        self.annot_filter_ref = self.add_module("lnc_rna.annot_filter")
        self.annot_class_ref = self.add_module("lnc_rna.annot_class")
        self.annot_mapdb = self.add_module("lnc_rna.annot_mapdb")
        self.annot_class_new = self.add_module("lnc_rna.annot_class")
        self.annot_orfpfam = self.add_module("lnc_rna.annot_orfpfam")
        self.annot_filter_new = self.add_module("lnc_rna.annot_filter")
        self.merge_annot =  self.add_tool("lnc_rna.annotation.annot_merge")
        self.merge_known_new = self.add_tool("lnc_rna.merge_known_new")
        self.filter_by_exp = self.add_tool("lnc_rna.filter_by_express")

        self.classify_quant = self.add_tool("lnc_rna.classify_quant")
        if self.option("map_assess_method") != None:
            self.map_qc = self.add_module("lnc_rna.map_assessment")

        self.target_cistrans = self.add_module("lnc_rna.target_cistrans")
        '''
        if self.option("is_target_cis").lower() == "true":
            self.target_cis_known = self.add_tool("lnc_rna.lnc_target_cis")
            self.target_cis_novel = self.add_tool("lnc_rna.lnc_target_cis")
        if self.option("is_target_trans").lower() == "true":
            self.target_trans_known = self.add_module("lnc_rna.lnc_target_trans")
            self.target_trans_novel = self.add_module("lnc_rna.lnc_target_trans")
        '''
        self.lnc_family = self.add_module("lnc_rna.lncrna_family")
        if self.mirpre:
            self.mir_pre = self.add_module("lnc_rna.mirna_precursor")


        if self.option("is_loc_cons").lower() == "true":
            self.loc_cons = self.add_tool("lnc_rna.structure.lncrna_loc_cons")

        if self.option("is_seq_cons").lower() == "true":
            self.seq_cons = self.add_tool("lnc_rna.structure.lncrna_seq_cons")

        if self.option("is_snp").lower() == "true":
            self.snp = self.add_module("lnc_rna.snp_indel")
        if self.option("is_as").lower() == "true":
            self.altersplicing = self.add_module("lnc_rna.rmats")
        self.express = self.add_module("lnc_rna.quant")
        self.exp_pca = self.add_tool("lnc_rna.exp_pca")
        self.exp_pca_t = self.add_tool("lnc_rna.exp_pca")
        self.exp_corr = self.add_tool("lnc_rna.exp_corr")
        self.exp_corr_t = self.add_tool("lnc_rna.exp_corr")
        # self.diffexpress = self.add_tool("lnc_rna.diffexp")
        # self.diffexpress_t = self.add_tool("lnc_rna.diffexp")
        self.diffexpress = self.add_module("lnc_rna.diffexp_batch_new")
        self.diffexpress_t = self.add_module("lnc_rna.diffexp_batch_new")
        self.gene_fa = self.add_tool("lnc_rna.gene_fa")
        self.transcripts = self.add_tool("ref_rna_v2.transcript_abstract")
        # self.ellipse = self.add_tool("graph.ellipse")
        self.ellipse = self.add_tool("ref_rna_v3.ellipse")
        self.ellipse_t = self.add_tool("graph.ellipse")

        #判断流程结束tool/module list
        if self.option("sample_num") == "multiple":
            self.final_tools = [self.qc_stat_after, self.annot_filter_ref, self.annot_class_ref, self.express, self.diffexpress, self.diffexpress_t, self.merge_annot, self.gene_fa]

            group_spname =  self.option("group_table").get_group_spname()
            group_snum = [len(group_spname[g]) for g in group_spname]
            self.samples_num = sum(group_snum)
            self.min_group_num = min(group_snum)
            if self.option("is_snp").lower() == "true":
                self.final_tools.append(self.snp)
            if self.option("group_table").prop["sample_number"] > 2:
                self.final_tools.append(self.exp_pca)
                if self.min_group_num >= 3:
                    self.final_tools.append(self.ellipse)
                    self.final_tools.append(self.ellipse_t)
            if self.option("is_as").lower() == "true":
                self.final_tools.append(self.altersplicing)
            self.logger.info(self.final_tools)
        else:
            self.samples_num = sum(group_snum)
            self.final_tools = [self.qc_stat_after, self.annot_filter_ref, self.annot_class_ref, self.express, self.merge_annot, self.gene_fa]
        if self.option("map_assess_method") != None:
            self.final_tools.append(self.map_qc)



        self.final_tools.append(self.assembly)
        self.final_tools.append(self.annot_mapdb)
        self.final_tools.append(self.annot_class_new)
        self.final_tools.append(self.annot_orfpfam)
        self.final_tools.append(self.annot_filter_new)
        self.final_tools.append(self.lncrna_stat)
        self.final_tools.append(self.lnc_family)

        if self.mirpre:
            self.final_tools.append(self.mir_pre)

        if self.option("sample_num") == "multiple" and self.samples_num >= 3:
            self.final_tools.append(self.target_cistrans)

        '''
        if self.option("is_target_cis").lower() == "true":
            self.final_tools.append(self.target_cis_known)
            self.final_tools.append(self.target_cis_novel)
        if self.option("is_target_trans").lower() == "true" and self.option("sample_num") == "multiple":
            self.final_tools.append(self.target_trans_known)
            self.final_tools.append(self.target_trans_novel)
        '''
        if self.option("knownlnc_fasta").is_set:
            self.final_tools.append(self.lncrna_new_annot)
        '''
        if self.option("is_loc_cons").lower() == "true":
            self.final_tools.append(self.loc_cons)

        if self.option("is_seq_cons").lower() == "true":
            self.final_tools.append(self.seq_cons)
        '''


        self.end_step = "all"


        # 添加step，显示在页面进度条
        self.step.add_steps("filecheck", "rna_qc", "qc_stat_before", "qc_stat_after", "mapping",
                            "express", "assembly", "map_qc", "merge_known_new", "filter_by_exp", "lnc_family", "classify_quant",
                            "lncrna_identify", "lncrna_new", "lncrna_stat",
                            "annot_mapdb", "annot_orfpfam", "merge_annot", "annot_class_new", "annot_class_ref", "annot_filter_new", "annot_filter_ref",
                            "transcripts_fa", "gene_fa")
        if self.option("is_snp").lower() == "true":
            self.step.add_steps("snp_rna")
        if self.option("sample_num") == "multiple":
            self.step.add_steps("diffexpress", "diffexpress_t", "exp_corr", "exp_corr_t")
            if self.samples_num >= 3:
                self.step.add_steps("target_cistrans")
            if self.option("is_as").lower() == "true":
                    self.step.add_steps("altersplicing")
            if self.option("group_table").prop["sample_number"] > 2:
                self.step.add_steps("exp_pca")
                self.step.add_steps("exp_pca_t")
                if self.min_group_num >= 3:
                    self.step.add_steps("ellipse")
                    self.step.add_steps("ellipse_t")
        '''
        if self.option("is_target_cis").lower() == "true":
            self.step.add_steps("target_cis_known", "target_cis_novel")
        if self.option("is_target_trans").lower() == "true" and self.option("sample_num") == "multiple":
            self.step.add_steps("target_trans_known", "target_trans_novel")
        '''
        if self.option("is_mirna_pre").lower()== "true":
            self.step.add_steps("mir_pre")
        if self.option("knownlnc_fasta").is_set:
            self.step.add_steps("lncrna_new_annot")
        self.logger.info("init end")
        '''
        if self.option("is_loc_cons"):
            self.step.add_steps("loc_cons")
        if self.option("is_seq_cons"):
            self.step.add_steps("seq_cons")
        '''

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        # data = os.path.join(self.work_dir, 'data.json')
        # if os.path.exists(data):
        #     with open(data, 'r') as load_f:
        #         load_dict = json.load(load_f)
        #         if 'rerun' in load_dict and load_dict['rerun']:
        #             self.logger.info("该项目重运行中，先删除mongo库中已有数据")
        #             self.delete_mongo_data()
        if self._sheet.rerun:
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'lnc_rna')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))
        delete = DeleteDemoMongo(self.task_id, 'lnc_rna')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'ref_rna_v2')
        # code = os.system(cmd)
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")


    def check_options(self):
        group_dict = self.option("group_table").prop['group_dict']
        dup = False
        for g in group_dict.keys():
            if len(group_dict[g]) > 1:
                dup = True
                pass
        if self.option("is_duplicate") != dup:
            self.logger.info(self.option("is_duplicate"))
            self.logger.info(dup)
            raise OptionError("生物学重复参数选择和分组方案不匹配，请检查输入是否有误")
        if self.option('sample_num') == 'single':
            self.option('is_duplicate', False)
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE")
        try:
            nr_evalue = float(self.option("nr_evalue"))
            cog_evalue = float(self.option("cog_evalue"))
            kegg_evalue = float(self.option("kegg_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("cog_blast_evalue", cog_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围")
        if not self.option("cog_blast_evalue") > 0 and not self.option("cog_blast_evalue") < 1:
            raise OptionError("Cog比对的E值超出范围")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围")
        if not self.option("align_method").lower() in ["tophat", "hisat"]:
            raise OptionError("比对软件应在Tophat与Hisat中选择", code = "13700316")
        for i in self.option('map_assess_method').split(','):
            if i.lower() not in ["rrna", "saturation", "distribution", "coverage", "chr_stat"]:
                raise OptionError("比对质量评估分析没有%s，请检查", variables = (i))
        if self.option("ref_genome") not in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
            # 参数隐藏
            self.option("cpat", "False")
            # if self.option("cpat").lower() == "true":
            #     raise OptionError("该物种不能使用cpat预测lincRNA，请检查")

        if self.option('sample_num') == 'multiple':
            group_size = list()
            group_dict = self.option("group_table").prop['group_dict']
            for key in group_dict:
                group_size.append(len(group_dict[key]))
            group_size.sort()
            if group_size[0] == group_size[-1] == 1:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "DEGseq")
                    self.option("diff_fdr_ci", "0.001")
                    self.logger.info("该项目没有生物学重复,不可以使用DESeq2,修改方法为DEGseq,阈值设置为0.001")
            elif group_size[0] == 1 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "edgeR")
                    self.option("diff_fdr_ci", "0.05")
                    self.logger.info("该项目部分组别有生物学重复,部分组别无生物学重复,不可以使用DESeq2,修改方法为edgeR,阈值设置为0.05")
            elif group_size[0] >= 2 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "degseq":
                    self.option("diff_method", "DESeq2")
                    self.option("diff_fdr_ci", "0.05")
                    self.logger.info("该项目有生物学重复,不可以使用DEGseq,DESeq2,阈值设置为0.05")
            else:
                pass

        if self.option('mapping_sample_percent') or self.option('rrna_sample_percent'):
            try:
                mapping_sample_percent = float(self.option("mapping_sample_percent"))
                mapping_ratio = float(self.option("mapping_ratio"))
                rrna_sample_percent = float(self.option("rrna_sample_percent"))
                rrna_ratio = float(self.option("rrna_ratio"))
            except:
                raise OptionError("判断终止那边不要传入非数字形式")
            else:
                if mapping_ratio and not mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整")
                if not mapping_ratio and mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整")
                if rrna_sample_percent and not rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整")
                if not rrna_sample_percent and rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整")
        return True
        '''
        if self.option("is_assemble") == True:
            if self.option("assemble_method").lower() not in ["cufflinks", "stringtie"]:
                raise OptionError("拼接软件应在cufflinks和stringtie中选择")
        '''
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_json(self):
        pass
        # f = open(self.json_path, "r")
        # json_dict = json.loads(f.read())
        # return json_dict

    def run(self):
        """
        ref-rna workflow run方法
        :return:
        """
        if self.option("fastq_dir").is_set:
            self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计

        if self.option("qc_dir").is_set:
            self.filecheck.on('end', self.run_mapping)  # 质控前统计
            self.filecheck.on('end', self.run_qc_stat, True)  # 质控后统计
        else:
            self.filecheck.on('end', self.run_qc)
            self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
            self.qc.on('end', self.run_mapping)
        self.mapping.on('end', self.run_assembly)

        if self.option("map_assess_method") != None:
            self.mapping.on('end', self.run_map_assess)

        self.on_rely([self.qc_stat_after, self.map_qc], self.run_check)
        self.assembly.on('end', self.run_new_lncrna)
        self.lncrna_new.on('end', self.run_merge_new)
        self.merge_known_new.on('end', self.run_express)
        self.express.on('end', self.run_lncrna_identify)
        self.express.on('end', self.run_filter_by_exp)
        self.filter_by_exp.on('end', self.run_annot_filter_ref)
        self.annot_filter_ref.on('end', self.run_annot_class_ref)

        self.filter_by_exp.on('end', self.run_annot_mapdb)
        self.filter_by_exp.on('end', self.run_annot_orfpfam)
        self.filter_by_exp.on('end', self.run_classify_quant)
        self.filter_by_exp.on('end', self.run_lnc_family)

        if self.option("is_snp").lower()== "true":
            self.filter_by_exp.on('end', self.run_snp)
        if self.option("is_as").lower()== "true" and self.option("sample_num") == "multiple":
            self.filter_by_exp.on('end', self.run_altersplicing)

        self.on_rely([self.lncrna_identify, self.filter_by_exp], self.run_lncrna_stat)

        if self.mirpre:
            self.filter_by_exp.on('end', self.run_mir_pre)
        if self.option("knownlnc_fasta").is_set:
            self.filter_by_exp.on('end', self.run_new_annot)

        self.on_rely([self.annot_mapdb,self.annot_orfpfam], self.run_annot_filter_new)
        self.annot_filter_new.on('end', self.run_annot_class_new)
        self.on_rely([self.annot_class_new, self.annot_class_ref], self.run_merge_annot)
        '''
        if self.option("is_target_cis").lower() == "true":
            self.merge_annot.on('end', self.run_target_cis)
        '''
        if self.samples_num >= 3:
            self.merge_annot.on('end', self.run_target_cistrans)

        self.lncrna_identify.on('end', self.run_gene_fa)

        if self.option("sample_num") == "multiple":
            self.diffexpress.on('end', self.run_exp_corr)
            '''
            if self.option("is_target_trans").lower() == "true":
                self.on_rely([self.diffexpress_t, self.merge_annot], self.run_target_trans)
            '''
            self.filter_by_exp.on('end', self.run_diffexpress)
            if self.option("group_table").prop["sample_number"] > 2:
                self.diffexpress.on('end', self.run_exp_pca)
                if self.min_group_num >= 3:
                    self.exp_pca.on('end', self.run_ellipse)
                    self.exp_pca_t.on('end', self.run_ellipse_t)

        for i in self.final_tools:
            self.logger.info("final tool {}".format(i))
        self.on_rely(self.final_tools, self.end)



        if self.option("qc_dir").is_set:
            self.prepare_file()
        self.run_filecheck()
        super(LncRnaWorkflow, self).run()

    def prepare_file(self):
        fq_list = os.path.join(self.work_dir, "fq_list.txt")
        fq_dict = dict()
        fq_path = self.option("qc_dir").prop['path'] + "/"
        with open(fq_path + "list.txt", 'r') as qc_f:
            for line in qc_f:
                if len(line.split("\t")) == 3:
                    file_dir, sample, strand = line.strip().split("\t")
                    if sample in fq_dict:
                        fq_dict[sample].update({strand: file_dir})
                    else:
                        fq_dict[sample] = {strand: file_dir}
                elif len(line.split("\t")) == 2:
                    file_dir, sample = line.strip().split("\t")
                    fq_dict[sample] = {"l": file_dir}
        with open(fq_list, "w") as w:
            for sample in fq_dict:
                if "r" in fq_dict[sample]:
                    w.write("{}\t{}\t{}\n".format(sample, fq_path + fq_dict[sample]["l"], fq_path + fq_dict[sample]["r"]))
                else :
                    w.write("{}\t{}\n".format(sample, fq_path + fq_dict[sample]["l"]))
        self.qc_list = fq_list


    def run_new_lncrna(self):
        if self.genome_info["name"] in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
            species2cpat = {
                'Homo_sapiens': "Human",
                'Mus_musculus': "Mouse",
                'Danio_rerio': "Zebrafish",
                'Drosophila_melanogaster': "Fly"
            }
            hexamer_dat = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_Hexamer.tsv".format(species2cpat[self.genome_info["name"]])
            logit_model = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_logitModel.RData".format(species2cpat[self.genome_info["name"]])
        else:
            pass

        opts = {
            "cpc" : self.option("cpc").lower()== "true",
            "cnci" : self.option("cnci").lower()== "true",
            "cpat" : self.option("cpat").lower()== "true",
            "pfamscan" : self.option("pfamscan").lower()== "true",
            "new_fasta": os.path.join(self.assembly.output_dir, "NewTranscripts/new_transcripts.fa"),
            "new_gtf": os.path.join(self.assembly.output_dir, "NewTranscripts/new_transcripts.gtf"),
            "mrna_gtf": self.mrna_gtf,
            "lnc_db_gtf": self.lnc_gtf,
            "identify_num" : self.option("lnc_identify_num"),
            "transcript_len" : self.option("transcript_len"),
            "exon_num" : self.option("exon_num"),
            "orf_len" : self.option("orf_len"),
            "cpc_score" : self.option("cpc_score"),
            "cnci_score" : self.option("cnci_score"),
            "taxonmy" : self.option("taxonmy"),
            # "lnc_db_gtf" : self.lnc_gtf,
            "biomart" : self.des,
            "biomart_type" : self.des_type,
        }

        if self.option("cpat").lower()== "true":
            if self.genome_info["name"] in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
                hexamer_dat = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_Hexamer.tsv".format(species2cpat[self.genome_info["name"]])
                logit_model = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_logitModel.RData".format(species2cpat[self.genome_info["name"]])
                opts.update({
                    "cpat_score" : self.option("cpat_score"),
                    "hexamer_dat" : hexamer_dat,
                    "logit_model" : logit_model
                })
            else:
                pass

        self.lnc_new_opts = opts
        self.lncrna_new.set_options(opts)
        self.lncrna_new.on('end', self.set_output, 'lncrna_new')
        self.lncrna_new.on('start', self.set_step, {'start': self.step.lncrna_new})
        self.lncrna_new.on('end', self.set_step, {'end': self.step.lncrna_new})
        self.lncrna_new.run()


    def run_merge_new(self):
        opts = {
            "all_known_gtf": self.filecheck.option("all_known_gtf"),
            "all_known_fa": self.filecheck.option("all_known_fa"),
            "new_mrna_gtf": os.path.join(self.lncrna_new.output_dir,  'novel_mrna.gtf'),
            "new_mrna_fa":  os.path.join(self.lncrna_new.output_dir,  'novel_mrna.fa'),
            "new_lncrna_gtf":  os.path.join(self.lncrna_new.output_dir,  'novel_lncrna.gtf'),
            "new_lncrna_fa":  os.path.join(self.lncrna_new.output_dir,  'novel_lncrna.fa'),
        }
        self.merge_known_new.set_options(opts)
        self.merge_known_new.on('end', self.set_output, 'merge_known_new')
        self.merge_known_new.on('start', self.set_step, {'start': self.step.merge_known_new})
        self.merge_known_new.on('end', self.set_step, {'end': self.step.merge_known_new})
        self.merge_known_new.run()

    def run_new_annot(self):
        opts = {
            "reference": self.option("knownlnc_fasta"),
            "query": self.filter_by_exp.output_dir + "/filtered_lncnovel/novel_lncrna.fa",
            "evalue": self.option("knownlnc_evalue"),
            "qcov": self.option("knownlnc_qcov"),
            "scov": self.option("knownlnc_scov"),
        }
        self.lncrna_new_annot.set_options(opts)
        self.lncrna_new_annot.on('end', self.set_output, 'lncrna_new_annot')
        self.lncrna_new_annot.on('start', self.set_step, {'start': self.step.lncrna_new_annot})
        self.lncrna_new_annot.on('end', self.set_step, {'end': self.step.lncrna_new_annot})
        self.lncrna_new_annot.run()

    def run_lncrna_identify(self):
        opts = {
            "exp_matrix" : self.express.output_dir + "/transcript.tpm.matrix",
            "lnc_db_gtf" : self.lnc_gtf,
            "lnc_db_fa" : self.lnc_fa,
            "mrna_gtf" : self.mrna_gtf,
            "mrna_fa" : self.mrna_fa,
            "biomart" : self.des,
            "biomart_type" : self.des_type,
            "ids_mapping" : self.lnc_ids,
            "database" : "ensembl"
        }

        self.lncrna_identify.set_options(opts)
        self.lncrna_identify.on('end', self.set_output, 'lncrna_identify')
        self.lncrna_identify.on('start', self.set_step, {'start': self.step.lncrna_identify})
        self.lncrna_identify.on('end', self.set_step, {'end': self.step.lncrna_identify})
        self.lncrna_identify.run()


    def run_lncrna_stat(self):
        opts = {
            "novel_lncrna_detail": self.filter_by_exp.output_dir + "/filtered_lncnovel/novel_lncrna_predict_detail.xls",
            "known_lncrna_detail": self.lncrna_identify.output_dir + "/known_lncrna_detail.xls",
            "exp_matrix" : self.express.output_dir + "/transcript.tpm.matrix",
            "gene_type": self.filter_by_exp.output_dir + "/filtered_file/gene_type.xls",
            "new_gene_gtf" : self.assembly.output_dir + "/NewTranscripts/new_genes.gtf",
        }

        self.lncrna_stat.set_options(opts)
        self.lncrna_stat.on('end', self.set_output, 'lncrna_stat')
        self.lncrna_stat.on('start', self.set_step, {'start': self.step.lncrna_stat})
        self.lncrna_stat.on('end', self.set_step, {'end': self.step.lncrna_stat})
        self.lncrna_stat.run()

    def run_filter_by_exp(self):
        novel_dir = self.lncrna_new.output_dir
        # known_dir =  self.lnc_dir
        opts = {
            "novel_mrna_gtf" :  novel_dir + "/novel_mrna.gtf",
            "novel_mrna_fa" :  novel_dir + "/novel_mrna.fa",
            "novel_lncrna_gtf" :  novel_dir +"/novel_lncrna.gtf",
            "novel_lncrna_fa" :  novel_dir +"/novel_lncrna.fa",
            "known_mrna_gtf" : self.mrna_gtf ,
            "known_mrna_fa" :  self.mrna_fa,
            "known_lncrna_gtf" :  self.lnc_gtf,
            "known_lncrna_fa" :  self.lnc_fa,
            "tpm_matrix" :  self.express.output_dir + "/transcript.tpm.matrix",
            "count_matrix" :  self.express.output_dir + "/transcript.count.matrix",
            "tpm_matrix_g" :  self.express.output_dir + "/gene.tpm.matrix",
            "count_matrix_g" :  self.express.output_dir + "/gene.count.matrix",
            "lnc_new_dir" :  novel_dir,
        }
        if os.path.exists(self.express.output_dir + "/gene.fpkm.matrix"):
            opts.update({"fpkm_matrix_g" :  self.express.output_dir + "/gene.fpkm.matrix"})
        if os.path.exists(self.express.output_dir + "/transcript.fpkm.matrix"):
            opts.update({"fpkm_matrix" :  self.express.output_dir + "/transcript.fpkm.matrix"})
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        if self.known_ko != "" and os.path.exists(os.path.join(db_path, self.known_ko)):
            opts.update({
                "known_ko": os.path.join(db_path, self.known_ko)
            })
        self.filter_by_exp.set_options(opts)
        self.filter_by_exp.on('end', self.set_output, 'filter_by_exp')
        self.filter_by_exp.on('start', self.set_step, {'start': self.step.filter_by_exp})
        self.filter_by_exp.on('end', self.set_step, {'end': self.step.filter_by_exp})
        self.filter_by_exp.run()

    def run_mir_pre(self):
        opts = {
            'lncrna_fa': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/all_lncrna.fa'),
            'species': self.mirna_info["organism"],
            'known_list': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/known_lncrna_ids.list'),
            'novel_list': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/novel_lncrna_ids.list'),
            't2g': os.path.join(self.assembly.output_dir,  'NewTranscripts/trans2gene')
        }
        self.mir_pre.set_options(opts)
        self.mir_pre.on('end', self.set_output, 'mir_pre')
        self.mir_pre.on('start', self.set_step, {'start': self.step.mir_pre})
        self.mir_pre.on('end', self.set_step, {'end': self.step.mir_pre})
        self.mir_pre.run()

    def run_lnc_family(self):
        opts = {
            'lncrna_fa': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/all_lncrna.fa'),
            'known_list': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/known_lncrna_ids.list'),
            'novel_list': os.path.join(self.filter_by_exp.output_dir,  'filtered_file/novel_lncrna_ids.list'),
            't2g': os.path.join(self.assembly.output_dir,  'NewTranscripts/trans2gene')
        }
        self.lnc_family.set_options(opts)
        self.lnc_family.on('end', self.set_output, 'lnc_family')
        self.lnc_family.on('start', self.set_step, {'start': self.step.lnc_family})
        self.lnc_family.on('end', self.set_step, {'end': self.step.lnc_family})
        self.lnc_family.run()

    def run_gene_fa(self):
        if self.option("is_assemble") == True:
            opts = {
                "ref_new_gtf": self.merge_known_new.output_dir + '/known_and_new.gtf',
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

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'mrna_gtf': self.mrna_gtf,
            'lncrna_gtf': self.lnc_gtf,
            'mrna_fa': self.mrna_fa,
            'lncrna_fa': self.lnc_fa,
            "in_gtf": self.ref_gtf,
            'sample_num': self.option('sample_num'),
            'is_duplicate': self.option('is_duplicate')
        }
        '''
        if os.path.exists(self.mrna_fa) and os.path.exists(self.mrna_gtf):
            opts.update({
                'mrna_fa': self.mrna_fa,
                'mrna_gtf': self.mrna_gtf,
            })
        else:
            opts.update({
                'mrna_fa': self.mrna_fa,
                'mrna_gtf': self.mrna_gtf,
            })

        if os.path.exists(self.lnc_fa) and os.path.exists(self.lnc_gtf):
            opts.update({
                'lncrna_fa': self.lnc_fa,
                'lncrna_gtf': self.lnc_gtf,
            })
        else:
            pass
        '''


        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        if self.option("control_file").is_set:
            opts.update({'control_file': self.option('control_file')})
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        if self.option('qc_soft') == "fastp":
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
            if self.option("qc_dir").is_set:
                self.qc_stat_after.set_options({
                    'fastq_dir': self.option('qc_dir'),
                    'fq_type': self.option('fq_type'),
                    'quality': quality,
                    'rfam': True,
                    'dup': True
                })
            else:
                self.qc_stat_after.set_options({
                    'fastq_dir': self.qc.option('sickle_dir'),
                    'fq_type': self.option('fq_type'),
                    'quality': quality,
                    'rfam': True,
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
        if self.option("qc_dir").is_set:
            qc_dir = self.option("qc_dir")
        else:
            qc_dir = self.qc.option("sickle_dir")
        opts = {
            "ref_genome": self.genome_info["organism_name"],
            "genome_version": self.genome_info["assembly"], # 参考基因组版本
            "genome_annot_version": self.genome_info["annot_version"],  # 参考基因组注释版本
            "mapping_method": self.option("align_method").lower(),  # 比对软件
            "seq_method": self.option("fq_type"),   # PE or SE
            "fastq_dir": qc_dir,
            "assemble_method": self.option("assemble_method").lower(),
        }
        if self.option("strand_specific"):
            opts.update(
                {
                    "strand_specific":True,
                }
            )
            strand_dir = "firststrand"

            if self.option("strand_dir").lower() in ["rf", 'r', "forward"]:
                strand_dir = "firststrand"
            elif self.option("strand_dir").lower() in ["fr", 'f', "backward"]:
                strand_dir = "secondstrand"
            opts.update({
                "strand_direct": strand_dir,
            })
        self.mapping.set_options(opts)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_assembly(self):
        self.logger.info("开始运行拼接步骤")
        bam_locfile = os.path.join(self.mapping.output_dir, "bam_loc.txt")
        opts = {
            'bam_loc': bam_locfile,
            "assemble_method": self.option("assemble_method").lower(),
            "ref_gtf": self.filecheck.option("all_known_gtf"),
            "ref_fa": self.ref_genome,
        }
        # 如果具有链特异性
        if self.option("strand_specific"):
            strand_dir = "firststrand"
            if self.option("strand_dir").lower() in ["rf", 'r', "forward"]:
                strand_dir = "firststrand"
            elif self.option("strand_dir").lower() in ["fr", 'f', "backward"]:
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
            "keep_list": os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_mrna_ids.list"),
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
            "keep_list": os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_mrna_ids.list"),
        }
        self.annot_filter_new.set_options(options)
        self.annot_filter_new.on('start', self.set_step, {'start': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_step, {'end': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_output, "annot_filter_new")
        self.annot_filter_new.run()

    def run_annot_class_ref(self):
        ref_filter_dir = self.annot_filter_ref.output_dir
        options = {
            'taxonomy': self.taxon_dict[self.option("kegg_database")][1],
            'blast_nr_xml': ref_filter_dir + "/nr/blast.filter.xml",
            'blast_kegg_xml': ref_filter_dir + "/kegg/blast.filter.xml",
            'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.filter.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.filter.xml",
            'pfam_domain': ref_filter_dir + "/pfam/pfam.filter.xls",
            'blast2go_annot': ref_filter_dir + "/go/blast2go.filter.xls",
            'link_bgcolor': "yellow",
            'png_bgcolor': '#FFFF00',
            'g2t2p': self.g2t2p,
            'gtf': os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_mrna.gtf"),
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'type': "ref",
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
        }

        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        if self.known_ko != "" and os.path.exists(os.path.join(db_path, self.known_ko)):
            options.update({
                "known_ko":  os.path.join(self.filter_by_exp.output_dir,  'filtered_file/known_ko.xls'),
            })

        self.annot_class_ref.set_options(options)
        self.annot_class_ref.on('start', self.set_step, {'start': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_step, {'end': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_output, "annot_class_ref")
        self.annot_class_ref.run()

    def run_annot_class_new(self):
        new_filter_dir = self.annot_filter_new.output_dir
        options = {
            'taxonomy': self.taxon_dict[self.option("kegg_database")][1],
            'blast_nr_xml': new_filter_dir + "/nr/blast.filter.xml",
            'blast_kegg_xml': new_filter_dir + "/kegg/blast.filter.xml",
            'blast_eggnog_xml': new_filter_dir + "/eggnog/blast.filter.xml",
            'blast_swissprot_xml': new_filter_dir + "/swissprot/blast.filter.xml",
            'pfam_domain': new_filter_dir + "/pfam/pfam.filter.xls",
            'blast2go_annot': new_filter_dir + "/go/blast2go.filter.xls",
            'gene2trans': self.annot_orfpfam.output_dir + "/all_tran2gen.txt",
            'link_bgcolor': "green",
            'png_bgcolor': '#00CD00',
            'type': "new",
            'gtf': self.assembly.option("new_transcripts_gtf"),
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
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

    def run_classify_quant(self):
        lncrna_dir = os.path.join(self.filter_by_exp.output_dir, "filtered_file")
        quant_dir = self.express.output_dir
        opts = {
            "known_lnc_rna" : lncrna_dir + "/" + "known_lncrna_ids.list",
            "new_lnc_rna" : lncrna_dir + "/" + "novel_lncrna_ids.list",
            "known_mrna" : lncrna_dir + "/" + "known_mrna_ids.list",
            "new_mrna" : lncrna_dir + "/" + "novel_mrna_ids.list",
            "quant_results" : quant_dir,
            "t2g" : os.path.join(self.merge_known_new.output_dir, "t2g.txt")
        }
        self.classify_quant.set_options(opts)
        self.classify_quant.on("end", self.set_output, "classify_quant")
        self.classify_quant.on('start', self.set_step, {'start': self.step.classify_quant})
        self.classify_quant.on('end', self.set_step, {'end': self.step.classify_quant})
        self.classify_quant.run()

    def run_express(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option("strand_dir").lower() in ["fr", 'f', "backward"]:
                if self.option("fq_type") == "PE":
                    self.libtype = "fr"
                else:
                    self.libtype = "f"
            else:
                if self.option("fq_type") == "PE":
                    self.libtype = "rf"
                else:
                    self.libtype = "r"
        else:
            self.libtype = None
        if self.option("qc_dir").is_set:
            fq_list = self.qc_list
        else:
            fq_list = self.qc.option("fq_list").prop['path']

        if self.option("is_assemble") == True:
            transcriptome = os.path.join(self.merge_known_new.output_dir, "known_and_new.fa")
        else:
            transcriptome = os.path.join(self.merge_known_new.output_dir, "known_and_new.fa")

        opts = {
            "fastq" : fq_list,
            "method" : self.option("express_method"),
            "tpm_threshold" : self.option("filter_tpm"),
            "libtype" : self.libtype,
            "transcriptome" : transcriptome,
            "t2g" : os.path.join(self.merge_known_new.output_dir, "t2g.txt")
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
                "exp" : self.known_fpkm
            }
        else:
            opts = {
                "exp" : self.known_tpm
            }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.run()

        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : self.known_fpkm_t
            }
        else:
            opts = {
                "exp" : self.known_tpm_t
            }
        self.exp_pca_t.set_options(opts)
        self.exp_pca_t.on("end", self.set_output, "exp_pca_t")
        self.exp_pca_t.on('start', self.set_step, {'start': self.step.exp_pca_t})
        self.exp_pca_t.on('end', self.set_step, {'end': self.step.exp_pca_t})
        self.exp_pca_t.run()

    def run_ellipse(self):
        self.logger.info("开始运行pca ellipse")
        pc_map = {'pca':"/PCA.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        opts = {
            'analysis': "pca",
            'pc_table': self.exp_pca.output_dir + pc_map['pca'],
            'group_table':self.option('group_table').prop['path']
        }

        self.ellipse.set_options(opts)
        self.ellipse.on("end", self.set_output, "ellipse")
        self.ellipse.on('start', self.set_step, {'start': self.step.ellipse})
        self.ellipse.on('end', self.set_step, {'end': self.step.ellipse})
        self.ellipse.run()

    def run_ellipse_t(self):
        self.logger.info("开始运行pca ellipse")
        pc_map = {'pca':"/PCA.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        opts = {
            'analysis': "pca",
            'pc_table': self.exp_pca_t.output_dir + pc_map['pca'],
            'group_table':self.option('group_table').prop['path']
        }

        self.ellipse_t.set_options(opts)
        self.ellipse_t.on("end", self.set_output, "ellipse_t")
        self.ellipse_t.on('start', self.set_step, {'start': self.step.ellipse_t})
        self.ellipse_t.on('end', self.set_step, {'end': self.step.ellipse_t})
        self.ellipse_t.run()

    def run_exp_corr(self):
        self.logger.info("开始运行聚类分析")

        if os.path.exists(self.work_dir + "/gene.fpkm.matrix"):
            known_fpkm = self.work_dir + "/gene.fpkm.matrix"
        else:
            known_fpkm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/gene.tpm.matrix"):
            known_tpm = self.work_dir + "/gene.tpm.matrix"
        else:
            known_tpm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_tpm_known.xls")

        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : known_fpkm
            }
        else:
            opts = {
                "exp" : known_tpm
            }
        self.exp_corr.set_options(opts)
        self.exp_corr.on("end", self.set_output, "exp_corr")
        self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.run()

        novel_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_lncrna.gtf")

        if os.path.exists(self.work_dir + "/transcript.fpkm.matrix"):
            known_fpkm = self.work_dir + "/transcript.fpkm.matrix"
        else:
            known_fpkm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/transcript.tpm.matrix"):
            known_tpm = self.work_dir + "/transcript.tpm.matrix"
        else:
            known_tpm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_tpm_known.xls")

        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp" : known_fpkm
            }
        else:
            opts = {
                "exp" : known_tpm
            }

        self.exp_corr_t.set_options(opts)
        self.exp_corr_t.on("end", self.set_output, "exp_corr_t")
        self.exp_corr_t.on('start', self.set_step, {'start': self.step.exp_corr_t})
        self.exp_corr_t.on('end', self.set_step, {'end': self.step.exp_corr_t})
        self.exp_corr_t.run()

    def run_target_cistrans(self):
        known_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_lncrna.gtf")
        novel_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_lncrna.gtf")
        mrna_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/all_mrna.gtf")
        annotation =  os.path.join(self.merge_annot.output_dir , "allannot_class/all_annot.xls")


        options = {
            "known" : known_gtf,
            "novol" : novel_gtf,
            "mrna_gtf" : mrna_gtf,
            "annotation" : annotation,
            "up_dis": self.option("up_dis"),
            "down_dis": self.option("down_dis"),
            "cor_cutoff": float(self.option("target_cor_cutoff")),
            "corr_way": self.option("target_corr_way"),
            'padjust_way': self.option("target_padjust_way").lower(),
            'exp_matrix_lnc': os.path.join(self.classify_quant.output_dir, 'transcript.tpm.matrix'),
            'exp_matrix_target': os.path.join(self.classify_quant.output_dir, 'gene.tpm.matrix')
        }
        if self.option("target_pvalue_type") == "pvalue":
            options.update({'pvalue_cutoff': float(self.option("target_pvalue_cutoff"))})
        else:
            options.update({'qvalue_cutoff': float(self.option("target_pvalue_cutoff"))})
        self.target_cistrans.set_options(options)
        self.target_cistrans.on("end", self.set_output, "target_cistrans")
        self.target_cistrans.on('start', self.set_step, {'start': self.step.target_cistrans})
        self.target_cistrans.on('end', self.set_step, {'end': self.step.target_cistrans})
        self.target_cistrans.run()

    def run_target_cis(self):
        known_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_lncrna.gtf")
        novel_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_lncrna.gtf")
        mrna_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/all_mrna.gtf")
        annotation =  os.path.join(self.merge_annot.output_dir , "allannot_class/all_annot.xls")

        options = {
            "lncrna_gtf" : known_gtf,
            "mrna_gtf" : mrna_gtf,
            "annotation" : annotation,
            "up_dis": 10,
            "down_dis": 10,
        }
        self.target_cis_known.set_options(options)
        self.target_cis_known.on("end", self.set_output, "target_cis_known")
        self.target_cis_known.on('start', self.set_step, {'start': self.step.target_cis_known})
        self.target_cis_known.on('end', self.set_step, {'end': self.step.target_cis_known})
        self.target_cis_known.run()

        options = {
            "lncrna_gtf" : novel_gtf,
            "mrna_gtf" : mrna_gtf,
            "annotation" : annotation,
            "up_dis": 10,
            "down_dis": 10,
        }
        self.target_cis_novel.set_options(options)


        self.target_cis_novel.on("end", self.set_output, "target_cis_novel")
        self.target_cis_novel.on('start', self.set_step, {'start': self.step.target_cis_novel})
        self.target_cis_novel.on('end', self.set_step, {'end': self.step.target_cis_novel})
        self.target_cis_novel.run()

    def run_target_trans(self):
        known_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_lncrna.gtf")
        novel_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_lncrna.gtf")
        mrna_gtf = os.path.join(self.filter_by_exp.output_dir, "filtered_file/all_mrna.gtf")
        known_fa = os.path.join(self.filter_by_exp.output_dir, "filtered_file/known_lncrna.fa")
        novel_fa = os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_lncrna.fa")
        target_fa = os.path.join(self.filter_by_exp.output_dir, "filtered_file/all_mrna.fa")
        annotation = os.path.join(self.merge_annot.output_dir , "allannot_class/all_annot.xls")
        diff_summary = glob.glob(os.path.join(self.diffexpress_t.output_dir , "*_diff_summary.xls"))[0]
        options = {
            "target" : target_fa,
            "query": known_fa,
            "method": self.option("target_trans_method"),
            "gtf": mrna_gtf,
            "lnc_gtf": known_gtf,
            "annotation": annotation,
            "diff_summary": diff_summary,
            'query_list_type': "T",
            'target_list_type': "T"
        }
        self.target_trans_known.set_options(options)
        self.target_trans_known.on("end", self.set_output, "target_trans_known")
        self.target_trans_known.on('start', self.set_step, {'start': self.step.target_trans_known})
        self.target_trans_known.on('end', self.set_step, {'end': self.step.target_trans_known})
        self.target_trans_known.run()

        options = {
            "target" : target_fa,
            "query": novel_fa,
            "method": self.option("target_trans_method"),
            "gtf": mrna_gtf,
            "lnc_gtf": novel_gtf,
            "annotation": annotation,
            "diff_summary": diff_summary,
            'query_list_type': "T",
            'target_list_type': "T"
        }
        self.target_trans_novel.set_options(options)

        self.target_trans_novel.on("end", self.set_output, "target_trans_novel")
        self.target_trans_novel.on('start', self.set_step, {'start': self.step.target_trans_novel})
        self.target_trans_novel.on('end', self.set_step, {'end': self.step.target_trans_novel})
        self.target_trans_novel.run()


    def run_diffexpress(self):
        self.logger.info("开始运行基因差异表达分析")
        if os.path.exists(self.work_dir + "/gene.fpkm.matrix"):
            self.known_fpkm = self.work_dir + "/gene.fpkm.matrix"
        else:
            self.known_fpkm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/gene.tpm.matrix"):
            self.known_tpm = self.work_dir + "/gene.tpm.matrix"
        else:
            self.known_tpm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_tpm_known.xls")


        if os.path.exists(self.work_dir + "/transcript.fpkm.matrix"):
            self.known_fpkm_t = self.work_dir + "/transcript.fpkm.matrix"
        else:
            self.known_fpkm_t = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/transcript.tpm.matrix"):
            self.known_tpm_t = self.work_dir + "/transcript.tpm.matrix"
        else:
            self.known_tpm_t = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_tpm_known.xls")

        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            '''
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
            '''
            count_file = self.express.output_dir + "/gene.count.matrix"
            fpkm_file = self.express.output_dir + "/gene.fpkm.matrix"
            opts = {
                "count" : count_file,
                "exp" : fpkm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                # "pvalue_padjust" : self.option("pvalue_padjust"),
                # "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                # "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                opts.update({
                    'pvalue_padjust': self.option('pvalue_padjust'),
                    'pvalue': float(self.option('diff_fdr_ci')),
                    'padjust_way': self.option('padjust_way'),
                })
            if self.option('diff_method').lower() == 'noiseq':
                opts.update({
                    self.option('pvalue_padjust'): float(self.option('diff_fdr_ci')),
                })
        else:
            '''
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
            '''
            count_file = self.express.output_dir + "/gene.count.matrix"
            tpm_file = self.express.output_dir + "/gene.tpm.matrix"
            opts = {
                "count" : count_file,
                "exp" : tpm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                # "pvalue_padjust" : self.option("pvalue_padjust"),
                # "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                # "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                opts.update({
                    'pvalue_padjust': self.option('pvalue_padjust'),
                    'pvalue': float(self.option('diff_fdr_ci')),
                    'padjust_way': self.option('padjust_way'),
                })
            if self.option('diff_method').lower() == 'noiseq':
                opts.update({
                    self.option('pvalue_padjust'): float(self.option('diff_fdr_ci')),
                })
            if self.option("filter_method") != "None" and self.option("tpm_filter_threshold") != "NA":
                opts["filter_method"]=self.option("filter_method")
                opts["tpm_filter_threshold"]=self.option("tpm_filter_threshold")
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            count_file = self.express.output_dir + "/transcript.count.matrix"
            fpkm_file = self.express.output_dir + "/transcript.fpkm.matrix"
            opts = {
                "count" : count_file,
                "exp" : fpkm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                # "pvalue_padjust" : self.option("pvalue_padjust"),
                # "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                # "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                opts.update({
                    'pvalue_padjust': self.option('pvalue_padjust'),
                    'pvalue': float(self.option('diff_fdr_ci')),
                    'padjust_way': self.option('padjust_way'),
                })
            if self.option('diff_method').lower() == 'noiseq':
                opts.update({
                    self.option('pvalue_padjust'): float(self.option('diff_fdr_ci')),
                })
        else:
            count_file = self.express.output_dir + "/transcript.count.matrix"
            tpm_file = self.express.output_dir + "/transcript.tpm.matrix"
            opts = {
                "count" : count_file,
                "exp" : tpm_file,
                "group" : self.option("group_table"),
                "cmp" : self.option("control_file"),
                # "pvalue_padjust" : self.option("pvalue_padjust"),
                # "pvalue" : float(self.option("diff_fdr_ci")),
                "fc" : float(self.option("fc")),
                # "padjust_way" : self.option("padjust_way"),
                "method" : self.option("diff_method"),
                "exp_type": self.option("exp_way")
            }
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                opts.update({
                    'pvalue_padjust': self.option('pvalue_padjust'),
                    'pvalue': float(self.option('diff_fdr_ci')),
                    'padjust_way': self.option('padjust_way'),
                })
            if self.option('diff_method').lower() == 'noiseq':
                opts.update({
                    self.option('pvalue_padjust'): float(self.option('diff_fdr_ci')),
                })
            if self.option("filter_method") != "None"  and self.option("tpm_filter_threshold") != "NA":
                opts["filter_method"]=self.option("filter_method")
                opts["tpm_filter_threshold"]=self.option("tpm_filter_threshold")
        self.diffexpress_t.set_options(opts)
        self.diffexpress_t.on("end", self.set_output, "diffexpress_t")
        self.diffexpress_t.on('start', self.set_step, {'start': self.step.diffexpress_t})
        self.diffexpress_t.on('end', self.set_step, {'end': self.step.diffexpress_t})
        self.diffexpress_t.run()


    def run_annot_mapdb(self, event):
        self.logger.info("开始运行diamond注释")
        self.logger.info("开始运行diamond 数据库为{}".format(self.option("nr_database")))

        opts = {
            "query" : os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_mrna.fa"),
            "nr_db" : self.taxon_dict[self.option("nr_database")][0],
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version']

        }
        self.annot_mapdb.set_options(opts)
        self.annot_mapdb.on("end", self.set_output, "annot_mapdb")
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.logger.info("开始运行pfam注释")
        opts = {
            "fasta" : os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_mrna.fa"),
            "gtf" : os.path.join(self.filter_by_exp.output_dir, "filtered_file/novel_mrna.gtf"),
            "pfam_version": self.annot_config_dict['pfam']['version'],
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
        options.update({
            "annot_class_new": self.annot_class_new.output_dir,
            "annot_db": self.annot_mapdb.output_dir
        })
        self.merge_annot.set_options(options)
        self.merge_annot.on('start', self.set_step, {'start': self.step.merge_annot})
        self.merge_annot.on('end', self.set_step, {'end': self.step.merge_annot})
        self.merge_annot.on('end', self.set_output, "merge_annot")
        self.merge_annot.run()

    def run_snp(self):
        self.logger.info("开始运行snp步骤")
        # "ref_gtf": self.filecheck.option("all_known_gtf"),

        opts = {
            "ref_fa": self.ref_genome,
            "ref_genome":  self.species_name,
            "ref_gtf": os.path.join(self.filter_by_exp.output_dir, "filtered_file/all.gtf"),
            "bam_list": self.mapping.option("bamlist"),
            'des': self.des,
            'des_type': self.des_type,
            'method_type': self.option("snp_method").lower()
        }
        self.snp.set_options(opts)
        self.snp.on("start", self.set_step, {"start": self.step.snp_rna})
        self.snp.on("end", self.set_step, {"end": self.step.snp_rna})
        self.snp.on("end", self.set_output, "snp")
        self.snp.run()

    def run_map_assess(self):
        opts = {
            "bam": self.mapping.option("bam_output"),
            "bed": self.filecheck.option("bed").prop['path'],
            "analysis": self.option("map_assess_method")
        }
        self.map_qc.set_options(opts)
        self.map_qc.on("start", self.set_step, {"start": self.step.map_qc})
        self.map_qc.on("end", self.set_step, {"end": self.step.map_qc})
        self.map_qc.on("end", self.set_output, "map_qc")
        self.map_qc.run()

    def run_altersplicing(self):
        if self.option("strand_specific"):
            lib_type = "fr-firststrand"
        else:
            lib_type = "fr-unstranded"
        # "ref_gtf": self.filecheck.option("all_known_gtf"),
        opts = {
            "bam_loc": os.path.join(self.mapping.output_dir, "bam_loc.txt"),
            "lib_type": lib_type,
            "ref_gtf": os.path.join(self.filter_by_exp.output_dir, "filtered_file/all.gtf"),
            "group_table": self.option("group_table"),
            "control_table": self.option("control_file")
        }
        if self.option("fq_type") == "PE":
            opts.update({"read_type": "paired"})
        else:
            opts.update({"read_type": "single"})
        self.altersplicing.set_options(opts)
        self.altersplicing.on("start", self.set_step, {"start": self.step.altersplicing})
        self.altersplicing.on("end", self.set_step, {"end": self.step.altersplicing})
        self.altersplicing.on("end", self.set_output, "altersplicing")
        self.altersplicing.run()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code= "13700319")
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
        elif event['data'] == 'qc_stat_before':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
            self.logger.info("开始设置qc的输出目录")
        elif event['data'] == 'qc_stat_after':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
            # self.check_result('qc_stat_after')
        elif event['data'] == 'transcripts':
            self.move2outputdir(obj.output_dir, 'abstract_transcripts')
        else:
            self.move2outputdir(obj.output_dir, event['data'])
            # .check_result(event['data'])

    def run_check(self):
        is_qc_right = self.check_qc()
        is_mapping_right = self.check_mapping()
        if is_qc_right and is_mapping_right:
            pass
        else:
            error_log = "异常结束 "
            if is_qc_right == False:
                error_log += "rrna含量异常，不符合参数设置"
            if is_mapping_right == False:
                error_log += "比对率异常，不符合参数设置"
            self.set_error_task(error_log)
            self.end_step = "mapping"
            self.end()

    def set_error_task(self, error_log=""):
        db = Config().get_mongo_client(mtype="project")[Config().get_mongo_dbname("project")]
        email = db['sg_task_email']
        email.update({"task_id": self.task_id}, {"$set": {"status": 5, "error_log": error_log}},
                     upsert=True)


    def check_qc(self):
        self.logger.info("检查样本核糖体占比")
        rrna_files = os.path.join(self.qc_stat_after.output_dir, 'stat_results')
        sample_num = self.option("group_table").prop["sample_number"]
        rrna_failed = 0
        with open(rrna_files, "r") as rrna_r:
            _ = rrna_r.readline()
            for line in rrna_r:
                rrna_ratio = line.strip().split('\t')[-1]
                self.logger.info("占比 {} 阈值 {}".format(rrna_ratio, self.option('rrna_ratio')))
                if float(rrna_ratio) > float(self.option('rrna_ratio')):
                    # self.logger.info("样本 {} 比对率不满足".format(rrna_r.readline().strip().split('\t')[0]))
                    rrna_failed += 1
        if self.option('rrna_stop') == "True" and float(rrna_failed)/sample_num*100 > float(self.option('rrna_sample_percent')):
            self.logger.info("样本 {} 比对率不满足条件， 程序中断".format(rrna_failed))
            return False
        else:
            return True

    def check_mapping(self):
        mapping_list = []
        unmap_list = []
        with open(os.path.join(self.output_dir, "mapping/Comparison_results"), "r") as f:
            f.readline()
            for line in f:
                cols = line.split("\t")
                mapped = float(cols[2].split("(")[0])
                total = float(cols[1])
                mapping_list.append(cols[0])
                if mapped/total < float(self.option("mapping_ratio"))/100:
                    unmap_list.append(cols[0])
        if self.option("mapping_stop") == "True" and float(len(unmap_list))/float(len(mapping_list)) > float(self.option("mapping_sample_percent"))/100:
            self.logger.info("样本 {} 比对率不满足条件， 程序中断".format(unmap_list))
            return False
        else:
            return True

    def check_result(self, step):
        if step == "qc_stat_after":
            if not self.check_qc():
                self.end_step = step
                # self.end()
        elif step == "mapping":
            if not self.check_mapping():
                self.end_step = step
                # self.end()

        else:
            pass

    def end(self, step=None):
        self.run_api()


        self.logger.info("更新主表 sg_task ")
        self.logger.info("更新主表 step {}".format(step))

        read_type_mapping = {
            "PE": "paired",
            "SE": "single"
        }
        lib_type_mapping = {
            "forward": "firststrand",
            "backward": "secondstrand",
            "RF": "firststrand",
            "FR": "secondstrand",
            "R": "firststrand",
            "F": "secondstrand",
            "U": "unstranded"
        }

        db = Config().get_mongo_client(mtype="lnc_rna")[Config().get_mongo_dbname("lnc_rna")]
        sg_task_col = db["sg_task"]
        if self.end_step in ["qc_stat_after", "mapping"]:
            task_updata_dict = {
                "ref_gtf" : self.ref_gtf,
                "ref_genome" : self.ref_genome,
                "genome_id" : self.genome_id,
                "organism_name" : self.genome_info["organism_name"],
                "annot_version" : self.genome_info["annot_version"],
                "assembly" : self.genome_info["assembly"],
                "des" : self.des,
                "des_type" : self.des_type,
                "read_type" : read_type_mapping[self.option("fq_type")],
                "lib_type" : lib_type_mapping[self.option("strand_dir")],
                "fastq" : os.path.join(self.workflow_output, "02 Basic_Analysis/01 QC/fq_list.txt"),
                "ref_fa_fai" : self.ref_genome + ".fai",
                "version": "v1.1",
                "samples_num": self.samples_num
            }
            self.logger.info("更新主表 sg_task {}".format(task_updata_dict))
            sg_task_col.update({"task_id" : self.task_id}, {"$set": task_updata_dict}, upsert=True)
            pass
        else:
            self.logger.info("更新主表 sg_task ")
            ## 更新一系列主表的字段，用于页面交互分析
            task_updata_dict = {
                "ref_gtf" : self.ref_gtf,
                "ref_genome" : self.ref_genome,
                "genome_id" : self.genome_id,
                "organism_name" : self.genome_info["organism_name"],
                "annot_version" : self.genome_info["annot_version"],
                "assembly" : self.genome_info["assembly"],
                "des" : self.des,
                "des_type" : self.des_type,
                "read_type" : read_type_mapping[self.option("fq_type")],
                "lib_type" : lib_type_mapping[self.option("strand_dir")],
                "fastq" : os.path.join(self.workflow_output, "02 Basic_Analysis/01 QC/fq_list.txt"),
                "known_lnc_gtf" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/known_lncrna.gtf'),
                "novel_lnc_gtf" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/novel_lncrna.gtf'),
                "all_lnc_gtf" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/all_lncrna.gtf'),
                "all_lnc_type" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/trans_type.xls'),
                "all_gene_type" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/gene_type.xls'),
                "all_trans_type" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/trans_type.xls'),
                "known_lnc_fa" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/known_lncrna.fa'),
                "novel_lnc_fa" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/novel_lncrna.fa'),
                "all_lnc_fa" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/all_lncrna.fa'),
                "assemble_fa" : os.path.join(self.workflow_output, '02 Basic_Analysis/03 Assemble/all_transcripts.fa'),
                "annot" : os.path.join(self.workflow_output, '03 Annotation/allannot_class/all_annot.xls'),
                "refrna_seqdb" : os.path.join(self.workflow_output, '11 Other/Sequence_database/refrna_seqs.db'),
                "assemble_t2g" : os.path.join(self.workflow_output, '02 Basic_Analysis/03 Assemble/trans2gene'),
                "mrna_fa" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/all_mrna.fa'),
                "mrna_gtf" : os.path.join(self.workflow_output, '11 Other/Filtered_result/filtered_file/all_mrna.gtf'),
                "ref_and_new_gtf" : os.path.join(self.workflow_output, '02 Basic_Analysis/03 Assemble/ref_and_new.gtf'),
                "ref_fa_fai" : self.ref_genome + ".fai",
                "version": "v1.1",
                "samples_num": self.samples_num
            }
            self.logger.info("更新主表 sg_task {}".format(task_updata_dict))

            sg_task_col.update({"task_id" : self.task_id}, {"$set": task_updata_dict}, upsert=True)
            annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
            if self.annot_config_dict['kegg']['version'] > "2020":
                pass
                # if self.option("kegg_org") not in [None, ""]:
                #     annot_version_dict['kegg'] += "_spe"
            else:
                del annot_version_dict['kegg']
            sg_task_col.update({'task_id': self.task_id},
                               {'$set': {'database_version': annot_version_dict,
                                         'annot_group': self.option("annot_group")}}, upsert=True)
            col1 = db["sg_annotation_stat"]
            col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": os.path.join(self.workflow_output, "03 Annotation")}}, upsert=True)
            col2 = db["sg_exp"]
            col2.update({"task_id" : self.task_id, "exp_level" : "G"}, {"$set": {"count_file": os.path.join(self.workflow_output, "05 Express/01 ExpAnnalysis/gene.count.matrix.xls")}}, upsert=True)
            col2.update({"task_id" : self.task_id, "exp_level" : "T"}, {"$set": {"count_file": os.path.join(self.workflow_output, "05 Express/01 ExpAnnalysis/transcript.count.matrix.xls")}}, upsert=True)
            if self.option("sample_num") == "multiple":
                if self.option("is_as").lower() == "true":
                    col3 = db["sg_splicing_rmats"]
                    for output_file in os.listdir(self.altersplicing.output_dir):
                        test_group, ctrl_group = output_file.split('_vs_')
                        # compare_plan = '{}|{}'.format(ctrl_group, test_group)
                        compare_plan = '{}|{}'.format(test_group, ctrl_group)
                        col3.update({
                            "task_id" : self.task_id,
                            "compare_plan": compare_plan},
                            {"$set": {"result_dir": os.path.join(self.workflow_output, '09 AS', output_file)}},
                            upsert=True
                        )
            self.merge_annotation_exp_matrix() # 表达量表增加注释信息
            if self.option("sample_num") == "multiple":
                self.merge_annotation_diffexp_matrix() # 差异表达量表增加注释信息

            self.modify_output()
        super(LncRnaWorkflow, self).end()

    def modify_output(self):
        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"

        tool_dir_dict = {
            'Sequence_database': '11 Other/Sequence_database',
            'QC_stat/before_qc/fastq_stat.xls': '02 Basic_Analysis/01 QC/rawdata_statistics.xls',
            'QC_stat/after_qc/fastq_stat.xls': '02 Basic_Analysis/01 QC/cleandata_statistics.xls',
            'QC_stat/sickle_dir/fq_list.txt': '02 Basic_Analysis/01 QC/fq_list.txt',
            'mapping/bam': '02 Basic_Analysis/02 Align/AlignBam',
            'mapping/stat': '02 Basic_Analysis/02 Align/AlignStat',
            'map_qc/chr_stat': '02 Basic_Analysis/02 Align/QualityAssessment',
            'map_qc/coverage': '02 Basic_Analysis/02 Align/QualityAssessment',
            'map_qc/distribution': '02 Basic_Analysis/02 Align/QualityAssessment',
            'map_qc/saturation': '02 Basic_Analysis/02 Align/QualityAssessment',
            'assembly/NewTranscripts': '02 Basic_Analysis/03 Assemble',
            'merge_annot': '03 Annotation',
            'classify_quant': '05 Express/01 ExpAnnalysis',
            'diffexpress': '06 Diff_Express',
            'diffexpress_t': '06 Diff_Express',
            'exp_corr/sample_correlation.xls': '05 Express/02 ExpCorr/G_sample_correlation.xls',
            'exp_pca': '05 Express/03 ExpPCA',
            'exp_corr_t/sample_correlation.xls': '05 Express/02 ExpCorr/T_sample_correlation.xls',
            'exp_pca_t': '05 Express/03 ExpPCA',
            'filter_by_exp': '11 Other/Filtered_result',
            'lncrna_identify': '04 LncRNA_Analysis/01 Known_LncRNA',
            'lncrna_new': '04 LncRNA_Analysis/02 Novel_LncRNA',
            'lncrna_stat': '04 LncRNA_Analysis/03 LncRNA_stat',
            'mir_pre': '08 LncRNA_Structure/02 miRNA_Pre',
            'lnc_family': '08 LncRNA_Structure/01 LncRNA_Family',
            '../snp_tmp': '10 SNP_InDel',
            'altersplicing': '09 AS',
            'target_cistrans': '07 LncRNA_Target'
        }

        '''
        'target_cis_novel/lnc_rna_cistarget.annot.xls': 'LncRNA_target_cis/novel_target.xls',
        'target_cis_known/lnc_rna_cistarget.annot.xls': 'LncRNA_target_cis/known_target.xls',
        'target_trans_novel/all_merge_out.annot.xls': 'LncRNA_target_trans/novel_target.xls',
        'target_trans_known/all_merge_out.annot.xls': 'LncRNA_target_trans/known_target.xls',
        '''
        if self.option("qc_dir").is_set:
            tool_dir_dict.update({"../fq_list.txt": '02 Basic_Analysis/01 QC/fq_list.txt'})

        rename_list = [
            ('02 Basic_Analysis/03 Assemble/trans2gene', ['trans2gene', 'trans2gene.txt']),
            ('02 Basic_Analysis/02 Align/AlignStat/*.stat', ['.stat', '_align_stat.txt']),
            ('02 Basic_Analysis/02 Align/QualityAssessment/*bam_chr_stat.xls',
             ['bam_chr_stat.xls', 'chr_distribution.xls']),
            ('02 Basic_Analysis/02 Align/QualityAssessment/coverage_*', ['coverage_', '']),
            ('02 Basic_Analysis/02 Align/QualityAssessment/*.reads_distribution.txt',
             ['.reads_distribution.txt', '.region_distribution.xls']),
            ('02 Basic_Analysis/02 Align/QualityAssessment/*.eRPKM.xls.cluster_percent.xls',
             ['.eRPKM.xls.cluster_percent.xls', '.region_distribution.xls']),
            ('02 Basic_Analysis/02 Align/QualityAssessment/satur_*', ['satur_', '']),
            ('05 Express/01 ExpAnnalysis/*.matrix', ['matrix', 'matrix.xls']),
            ('08 LncRNA_Structure/02 miRNA_Pre/*.tabular', ['.tabular', '.xls']),
            ('08 LncRNA_Structure/01 LncRNA_Family/*.tabular', ['.tabular', '.xls']),
        ]

        remove_list = ["03 Annotation/*/*/*/*.html.mark",
                       "03 Annotation/*/*/*/*.KOs.txt",
                       "03 Annotation/*/*/*/*.pdf",
                       "10 SNP_InDel/snp_annotation.xls",
                       "10 SNP_InDel/data_anno_pre.xls",
                       '02 Basic_Analysis/02 Align/QualityAssessment/*.R',
                       '02 Basic_Analysis/03 Assemble/old_*',
                       '05 Express/01 ExpAnnalysis/alignment_rate.txt',
                       '05 Express/01 ExpAnnalysis/trans_ids.list'
                       ]
        if self.option("is_as").lower() == "true":
            # move AS files generated while running api.database
            as_files = glob.glob(os.path.join(self.altersplicing.output_dir, '*/*splicing*stats*'))
            for item in as_files:
                folders = item.split('/')
                newpath = os.path.join(origin_dir, 'altersplicing', folders[-2], folders[-1])
                CopyFile().linkdir_test(item, newpath)

        for k, v in tool_dir_dict.items():
            CopyFile().linkdir_test(os.path.join(origin_dir, k), os.path.join(target_dir, v))
        for rename in rename_list:
            CopyFile().renamefile(os.path.join(target_dir, rename[0]), rename[1])
        for files in remove_list:
            CopyFile().remove(os.path.join(target_dir, files))

        CopyFile().linkdir_test(self.genome_stat, target_dir + "/01 Background/genome_info.xls")




        # Software Info
        software_info = os.path.join(target_dir, '01 Background', "software_info.xls")
        db = Config().get_mongo_client(mtype='lnc_rna', dydb_forbid=True)[Config().get_mongo_dbname('lnc_rna', dydb_forbid=True)]
        my_collection = db['sg_software_database']
        my_results = my_collection.find({})
        with open(software_info, "w") as w:
            w.write("\t".join(["Soft/Database", "Version", "Analysis", "Source"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["software_database"]), str(collection["version"]), str(collection["usage"]),
                     str(collection["source"])]) + "\n")

        # Sample info
        sample_info = os.path.join(target_dir, '01 Background', "sample_info.xls")
        if self.option('sample_num') == 'multiple':
            group_table = self.option('group_table').prop["path"]
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
        productive_names = {}
        mj_names = {}
        if self.option("productive_table").is_set:
            with open(self.option("productive_table").path, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    items = line.strip().split("\t")
                    if len(items) >= 2:
                        productive_names[items[0]] = items[1]
                    if len(items) >= 3:
                        mj_names[items[0]] = items[2]
        with open(sample_info, "w") as w, open(group_table, "r") as f:
            if self.option("productive_table").is_set:
                if mj_names:
                    w.write(
                        "\t".join(["Species Latin Name", "Sample Productive Name", "MJ_name", "Sample Name",
                                   "Group Name"]) + "\n")
                else:
                    w.write(
                        "\t".join(["Species Latin Name", "Sample Productive Name", "Sample Name", "Group Name"]) + "\n")
            else:
                w.write("\t".join(["Species Latin Name", "Sample Name", "Group Name"]) + "\n")
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    w.write("\t".join(
                        [self.option("ref_genome"), line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")

        if self.end_step in ["qc_stat_after", "mapping"]:
            pass
        else:
            if self.option('is_assemble'):
                os.system("cat {} {} > {}".format(
                    self.known_pep,
                    os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.pep'),
                    os.path.join(target_dir, '02\ Basic_Analysis/03\ Assemble/all_pep.fa')
                ))
                os.system("cat {} {} > {}".format(
                    self.known_cds,
                    os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.cds'),
                    os.path.join(target_dir, '02\ Basic_Analysis/03\ Assemble/all_cds.fa')
                ))
                with open(os.path.join(target_dir, '02 Basic_Analysis/03 Assemble/all_id.xls'), 'w') as fo, \
                        open(os.path.join(self.merge_annot.output_dir, "allannot_class/all_tran2gene.txt"), 'r') as fin:
                    fo.write("gene_id\ttranscript_id\tprotein_id\n")
                    for line in fin:
                        cols = line.strip("\n").split("\t")
                        fo.write("{}\t{}\t{}\n".format(cols[1], cols[0], cols[4]))

            else:
                os.link(self.known_cds,
                        os.path.join(target_dir, '02 Basic_Analysis/03 Assemble/all_cds.fa'))
                os.link(self.known_pep,
                        os.path.join(target_dir, '02 Basic_Analysis/03 Assemble/all_pep.fa'))

                with open(os.path.join(target_dir, '02 Basic_Analysis/03 Assemble/all_id.xls'), 'w') as fo, \
                        open(os.path.join(self.merge_annot.output_dir, "refannot_class/all_tran2gene.txt"), 'r') as fin:
                    fo.write("gene_id\ttranscript_id\tprotein_id\n")
                    for line in fin:
                        cols = line.strip("\n").split("\t")
                        fo.write("{}\t{}\t{}\n".format(cols[1], cols[0], cols[4]))

        sdir = self.add_upload_dir(target_dir)
        sdir.add_relpath_rules([[".", "", "工作流结果目录", 0]])

        file_reg = LncFileDes()

        for out_dir in ["01 Background", "02 Basic_Analysis", "03 Annotation", "04 LncRNA_Analysis", "05 Express",
                        "06 Diff_Express", "07 LncRNA_Target", "08 LncRNA_Structure", "09 AS", "10 SNP_InDel", "11 Other"]:
            dir_word = out_dir.split(" ")[1]
            if hasattr(file_reg, dir_word):
                repaths = getattr(file_reg, dir_word)
                for path in repaths:
                    if path[0] == "" or path[0] == ".":
                        path[0] = out_dir
                    else:
                        path[0] = out_dir + "/" + path[0]
                sdir.add_regexp_rules(repaths)
        print sdir._regexp_rules

    def run_api(self, step=None):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.lnc_rna')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.stop_timeout_check()
        self.export_genome_info()

        if self.end_step == "qc_stat_after" or self.end_step == "mapping":
            self.export_qc()
            self.export_productive_name()
            self.export_map_assess()
            self.logger.info("导表完成")
        else:
            self.export_qc()
            self.export_productive_name()
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
            self.export_lncrna()
            self.export_lncfamily()
            if self.mirpre:
                self.export_mirpre()
            if self.samples_num >= 3:
                self.export_target()
            self.export_gene_detail()
            self.logger.info("导表完成")



    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("lnc_rna.genome_info")
        file_path = self.genome_stat
        species_name = self.species_name
        species = self.species
        ref_anno_version = self.ref_anno_version
        hyperlink = self.hyperlink
        self.api_geno.add_genome_info(file_path=file_path,species_name=species_name, species=species, ref_anno_version=ref_anno_version,hyperlink=hyperlink)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("lnc_rna.lnc_rna_qc")
        qc_stat = self.qc_stat_before.output_dir
        fq_type = self.option("fq_type").lower()
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before")
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat"
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat"  # 将qc前导表加于该处
        if self.option("fastq_dir").is_set:
            self.api_qc.add_gragh_info(quality_stat_before, "before")
        self.api_qc.add_bam_path(self.mapping.output_dir + "/bam")
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

    def export_productive_name(self):
        api = self.api.api('lnc_rna.lnc_rna_qc')
        if self.option('productive_table').is_set:
            api.add_productive_name(samples=self.option('group_table').prop['sample'],
                                    productive_table=self.option('productive_table').path)
        else:
            pass

    @time_count
    def export_assembly(self):
        gevent.sleep()
        self.api_assembly = self.api.api("lnc_rna.assemble")
        if self.option("assemble_method").lower() == "cufflinks":
            all_gtf_path = self.assembly.output_dir + "/Cufflinks"
            merged_path = self.assembly.output_dir + "/Cuffmerge"
        else:
            all_gtf_path = self.assembly.output_dir + "/Stringtie"
            merged_path = self.assembly.output_dir + "/StringtieMerge"
        params = self.option("assemble_method").lower()
        self.api_assembly.add_assemble_result(all_gtf_path=all_gtf_path, params=params, merged_path=merged_path, statistics_path=self.assembly.output_dir + "/Statistics")


    @time_count
    def export_map_assess(self, only_mapping_rate=False):
        gevent.sleep()
        self.api_map = self.api.api("lnc_rna.lnc_rna_qc")
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
        self.api_annotation = self.api.api("lnc_rna.annotation")
        annot_dir = self.merge_annot.output_dir
        if self.option("is_assemble") == False:
            self.api_annotation.has_new = False
            trans2gene  = None
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        else:
            trans2gene = annot_dir + "/newannot_class/all_tran2gene.txt"
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        self.api_annotation.species_name = self.option("ref_genome")
        params = {
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
        query_id = self.api_annotation.run(annotation_dir=annot_dir,
                                trans2gene_new=trans2gene,
                                trans2gene_ref=trans2gene_ref,
                                params_dict=params,
                                taxon=self.option("kegg_database"),
                                exp_level=self.option("level").lower())

        self.api_annotation_lnc = self.api.api("lnc_rna.supplement")
        self.api_annotation_lnc.perfect_annotation_query_detail(
            query_id=query_id,
            known_tsv=os.path.join(self.lncrna_identify.output_dir, "known_lncrna_detail.xls"),
            novel_tsv=os.path.join(self.filter_by_exp.output_dir, "filtered_lncnovel/novel_lncrna_predict_detail.xls"),
            mrna_list=os.path.join(self.filter_by_exp.output_dir, "filtered_file/all_mrna_ids.list"),
            t2g_file=os.path.join(self.merge_known_new.output_dir, "t2g.txt"),
            entrez=self.entrez,
        )

    @time_count
    def export_target(self):
        self.lnc_target_api = self.api.api("lnc_rna.lnc_target")
        if self.option("group_table").is_set:
            group_dict = self.option('group_table').prop['group_dict']
            group_id = self.group_id
        else:
            group_dict = None
            group_id = None

        params_dict = {
            "up_dis": "10",
            "down_dis": "10",
            'group_dict':  group_dict,
            'group_id': str(group_id),
            "submit_location":"target_cistrans",
            'corr_cutoff': "0.9",
            'corr_way': "pearson",
            'padjust_way': "bh",
            "pvalue_type": "qvalue",
            "pvalue_cutoff": "0.05",
            "task_id": self.task_id,
            "task_type":2
        }

        self.lnc_target_api.import_target_cistrans(
            params_dict,
            os.path.join(self.target_cistrans.output_dir + '/cis_annot.xls'),
            os.path.join(self.target_cistrans.output_dir + '/trans_annot.xls')
        )


        '''
        params_cis = {
            "up_dis": "10",
            "down_dis": "10",
            "submit_location":"targe",
            "task_id": self.task_id,
            "task_type":2
        }
        self.lnc_target_api.species_name = self.ref_genome

        params_trans = {
            "lnc_geneset_id":"All",
            "m_geneset_id":"All",
            "method": self.option("target_trans_method").lower(),
            "submit_location":"trans",
            "task_id": self.task_id,
            "task_type":2

        }
        new_target_cis = self.target_cis_novel.output_dir + '/lnc_rna_cistarget.annot.xls'
        known_target_cis = self.target_cis_known.output_dir + '/lnc_rna_cistarget.annot.xls'

        new_target_trans = self.target_trans_novel.output_dir + '/all_merge_out.annot.xls'
        known_target_trans = self.target_trans_known.output_dir + '/all_merge_out.annot.xls'

        if self.option("is_target_cis").lower() == "true":
            self.lnc_target_api.import_target_cis(new_target_cis, known_target_cis, params_cis)
        if self.option("is_target_trans").lower() == "true" and self.option("sample_num") == "multiple":
            self.lnc_target_api.import_target_trans(new_target_trans, known_target_trans, params_trans)
        '''

    @time_count
    def export_lncrna(self):
        self.lncrna_api = self.api.api("lnc_rna.lnc_identify")
        params = self.lnc_new_opts

        if self.has_known_lnc:
            self.lncrna_api.known_lncrna_info(
                os.path.join(self.lncrna_identify.output_dir, "known_lncrna_detail.xls"))
        else:
            pass

        soft_list = [lnc_soft for lnc_soft in ['cpc', 'cnci', 'cpat', 'pfamscan'] if params[lnc_soft]]
        if self.option("knownlnc_fasta").is_set:
            annots = glob.glob(self.lncrna_new_annot.output_dir + '/*.xls')
            lnc_annot_path = annots[0]
        else:
            lnc_annot_path = None
        self.lncrna_api.new_lncrna_predict(
            os.path.join(self.filter_by_exp.output_dir, "filtered_lncnovel/novel_lncrna_predict_detail.xls"),
            os.path.join(self.filter_by_exp.output_dir, "filtered_lncnovel/novel_lncrna_stat.json"),
            tools = ",".join(soft_list).replace('pfamscan', 'pfam'),
            params = self.lnc_new_opts,
            lnc_annot_path = lnc_annot_path
        )


        self.lncrna_api.lncrna_stat(
            os.path.join(self.lncrna_stat.output_dir, 'lncrna_stat_in_sample.xls'),
            os.path.join(self.lncrna_stat.output_dir, 'lncrna_stat_in_category.xls')
        )

    @time_count
    def export_mirpre(self):
        self.mirpre_api = self.api.api("lnc_rna.mirna_precursor")

        self.mirpre_api.add_mirna_precursor(
            os.path.join(self.mir_pre.output_dir, 'miRNA_precursor.tabular')
        )

    @time_count
    def export_lncfamily(self):
        self.lncfamily_api = self.api.api("lnc_rna.lncrna_family")

        self.lncfamily_api.add_lncrna_family(
            os.path.join(self.lnc_family.output_dir, 'lncRNA_family.tabular')
        )


    @time_count
    def export_as(self):
        gevent.sleep()
        self.api_as = self.api.api("lnc_rna.rmats")
        gene_type_tsv = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_type.xls")
        for output_file in os.listdir(self.altersplicing.output_dir):
            test_group, ctrl_group = output_file.split('_vs_')
            # compare_plan = '{}|{}'.format(ctrl_group, test_group)
            compare_plan = '{}|{}'.format(test_group, ctrl_group)
            test_member = self.option('group_table').prop['group_dict'][test_group]
            ctrl_member = self.option('group_table').prop['group_dict'][ctrl_group]
            group_dict = {test_group: test_member, ctrl_group: ctrl_member}
            params = {
                'compare_plan': compare_plan,
                'control_id': str(self.control_id),
                'group_dict': group_dict,
                'group_id': str(self.group_id),
                'submit_location': 'splicingrmats',
                'task_id': self.task_id,
                'task_type': 2
            }
            self.logger.info(params)
            outpath = os.path.join(self.altersplicing.output_dir, output_file)
            self.api_as.add_sg_splicing_rmats(
                splicing_id=None,
                outpath=outpath,
                compare_plan=compare_plan,
                params=params,
                s3_output=self._sheet.output,
                gene_type_tsv=gene_type_tsv
            )

    @time_count
    def export_snp(self):
        gevent.sleep()
        self.logger.info("开始进行Snpfinal的导表")
        task_id = self.task_id
        project_sn = self.project_sn
        gene_type_tsv = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_type.xls")
        if self.option("snp_method").lower() == "gatk":
            params = dict(
                task_id=task_id,
                submit_location="snp",
                task_type=2,
                method_type="gatk"
            )
            self.api_snp = self.api.api("lnc_rna.snp_indel")
            snp_anno = self.snp.output_dir
            if os.path.exists(self.work_dir + "/snp_tmp"):
                shutil.rmtree(self.work_dir + "/snp_tmp")
            os.mkdir(self.work_dir + "/snp_tmp")

            self.api_snp.add_snp(module_output=snp_anno, params=params,
                                 method_type="gatk", upload_dir=self.work_dir + "/snp_tmp", gene_type_tsv=gene_type_tsv)
        if self.option("snp_method").lower() == "samtools":
            params = dict(
                task_id=task_id,
                submit_location="snp",
                task_type=2,
                method_type="samtools"
            )
            self.api_snp = self.api.api("lnc_rna.snp_indel")
            snp_anno = self.snp.output_dir
            if os.path.exists(self.work_dir + "/snp_tmp"):
                shutil.rmtree(self.work_dir + "/snp_tmp")
            os.mkdir(self.work_dir + "/snp_tmp")
            # os.link(os.path.join(snp_anno, 'snp_anno.xls'), os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))
            self.api_snp.add_snp(module_output=snp_anno, params=params,
                                 method_type="samtools",upload_dir=self.work_dir + "/snp_tmp", gene_type_tsv=gene_type_tsv)
        if os.path.exists(os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls')):
            os.remove(os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))
        if os.path.exists(os.path.join(snp_anno, 'snp_annotation.xls')):

            os.link(os.path.join(snp_anno, 'snp_annotation.xls'), os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))
        elif os.path.exists(os.path.join(snp_anno, 'data_anno_pre.xls')):
            os.link(os.path.join(snp_anno, 'data_anno_pre.xls'), os.path.join(self.work_dir, "snp_tmp", 'snp_anno.xls'))

    @time_count
    def export_expression(self):
        gevent.sleep()
        if os.path.exists(self.work_dir + "/gene.fpkm.matrix"):
            self.known_fpkm = self.work_dir + "/gene.fpkm.matrix"
        else:
            self.known_fpkm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/gene.tpm.matrix"):
            self.known_tpm = self.work_dir + "/gene.tpm.matrix"
        else:
            self.known_tpm = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_tpm_known.xls")


        if os.path.exists(self.work_dir + "/transcript.fpkm.matrix"):
            self.known_fpkm_t = self.work_dir + "/transcript.fpkm.matrix"
        else:
            self.known_fpkm_t = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_fpkm_known.xls")
        if os.path.exists(self.work_dir + "/transcript.tpm.matrix"):
            self.known_tpm_t = self.work_dir + "/transcript.tpm.matrix"
        else:
            self.known_tpm_t = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_tpm_known.xls")


        all_exp = self.api.api("lnc_rna.all_exp")
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
            if self.option('strand_dir')  in ["fr", 'f', "backward"]:
                if self.option("fq_type") == "PE":
                    self.libtype = "fr"
                else:
                    self.libtype = "f"
            else:
                if self.option("fq_type") == "PE":
                    self.libtype = "rf"
                else:
                    self.libtype = "r"
        else:
            self.libtype = None
        if self.option("exp_way").lower() == "tpm":
            if self.option("level").lower() == "transcript":
                exp_matrix = os.path.join(self.classify_quant.output_dir, 'transcript.tpm.matrix')
                ref_exp_matrix = self.known_tpm_t
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
                        rna_type = "all",
                        type='ref',
                    )
                    all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='T',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)

            exp_matrix = os.path.join(self.classify_quant.output_dir, 'gene.tpm.matrix')
            ref_exp_matrix = self.known_tpm
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
                    rna_type = "all",
                    type='ref',
                )
                all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='G',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)
        else:
            if self.option("level").lower() == "transcript":
                exp_matrix = os.path.join(self.classify_quant.output_dir, 'transcript.fpkm.matrix')
                ref_exp_matrix = self.known_tpm_t
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
                        rna_type = "all",
                        type='ref',
                    )
                    all_exp.add_distribution(ref_exp_matrix, group_dict, params=params_distribution, exp_level='T',
                                      quant_method=quant_method, project_sn=project_sn, task_id=task_id)

            exp_matrix = os.path.join(self.classify_quant.output_dir, 'gene.fpkm.matrix')
            ref_exp_matrix = self.known_tpm
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
                    rna_type="all",
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
                log_base="10",
                rna_type = "all",
                type="ref"
            )
            all_exp.add_exp_corr2(corr_output, exp_level='G', quant_method=quant_method, params=params,
                                  project_sn=project_sn, task_id=task_id)

            corr_output_t = self.exp_corr_t.work_dir
            params_t = params.copy()
            params_t['exp_level'] = "T"
            params_t['exp_id'] = str(trans_exp_id)
            all_exp.add_exp_corr2(corr_output_t, exp_level='T', quant_method=quant_method, params=params_t,
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
                    rna_type = "all",
                    type="ref"
                    # quant_method=quant_method,
                )
                pca_main_id = all_exp.add_exp_pca2(pca_output, quant_method=quant_method, exp_id=gene_exp_id, exp_level="G",
                                     params=params, project_sn=project_sn, task_id=task_id)
                if self.min_group_num >= 3:
                    all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(pca_main_id))

                pca_output_t = self.exp_pca_t.output_dir
                params = params.copy()
                params["exp_level"] = "T"
                params['exp_id'] = str(trans_exp_id)
                pca_main_id = all_exp.add_exp_pca2(pca_output_t, quant_method=quant_method, exp_id=gene_exp_id, exp_level="T",
                                     params=params, project_sn=project_sn, task_id=task_id)
                if self.min_group_num >= 3:
                    all_exp.insert_ellipse_table(self.ellipse_t.work_dir + '/ellipse_out.xls', str(pca_main_id))

            # add gene diff
            diff_output = self.diffexpress.output_dir
            uniform_output = self.diffexpress.uniform.output_dir
            gene_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_type.xls")
            exp_id, exp_level = gene_exp_id, 'G'
            diff_method = self.option('diff_method')
            stat_type = self.option('pvalue_padjust')
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                params = dict(
                    task_id=task_id,
                    submit_location="diff_detail",
                    task_type=2,
                    exp_id=str(exp_id),
                    group_id=str(group_id),
                    control_id=str(control_id),
                    exp_level=exp_level,
                    group_dict=group_dict,
                    # fc='{:g}'.format(self.option('fc')),
                    fc=str(float(self.option('fc'))),
                    # correct_method=self.option('padjust_way'),
                    stat_type=stat_type,
                    stat_cutoff=self.option('diff_fdr_ci'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    filter_method=self.option("filter_method"),
                    # quant_method=quant_method,
                    diff_method=diff_method,
                    rna_type = "all",
                    type='all',
                    correct_method=self.option('padjust_way'),
                    is_batch="False",
                )
            if self.option('diff_method').lower() in ['noiseq']:
                params = dict(
                    task_id=task_id,
                    submit_location="diff_detail",
                    task_type=2,
                    exp_id=str(exp_id),
                    group_id=str(group_id),
                    control_id=str(control_id),
                    exp_level=exp_level,
                    group_dict=group_dict,
                    # fc='{:g}'.format(self.option('fc')),
                    fc=str(float(self.option('fc'))),
                    # correct_method=self.option('padjust_way'),
                    prob=self.option('diff_fdr_ci'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    filter_method=self.option("filter_method"),
                    # quant_method=quant_method,
                    diff_method=diff_method,
                    rna_type = "all",
                    type='all',
                )
            all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                    group_dict=group_dict, group_id=group_id,
                                    exp_level=exp_level, quant_method=quant_method,
                                    diff_method=diff_method,
                                    project_sn=project_sn, task_id=task_id, params=params,
                                    pvalue_padjust=stat_type, seq_type=gene_type
                                    )
            # all_exp.add_diffexp(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
            #                     exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
            #                     project_sn=project_sn, task_id=task_id, params=params,
            #                     pvalue_padjust=stat_type, seq_type=gene_type)

            diff_output = self.diffexpress_t.output_dir
            uniform_output = self.diffexpress_t.uniform.output_dir
            gene_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_type.xls")
            exp_id, exp_level = trans_exp_id, 'T'
            diff_method = self.option('diff_method')
            stat_type = self.option('pvalue_padjust')

            params_t = params.copy()
            params_t.update({
                "exp_level": exp_level,
                "exp_id": str(trans_exp_id)
            })
            all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                    group_dict=group_dict, group_id=group_id,
                                    exp_level=exp_level, quant_method=quant_method,
                                    diff_method=diff_method,
                                    project_sn=project_sn, task_id=task_id, params=params_t,
                                    pvalue_padjust=stat_type, seq_type=gene_type
                                    )
            # all_exp.add_diffexp(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
            #                     exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
            #                     project_sn=project_sn, task_id=task_id, params=params_t,
            #                     pvalue_padjust=stat_type, seq_type=gene_type)

    @time_count
    def export_gene_detail(self):
        """
        导入基因详情表
        :return:
        """
        gevent.sleep()
        self.api_gene_detail = self.api.api('lnc_rna.add_gene_detail')
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
        gene_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/filtered_file/gene_type.xls")
        trans_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/filtered_file/trans_type.xls")
        lnc_ids = self.lnc_ids
        if self.option("is_assemble") == True:
            new_cds = self.annot_orfpfam.output_dir + "/novel_mrna.fa.transdecoder.cds"
            new_pep = self.annot_orfpfam.output_dir + "/novel_mrna.fa.transdecoder.pep"
        else:
            new_cds = None
            new_pep = None

        gene_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/gene_type.xls")
        trans_type = os.path.join(self.filter_by_exp.output_dir, "filtered_file/trans_type.xls")
        self.api_gene_detail.add_gene_detail(db_path, gene_bed, transcript_bed, species_urls, biomart_file, biomart_type,
                                             known_cds, known_pep, new_cds, new_pep, transcript_fasta, gene_fasta, t2g, gene_type, trans_type, lnc_ids)

    # 添加注释信息
    def merge_annotation_exp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        exp_output = os.path.join(self.output_dir, "classify_quant")
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
        diff_output = os.path.join(self.output_dir, "diffexpress")

        duplicate_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls') ## 为了防止流程重跑的时候反复增加注释结果表
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + '/' + '*_vs_*.*.xls')
        for each in target_files:
            if each.endswith("sizeFactor.xls"):
                continue
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_pd, gene_annot_pd], join='inner', axis=1)
            gene_out = each.split('.xls')[0] + '.annot.xls'
            gene_result.to_csv(gene_out, header=True, index=True, sep='\t')


        diff_output = os.path.join(self.output_dir, "diffexpress_t")
        duplicate_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls') ## 为了防止流程重跑的时候反复增加注释结果表
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + '/' + '*_vs_*.*.xls')
        all_annot = pd.read_table(annot, header=0, index_col=1)
        for each in target_files:
            if each.endswith("sizeFactor.xls"):
                continue
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_pd, all_annot], join='inner', axis=1)
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
