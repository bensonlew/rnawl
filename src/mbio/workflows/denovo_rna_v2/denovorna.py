# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""无参转录组工作流"""

import functools
import glob
import json
import os
import re
import shutil
import time
import unittest
from collections import OrderedDict,defaultdict

# from gevent.monkey import patch_all
import gevent
import pandas as pd
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow

from mbio.packages.denovo_rna_v2.copy_file import CopyFile
from mbio.packages.denovo_rna_v2.extract_annotation_result import RefAnnotation
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
from mbio.packages.dna_evolution.send_email import SendEmail
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


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


class DenovornaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(DenovornaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 单样本 OR 多样本
            # 样本个数 ['multiple', 'single']
            {"name": "sample_num", "type": "string", "default": "multiple"},  # 测序类型，single OR multiple

            ##基础参数设置
            # 数据类型 ['rawdata', 'cleandata']
            {'name': 'datatype', 'type': 'string', 'default': 'rawdata'},
            # 测序类型 ['PE', 'SE']
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            # 质控软件 ['fastp', 'seqprep']
            {'name': 'qc_soft', 'type': 'string', 'default': 'fastp'},
            # 测序质量 ['phred 33', 'phred 64']
            {'name': 'quality_score_system', 'type': 'string', 'default': 'phred 33'},
            # 链特异性 [False, True]
            {"name": "strand_specific", "type": "bool", "default": False},  # 是否有链特异性, 默认是False, 无特异性
            # 链特异性方向 ['PE': {'RF', 'FR'}, 'SE': {'R', 'F'}]
            {'name': 'strand_dir', 'type': 'string', 'default': 'forward'},  # 当链特异性时为True时，正义链为forward，反义链为reverse
            # 生物学重复 [True, False]
            {"name": "is_duplicate", "type": "bool", "default": ""},
            # 原始序列文件
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},
            # 质控序列文件
            {'name': 'qc_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # 分组方案
            {"name": "group", "type": "infile", "format": "denovo_rna_v2.group_table"},
            # 上机名称
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 对照组文件
            {"name": "control", "type": "infile", "format": "denovo_rna_v2.compare_table"},
            # 终止后续分析
            {'name': 'mapping_stop', 'type': 'bool', 'default': True},
            # 大于 %的样本
            {'name': 'mapping_sample_percent', 'type': 'float', 'default': 50.0},
            # Mapping Ratio 小于
            {'name': 'mapping_ratio', 'type': 'float', 'default': 50.0},
            # 终止后续分析
            {'name': 'rrna_stop', 'type': 'bool', 'default': True},
            # 大于 %的样本
            {'name': 'rrna_sample_percent', 'type': 'float', 'default': 50.0},
            # rRNA Ratio 小于 %
            {'name': 'rrna_ratio', 'type': 'float', 'default': 15.0},

            ##高级参数设置
            # 从头组装
            # 是否进行从头组装 [True, False]
            {"name": "assemble", "type": "bool", "default": True},
            # 组装方式[	'total','group','sample']
            {"name": "assemble_method", "type": "string", "default": 'total'},
            # 组装软件['Trinity','Spades']
            {"name": "assemble_soft", "type": "string", "default": 'Trinity'},

            # Trinity参数
            # 设置Trinity软件的kmer值
            {"name": "k_mer", "type": "int", "default": 31},
            # 设置Trinity软件的kmer值
            {"name": "min_kmer_cov", "type": "int", "default": 5},
            # 是否进行reads均一化 [True, False]
            {"name": "normalize", "type": "bool", "default": True},
            # reads均一化最大覆盖倍数
            {"name": "normalize_max_read_cov", "type": "int", "default": 50},
            # 设置最短contig的长度
            {"name": "min_contig_length", "type": "int", "default": 200},
            # 分割高密度基因区间基因，建议物种为Fungi时True，其他物种时False
            {"name": "jaccard_clip", "type": "bool", "default": False},
            # 不参与组装的样本
            {"name": "remove_samples", "type": "infile", 'format': "denovo_rna_v2.common"},
            # spades参数
            # 设置spades参数的kmer值
            {"name": "spades_k", "type": "string", "default": 'auto'},
            # 合并去冗余
            {"name": "merge_soft", "type": "string", "default": 'TGICL'},
            # 覆盖区域最小identity
            {"name": "tgicl_p", "type": "float", "default": 94},
            # 最小覆盖长度  [10-1000000]
            {"name": "tgicl_c", "type": "int", "default": 40},
            # 不选择从头组装
            # trinity组装结果文件
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"},
            # 基因和转录本对应关系文件
            {"name": "gene_to_trans", "type": "infile", 'format': "denovo_rna_v2.common"},

            # 分析水平 ['transcript','unigene']
            {"name": "level", "type": "string", 'default': "transcript"},

            # 是否进行组装结果优化
            {"name": "optimize", "type": "bool", "default": True},
            # 是否选择TransRate进行组装结果评估
            {"name": "transrate", "type": "bool", "default": True},
            # 利用cd-hit软件进行过滤时identity的取值，值为1时不进行过滤
            {"name": "cd_hit", "type": "float", "default": 0.99},
            # 过滤掉tpm小于给定值的转录本
            {"name": "tpm", "type": "float", "default": 0},
            {'type': 'int', 'name': 'filter_200', 'default': 0},
            {'type': 'int', 'name': 'filter_500', 'default': 0},
            {'type': 'int', 'name': 'filter_1000', 'default': 0},
            # 注释
            {"name": "nr_evalue", "type": "string", "default": "1e-5"},
            {"name": "string_evalue", "type": "string", "default": "1e-5"},
            {"name": "kegg_evalue", "type": "string", "default": "1e-5"},
            {"name": "swissprot_evalue", "type": "string", "default": "1e-5"},
            {"name": "pfam_evalue", "type": "string", "default": "1e-5"},
            {"name": "tf_evalue", "type": "string", "default": "1e-5"},

            # NR（一级分类）['Animal','Plant','Fungi','Protist','All']
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            # NR（二级分类）动物['Reptilia','Mammalia','Invertebrate','Fishes','Aves','Amphibia']植物['Algae','Spermatophyta','OtherPlant']
            {"name": "nr_sub_database", "type": "string", "default": None},  # nr库类型
            # kegg库类型['Animal','Plant','Fungi','Protist','All']
            {"name": "kegg_database", "type": "string", "default": "All"},  # kegg注释库类型
            # tf_database注释库类型 ['Plant','Animal','Other']
            {"name": "tf_database", "type": "string", "default": ""},

            # z注释
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "tf_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg,swissprot,pfam'},
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_2019"},

            # 表达及差异表达
            # 选择表达量计算软件，RSEM,Kallisto,Salmon，定量指标为TPM
            {"name": "express_method", "type": "string", "default": "RSEM"},
            # 表达定量指标 ['tpm', 'fpkm']
            {'name': 'exp_way', 'type': 'string', 'default': 'tpm'},
            # 差异分析软件 DESeq2,DEGseq,edgeR
            {"name": "diff_method", "type": "string", "default": "DESeq2"},
            # 选择判断显著性水平的指标 ['padjust','pvalue']
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},
            # 显著性水平 [0-1]
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},
            # 多重检验校正方法 [ 'Bonferroni','Holm','BH','BY']
            {"name": "padjust_way", "type": "string", "default": "BH"},
            # Fold Change，设置差异倍数阈值
            {"name": "fc", "type": "string", "default": 2},

            # # SNP/InDel分析
            # # 流程是否进行SNP分析 [True,False]
            # {"name": "is_snp", "type": "bool", "default": True},

            # SNP分析 ['True', 'False', 'Skip']
            {'name': 'is_snp', 'type': 'string', 'default': 'True'},

            # call_snp方法
            {"name": "snp_method", "type": "string", "default": 'Sentieon'},
            # {"name": "qual", "type": "float", "default": 20},  #过滤低质量的SNP位点
            # {"name": "dp", "type": "int", "default": 1},  #过滤低质量的SNP位点
            # {"name": "fastp", "type": "int", "default": 0},  #是否选用fastp进行质控
            {'name': 'report_img', 'type': 'bool', 'default': False},
            {'name': 'get_run_log', 'type': 'bool', 'default': False},
        ]

        # 获取输出目录
        self.workflow_output_tmp = self._sheet.output
        '''
        ## 暂时将上传目录设置为对象存储目录， 如果框架修改、需要进行相应的改动
        region_bucket = Config().get_project_region_bucket(project_type="denovo_rna_v2")
        if re.match(r"^\w+://\S+/.+$", workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        elif workflow_output_tmp.startswith("/mnt/ilustre/tsanger-data/"):
            self.workflow_output = self.workflow_output_tmp.replace(region_bucket, '/mnt/ilustre/tsanger-data/')
        '''
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong", code="12000414")

        self.project_sn = self._sheet.project_sn  # 获取project_sn
        self.task_id = self._sheet.id  # 获取task_id
        self.add_option(options)
        self.set_options(self._sheet.options())

        # 添加tool/module
        self.filecheck = self.add_tool("denovo_rna_v2.filecheck_denovo")
        if self.option('qc_soft') == "fastp":
            self.qc = self.add_module("denovo_rna_v2.fastp_rna")
        else:
            self.qc = self.add_module("denovo_rna_v2.hiseq_qc")
        self.hiseq_reads_stat_raw = self.add_module("denovo_rna_v2.hiseq_reads_stat")
        self.hiseq_reads_stat_clean = self.add_module("denovo_rna_v2.hiseq_reads_stat")
        self.assemble = self.add_module("denovo_rna_v2.denovo_assemble3")
        self.assemble_filter = self.add_module("denovo_rna_v2.denovo_assemble2_filter")
        self.filter_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")
        self.unigene_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")
        self.filter_unigene_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")
        # self.diamond = self.add_module("denovo_rna_v2.diamond")
        # self.blast = self.add_module("denovo_rna_v2.blast")
        self.annotation = self.add_module("denovo_rna_v2.denovo_annotation")
        if self.option("express_method") == "RSEM":
            self.align = self.add_module("denovo_rna_v2.quant")
        else:
            self.align = self.add_module("denovo_rna_v2.quant")
            self.express = self.add_module("denovo_rna_v2.quant")
        if self.option('sample_num') == 'multiple':
            self.exp_pca = self.add_tool("denovo_rna_v2.exp_pca")
            self.exp_corr = self.add_tool("denovo_rna_v2.exp_corr")
            group_spname = self.option("group").get_group_spname()
            group_dict = self.option("group").prop['group_dict']
            group_snum = [len(group_spname[g]) for g in group_spname]
            min_group_num = min(group_snum)
            if min_group_num >= 3:
                # self.ellipse = self.add_tool('graph.ellipse')
                self.ellipse = self.add_tool('denovo_rna_v3.ellipse')
            if len(group_dict) > 1:
                self.exp_venn = self.add_tool('ref_rna_v2.exp_venn')
            # self.diffexpress_gene = self.add_tool("denovo_rna_v2.diffexp")
            # self.diffexpress_gene = self.add_tool("denovo_rna_v3.batch.diffexp_batch")
            self.diffexpress_gene = self.add_module("denovo_rna_v3.diffexp_batch_new")
            if self.option("level").lower() == "transcript":
                # self.diffexpress_trans = self.add_tool("denovo_rna_v2.diffexp")
                # self.diffexpress_trans = self.add_tool("denovo_rna_v3.batch.diffexp_batch")
                self.diffexpress_trans = self.add_module("denovo_rna_v3.diffexp_batch_new")
        self.ssr = self.add_tool("denovo_rna_v2.ssr")
        self.cds_predict = self.add_module("denovo_rna_v2.annot_orfpfam")
        if self.option("snp_method").lower() == "samtools":
            self.snp = self.add_module("denovo_rna_v3.snp")
        elif self.option("snp_method").lower() == "sentieon":
            self.snp = self.add_module("denovo_rna_v3.sentieon")
        else:
            self.snp = self.add_module("denovo_rna_v3.call_snp_indel_gatk")
        self.snpfinal = self.add_tool("denovo_rna_v3.snpfinal_new2")
        self.diamond = self.add_module("denovo_rna_v2.annot_mapdb")
        self.nr_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.swissprot_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.cog_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.kegg_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.pfam_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.annot_stat = self.add_module("denovo_rna_v2.annot_class_beta")
        self.sequence_detail = self.add_tool("denovo_rna_v3.sequence_detail")
        self.diff_geneset_analysis = self.add_module("denovo_rna_v2.workflow_diffgt.diff_geneset_all_pipline")
        # 判断流程结束tool/module list
        if self.option("sample_num") == "multiple":
            self.final_tools = [self.diffexpress_gene,
                                self.ssr,
                                self.annot_stat,
                                self.cds_predict,
                                self.exp_corr,
                                self.diff_geneset_analysis]
            if self.option("is_snp") == "True":
                self.final_tools.append(self.snp)
                self.final_tools.append(self.snpfinal)
            if self.option("level").lower() == "transcript":
                self.final_tools.append(self.diffexpress_trans)
            group_spname = self.option("group").get_group_spname()
            group_dict = self.option("group").prop['group_dict']
            group_snum = [len(group_spname[g]) for g in group_spname]
            min_group_num = min(group_snum)
            if min_group_num >= 3:
                self.final_tools.append(self.ellipse)
            if len(group_dict) > 1:
                self.final_tools.append(self.exp_venn)
            if self.option("group").prop["sample_number"] > 2:
                self.final_tools.append(self.exp_pca)
        else:
            self.final_tools = [self.ssr, self.annot_stat, self.cds_predict, self.align]
        if self.option("express_method").lower() != "rsem":
            self.final_tools.append(self.express)
        if self.option("optimize") == True:
            self.final_tools.append(self.filter_evaluation)
        if self.option("level").lower() == "transcript":
            self.final_tools.append(self.unigene_evaluation)
            self.final_tools.append(self.filter_unigene_evaluation)
        # 添加step，显示在页面进度条
        if self.option('datatype') == 'rawdata':
            self.step.add_steps('filecheck', 'qualitycontrol', 'hiseq_reads_stat_raw', 'hiseq_reads_stat_clean',
                                "annotation", "express", "ssr_analysis", "cds_predict")
            # if not self.option('assemble') == False:
            self.step.add_steps('assemble')
            self.final_tools.append(self.hiseq_reads_stat_raw)
            self.final_tools.append(self.hiseq_reads_stat_clean)
        else:
            self.step.add_steps('filecheck', 'qualitycontrol', 'hiseq_reads_stat_clean',
                                "annotation", "express", "ssr_analysis", "cds_predict")
            # if not self.option('assemble') == False:
            self.step.add_steps('assemble')
        if self.option("sample_num") == "multiple":
            self.step.add_steps('diff_geneset_analysis')
            if self.option("optimize") == True and self.option("is_snp") == "True":
                self.step.add_steps("evaluation", "diff_gene", 'unigene_evaluation', "snp_analysis")
                if self.option("level").lower() == "transcript":
                    self.step.add_steps("diff_trans", "unigene_evaluation", "filter_unigene_evaluation")
            elif self.option("optimize") == True and self.option("is_snp") == "False":
                self.step.add_steps("evaluation", "diff_gene", 'unigene_evaluation')
                if self.option("level").lower() == "transcript":
                    self.step.add_steps("diff_trans", "unigene_evaluation", "filter_unigene_evaluation")
            elif self.option("optimize") == False and self.option("is_snp") == "True":
                self.step.add_steps("diff_gene", "snp_analysis")
                if self.option("level").lower() == "transcript":
                    self.step.add_steps("diff_trans", "unigene_evaluation", "filter_unigene_evaluation")
            else:
                if self.option("level").lower() == "transcript":
                    self.step.add_steps("diff_gene", "diff_trans", "unigene_evaluation", "filter_unigene_evaluation")
        else:
            # if self.option("optimize") == True:
            self.step.add_steps("assemble", "evaluation", "unigene_evaluation", "filter_unigene_evaluation")

        if self.option('sample_num') == 'multiple':
            group_size = list()
            group_dict = self.option("group").prop['group_dict']
            for key in group_dict:
                group_size.append(len(group_dict[key]))
            group_size.sort()
            if group_size[0] == group_size[-1] == 1:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "DEGseq")
                    self.option("diff_fdr_ci", '0.001')
                    self.logger.info("该项目没有生物学重复,不可以使用DESeq2,修改方法为DEGseq,阈值设置为0.001")
            elif group_size[0] == 1 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "edgeR")
                    self.option("diff_fdr_ci", '0.05')
                    self.logger.info("该项目部分组别有生物学重复,部分组别无生物学重复,不可以使用DESeq2,修改方法为edgeR,阈值设置为0.05")
            elif group_size[0] >= 2 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "degseq":
                    self.option("diff_method", "DESeq2")
                    self.option("diff_fdr_ci", '0.05')
                    self.logger.info("该项目有生物学重复,不可以使用DEGseq,修改方法为DESeq2,阈值设置为0.05")
            else:
                pass
        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        data = os.path.join(self.work_dir, 'data.json')
        if self._sheet.rerun:
            self.delete_mongo_data()
        # if os.path.exists(data):
        #     with open(data, 'r') as load_f:
        #         load_dict = json.load(load_f)
        #         if 'rerun' in load_dict and load_dict['rerun']:
        #             self.logger.info("该项目重运行中，先删除mongo库中已有数据")
        #             self.delete_mongo_data()
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))
        self.logger.info("{}".format(self.annot_config_dict))

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self.task_id, 'denovo_rna_v2')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'denovo_rna_v2')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))

    def check_options(self):
        """
        检查选项
        """
        # 基础参数
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE", code="12000401")

        if self.option('sample_num') == 'single':
            self.option('is_duplicate', False)

        # 组装相关参数
        if self.option("assemble") == False:
            if not self.option("assembly_file").is_set:
                raise OptionError("不进行拼接的时候，必须输入组装结果文件", code="12000402")
            if not self.option("gene_to_trans").is_set:
                raise OptionError("不进行拼接的时候，必须输入基因和转录本对应关系文件", code="12000403")
        if not self.option("k_mer") >= 20 and not self.option("string_evalue") <= 32:
            raise OptionError("kmer超出范围，值的范围为[25-32]", code="12000404")
        if self.option("assemble") == True:
            if not self.option("assemble_soft") in ['Trinity', 'Spades']:
                raise OptionError("组装软件应为Trinity或Spades", code="12000413")

        # 注释相关参数
        try:
            nr_evalue = float(self.option("nr_evalue"))
            string_evalue = float(self.option("string_evalue"))
            kegg_evalue = float(self.option("string_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
            tf_evalue = float(self.option("tf_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范", code="12000405")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("string_blast_evalue", string_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
            self.option("tf_blast_evalue", tf_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围", code="12000406")
        if not self.option("string_blast_evalue") > 0 and not self.option("string_blast_evalue") < 1:
            raise OptionError("String比对的E值超出范围", code="12000407")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围", code="12000408")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围", code="12000409")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围", code="12000410")
        if not self.option("tf_blast_evalue") > 0 and not self.option("tf_blast_evalue") < 1:
            raise OptionError("tf分析的E值超出范围", code="12000411")
        if self.option("tf_database") == "":
            raise OptionError("必须选择转录因子预测的数据库", code="12000412")
        if self.option("kegg_database") == "All":
            self.option("kegg_database", "All")
        elif self.option("kegg_database") == "Animal":
            self.option("kegg_database", "Animals")
        elif self.option("kegg_database") == "Plant":
            self.option("kegg_database", "Plants")
        elif self.option("kegg_database") == "Protist":
            self.option("kegg_database", "Protists")
        if self.option("nr_database") == "All":
            self.option("nr_database", "nr")
        else:
            nr = self.option("nr_database").lower()
            self.option("nr_database", nr)
        if self.option('rrna_stop'):
            if not isinstance(self.option('rrna_sample_percent'), float):
                raise OptionError('核糖体判断中止条件的样本比例值应为浮点数', code='12000414')
            if not isinstance(self.option('rrna_ratio'), float):
                raise OptionError('核糖体判断中止条件的阈值比例值应为浮点数', code='12000415')
        if self.option('mapping_stop'):
            if not isinstance(self.option('mapping_sample_percent'), float):
                raise OptionError('比对结果判断中止条件的样本比例值应为浮点数', code='12000416')
            if not isinstance(self.option('mapping_ratio'), float):
                raise OptionError('比对结果判断中止条件的阈值比例值应为浮点数', code='12000417')
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_filecheck(self):
        self.logger.info("开始运行文件检查")
        opts = {
            # 'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'group_table': self.option('group'),
            'control_file': self.option('control'),
            'assembly_file': self.option('assembly_file'),
            'gene_to_trans': self.option('gene_to_trans'),
            'sample_num': self.option('sample_num'),
            'is_duplicate': self.option('is_duplicate')
        }
        if self.option('datatype') == 'rawdata':
            opts.update({'fastq_dir': self.option('fastq_dir')})
        else:
            opts.update({'fastq_dir': self.option('qc_dir')})
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        self.logger.info("开始运行质控")
        if self.option('qc_soft') == "fastp":
            with open(self.option('fastq_dir').prop['path'] + '/list.txt', 'r') as list_r, open(
                    self.option('fastq_dir').prop['path'] + '/tmp.txt', 'w') as tmp_w:
                for line in list_r:
                    tmp_w.write(self.option('fastq_dir').prop['path'] + '/' + line)
            self.qc.set_options({
                'sample_path': self.option('fastq_dir').prop['path'] + '/tmp.txt',
                'length_required': '30',
                'fq_type': self.option('fq_type')
            })
        else:
            self.qc.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type')
            })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.qualitycontrol})
        self.qc.on('end', self.set_step, {'end': self.step.qualitycontrol})
        self.qc.run()

    def run_hiseq_reads_stat_raw(self):
        self.logger.info("开始运行质控前统计")
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        self.hiseq_reads_stat_raw.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'quality': quality
        })
        self.hiseq_reads_stat_raw.on('start', self.set_step, {'start': self.step.hiseq_reads_stat_raw})
        self.hiseq_reads_stat_raw.on('end', self.set_step, {'end': self.step.hiseq_reads_stat_raw})
        self.hiseq_reads_stat_raw.on('end', self.set_output, 'hiseq_reads_stat_raw')
        self.hiseq_reads_stat_raw.run()

    def run_hiseq_reads_stat_clean(self):
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                fastq_dir = self.qc.option('sickle_dir')
            elif self.option('qc_soft') == 'seqprep':
                fastq_dir = self.qc.option('sickle_dir')
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir')
        self.hiseq_reads_stat_clean.set_options({
            'fastq_dir': fastq_dir,
            'fq_type': self.option('fq_type'),
            'quality': quality,
            'dup': True,
            'rfam': True
        })
        self.hiseq_reads_stat_clean.on('start', self.set_step, {'start': self.step.hiseq_reads_stat_clean})
        self.hiseq_reads_stat_clean.on('end', self.set_step, {'end': self.step.hiseq_reads_stat_clean})
        self.hiseq_reads_stat_clean.on('end', self.set_output, 'hiseq_reads_stat_clean')
        # if self.option('datatype') == 'cleandata':
        # self.hiseq_reads_stat_clean.on('end', self.check_rrna)
        self.hiseq_reads_stat_clean.run()

    def check_rrna(self):
        if self.option('rrna_stop'):
            df = pd.read_table(os.path.join(self.hiseq_reads_stat_clean.output_dir, 'stat_results'))
            rrna_sample_percent = len(
                [i for i in df.iloc[:, -1] if i > self.option('rrna_ratio')]
            ) / float(df.shape[0]) * 100
            if rrna_sample_percent > self.option('rrna_sample_percent'):
                self.stop('rrna')

    def run_assemble(self):
        self.logger.info("开始运行组装")
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "sample_fq_list": fq_list,
            "fq_type": self.option("fq_type"),
            "strand_direct": self.option("strand_dir"),
            "group": self.option("group"),
            "assemble_soft": self.option("assemble_soft").lower(),
            "assemble_method": self.option("assemble_method")
        }
        if self.option("assemble_soft").lower() == "trinity":
            opts.update({
                "filter_sample": self.option("remove_samples"),
                "min_contig_length": self.option("min_contig_length"),
                "kmer_size": self.option("k_mer"),
                "min_kmer_cov": self.option("min_kmer_cov"),
                "jaccard_clip": self.option("jaccard_clip"),
                "no_normalize_reads": not self.option("normalize"),
                "normalize_max_read_cov": self.option("normalize_max_read_cov"),
            })
        else:
            opts.update({
                "spades_k": self.option("spades_k"),
            })

        if self.option("assemble_method") != "total":
            opts.update({
                "tgicl_p": self.option("tgicl_p"),
                "tgicl_c": self.option("tgicl_c"),
            })
        self.assemble.set_options(opts)
        self.assemble.on("end", self.set_output, "assemble")
        self.assemble.on('start', self.set_step, {'start': self.step.assemble})
        self.assemble.on('end', self.set_step, {'end': self.step.assemble})
        self.assemble.run()

    def run_assemble_filter(self):
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
            #add by fwy 20201106用户上传clean data很可能是gz格式，不可以直接使用
            fq_list_new = os.path.join(self.option("qc_dir").prop['path'], "fq_list_new.txt")
            with open(fq_list_new,"w") as n,open(fq_list,"r") as r:
                for line in r.readlines():
                    line = line.strip().split()
                    line_new = [info.split(".gz")[0] for info in line]
                    n.write("\t".join(line_new) + "\n")
            os.rename(os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt"),os.path.join(self.option("qc_dir").prop['path'], "fq_list_old.txt"))
            os.rename(os.path.join(self.option("qc_dir").prop['path'], "fq_list_new.txt"),
                      os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt"))
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        if self.option("assemble") == True:
            opts = {
                # "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
                "sample_fq_list": fq_list,
                "fq_type": self.option("fq_type"),
                "assemble_fa": self.assemble.output_dir + '/assemble_raw.fasta',
                "assemble_g2t": self.assemble.output_dir + '/assemble_raw.gene_trans_map',
                "min_contig_length": self.option("min_contig_length"),
                "species": self.option("kegg_database"),
            }
        else:
            opts = {
                "sample_fq_list": fq_list,
                "assemble_fa": self.option("assembly_file").prop['path'],
                "assemble_g2t": self.option("gene_to_trans").prop['path'],
                "species": self.option("kegg_database"),
            }
        if self.option("optimize") == True:
            opts.update({
                "transrate_filter": self.option("transrate"),
                "cdhit_filter": True,
                "TPM_filter": True,
                "cdhit_identity_threshold": self.option("cd_hit"),
                "TPM_threshold": self.option("tpm"),
            })
        else:
            opts.update({
                "transrate_filter": False,
                "cdhit_filter": False,
                "TPM_filter": False,
            })
        self.assemble_filter.set_options(opts)
        self.assemble_filter.on("end", self.set_output, "assemble")
        self.assemble_filter.on('start', self.set_step, {'start': self.step.assemble})
        self.assemble_filter.on('end', self.set_step, {'end': self.step.assemble})
        self.assemble_filter.run()

    def run_filter_evaluation(self):
        self.rawt2g = glob.glob(self.assemble_filter.output_dir + "/*filter.gene_trans_map")[0]
        self.logger.info("开始运行组装优化后评估")
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "sample_fq_list": fq_list,
            "fq_type": self.option("fq_type"),
            "species": self.option("kegg_database"),
            "min_contig_length": self.option("min_contig_length"),
            "assemble_fa": self.assemble_filter.option("filter_fa"),
            "assemble_g2t": self.rawt2g
            # "assemble_g2t" : os.path.join(self.assemble_filter.output_dir,"Trinity.filter.gene_trans_map"),
        }
        self.filter_evaluation.set_options(opts)
        self.filter_evaluation.on("end", self.set_output, "filter_evaluation")
        self.filter_evaluation.on('start', self.set_step, {'start': self.step.evaluation})
        self.filter_evaluation.on('end', self.set_step, {'end': self.step.evaluation})
        self.filter_evaluation.run()

    def run_unigene_evaluation(self):
        self.logger.info("开始运行组装优化后评估")
        seq = glob.glob(self.assemble_filter.output_dir + '/trinity_stat/unigene.fasta')[0]
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "sample_fq_list": fq_list,
            "fq_type": self.option("fq_type"),
            "min_contig_length": self.option("min_contig_length"),
            "species": self.option("kegg_database"),
            "assemble_fa": seq,
        }
        self.unigene_evaluation.set_options(opts)
        self.unigene_evaluation.on("end", self.set_output, "unigene_evaluation")
        self.unigene_evaluation.on('start', self.set_step, {'start': self.step.unigene_evaluation})
        self.unigene_evaluation.on('end', self.set_step, {'end': self.step.unigene_evaluation})
        self.unigene_evaluation.run()

    def run_filter_unigene_evaluation(self):
        self.logger.info("开始运行组装优化后评估")
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        seq = glob.glob(self.assemble_filter.output_dir + '/Trinity.filter.unigene.fasta')[0]
        opts = {
            # "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "sample_fq_list": fq_list,
            "fq_type": self.option("fq_type"),
            "min_contig_length": self.option("min_contig_length"),
            "species": self.option("kegg_database"),
            "assemble_fa": seq,
        }
        self.filter_unigene_evaluation.set_options(opts)
        self.filter_unigene_evaluation.on("end", self.set_output, "filter_unigene_evaluation")
        self.filter_unigene_evaluation.on('start', self.set_step, {'start': self.step.filter_unigene_evaluation})
        self.filter_unigene_evaluation.on('end', self.set_step, {'end': self.step.filter_unigene_evaluation})
        self.filter_unigene_evaluation.run()

    def run_align(self):
        self.logger.info("开始运行样本比对和表达定量分析")
        self.rawt2g = glob.glob(self.assemble_filter.output_dir + "/*filter.gene_trans_map")[0]
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
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "fastq" : self.qc.option("fq_list"),
            "fastq": fq_list,
            "method": "RSEM",
            "libtype": self.libtype,
            "transcriptome": self.assemble_filter.option("filter_fa"),
            "t2g": self.rawt2g
            # "t2g" : os.path.join(self.assemble.output_dir,"Trinity.filter.gene_trans_map"),
        }
        self.align.set_options(opts)
        self.align.on("end", self.set_output, "align")
        self.align.on('start', self.set_step, {'start': self.step.express})
        self.align.on('end', self.set_step, {'end': self.step.express})
        self.align.on('end', self.check_mapping)
        self.align.run()

    def run_express(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option('strand_dir') == 'forward':
                if self.option("fq_type") == "PE":
                    self.libtype = "rf"
                else:
                    self.libtype = "r"
            else:
                if self.option("fq_type") == "SE":
                    self.libtype = "fr"
                else:
                    self.libtype = "f"
        else:
            self.libtype = None
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "fastq": self.qc.option("fq_list"),
            "fastq": fq_list,
            "method": self.option("express_method"),
            "libtype": self.libtype,
            "transcriptome": self.assemble_filter.option("filter_fa"),
            "t2g": self.rawt2g
            # "t2g" : os.path.join(self.assemble.output_dir,"Trinity.filter.gene_trans_map"),
        }
        self.express.set_options(opts)
        self.express.on("end", self.set_output, "express")
        # self.express.on('start', self.set_step, {'start': self.step.express})
        # self.express.on('end', self.set_step, {'end': self.step.express})
        self.express.run()

    def check_mapping(self):
        if self.option('mapping_stop'):
            with open(os.path.join(self.align.output_dir, "alignment_rate.txt"), "r") as a:
                n = 0.0
                for i, stat in enumerate(a.readlines()[1:]):
                    if float(stat.strip().split()[-1]) * 100 < float(self.option('mapping_ratio')):
                        n += 1
                else:
                    if n / (i + 1) * 100 > float(self.option('mapping_sample_percent')):
                        self.stop('mapping')

    def run_exp_pca(self):
        self.logger.info("开始运行pca")
        if self.option("express_method") == "RSEM":
            if self.option("exp_way").lower() == "fpkm":
                opts = {
                    "exp": self.align.option("gene_fpkm")
                }
            else:
                opts = {
                    "exp": self.align.option("gene_tpm")
                }
        else:
            if self.option("exp_way").lower() == "fpkm":
                opts = {
                    "exp": self.express.option("gene_fpkm")
                }
            else:
                opts = {
                    "exp": self.express.option("gene_tpm")
                }

        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        # self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        # self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        if hasattr(self, 'ellipse'):
            self.exp_pca.on('end', self.run_ellipse)
        self.exp_pca.run()

    def run_ellipse(self):
        self.ellipse.set_options({
            'analysis': 'pca',
            'group_table': self.option('group').prop['path'],
            'pc_table': os.path.join(self.exp_pca.output_dir, 'PCA.xls'),
        })
        # self.ellipse.on('start', self.set_step, {'start': self.step.exp_pca})
        # self.ellipse.on('end', self.set_step, {'end': self.step.exp_pca})
        self.ellipse.on('end', self.set_output, 'exp_pca')
        self.ellipse.run()

    def run_exp_corr(self):
        self.logger.info("开始运行聚类分析")
        if self.option("express_method") == "RSEM":
            if self.option("exp_way").lower() == "fpkm":
                opts = {
                    "exp": self.align.option("gene_fpkm")
                }
            else:
                opts = {
                    "exp": self.align.option("gene_tpm")
                }
        else:
            if self.option("exp_way").lower() == "fpkm":
                opts = {
                    "exp": self.express.option("gene_fpkm")
                }
            else:
                opts = {
                    "exp": self.express.option("gene_tpm")
                }

        self.exp_corr.set_options(opts)
        self.exp_corr.on("end", self.set_output, "exp_corr")
        # self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        # self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.run()

    def run_exp_venn(self):
        if self.option("express_method") == "RSEM":
            if self.option("exp_way").lower() == "fpkm":
                express_matrix = self.align.option('gene_fpkm').path
            else:
                express_matrix = self.align.option('gene_tpm').path
        else:
            if self.option('exp_way').lower() == 'fpkm':
                express_matrix = self.express.option('gene_fpkm').path
            else:
                express_matrix = self.express.option('gene_tpm').path
        group_dict = self.option('group').prop['group_dict']
        if len(group_dict) <= 6:
            group_table = self.option('group').path
        else:
            group_table = '{}.6'.format(self.option('group').path)
            lines = ['#sample\tgroup\n']
            for i, (g, ss) in enumerate(group_dict.items()):
                if i < 6:
                    lines.extend('{}\t{}\n'.format(s, g) for s in ss)
            else:
                open(group_table, 'w').writelines(lines)
        self.exp_venn.set_options({
            'express_matrix': express_matrix,
            'group_table': group_table,
            'threshold': 1
        })
        # self.exp_venn.on('start', self.set_step, {'start': self.step.exp_corr})
        # self.exp_venn.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_venn.on('end', self.set_output, 'exp_corr')
        self.exp_venn.run()

    def run_diffexpress_gene(self):
        self.logger.info("开始运行基因差异表达分析")
        if self.option("express_method") == "RSEM":
            opts = {
                "count": self.align.option("gene_count").path,
                "group": self.option("group").path,
                "cmp": self.option("control"),
                "fc": float(self.option("fc")),
                "method": self.option("diff_method")
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

            if self.option("exp_way").lower() == "fpkm":
                opts.update({"exp_type": "fpkm"})
                opts.update({"exp_matrix": self.align.option("gene_fpkm")})
            else:
                opts.update({"exp_matrix": self.align.option("gene_tpm")})
        else:
            opts = {
                "count": self.express.option("gene_count"),
                "exp_matrix": self.express.option("gene_tpm"),
                "group": self.option("group"),
                "cmp": self.option("control"),
                "fc": float(self.option("fc")),
                "method": self.option("diff_method")
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
        self.diffexpress_gene.set_options(opts)
        self.diffexpress_gene.on("end", self.set_output, "diffexpress_gene")
        self.diffexpress_gene.on('start', self.set_step, {'start': self.step.diff_gene})
        self.diffexpress_gene.on('end', self.set_step, {'end': self.step.diff_gene})
        self.diffexpress_gene.run()

    def run_diffexpress_trans(self):
        self.logger.info("开始运行转录本差异表达分析")
        if self.option("express_method") == "RSEM":
            opts = {
                "count": self.align.option("transcript_count"),
                "group": self.option("group"),
                "cmp": self.option("control"),
                "fc": float(self.option("fc")),
                "method": self.option("diff_method"),
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
            if self.option("exp_way").lower() == "fpkm":
                opts.update({"exp_type": "fpkm"})
                opts.update({"exp_matrix": self.align.option("transcript_fpkm")})
            else:
                opts.update({"exp_matrix": self.align.option("transcript_tpm")})
        else:
            opts = {
                "count": self.express.option("transcript_count"),
                "exp_matrix": self.express.option("transcript_tpm"),
                "group": self.option("group"),
                "cmp": self.option("control"),
                "fc": float(self.option("fc")),
                "method": self.option("diff_method"),
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
        self.diffexpress_trans.set_options(opts)
        self.diffexpress_trans.on("end", self.set_output, "diffexpress_trans")
        self.diffexpress_trans.on('start', self.set_step, {'start': self.step.diff_trans})
        self.diffexpress_trans.on('end', self.set_step, {'end': self.step.diff_trans})
        self.diffexpress_trans.run()

    def run_cds_predict(self):
        self.logger.info("开始运行cds预测")
        t2g2u = glob.glob(self.assemble_filter.output_dir + "/*filter.gene_trans_map")[0]
        nrxml = os.path.join(self.diamond.output_dir, "nr", "blast.xml")
        swissxml = os.path.join(self.diamond.output_dir, "swissprot", "blast.xml")
        opts = {
            "fasta": self.assemble_filter.option("filter_fa").prop['path'],
            "species_type": self.option("tf_database"),
            "pfam_version": self.annot_config_dict['pfam']['version'],
            "blast_nr_xml": nrxml,
            "blast_swissprot_xml": swissxml,
            "isoform_unigene": t2g2u,
            # "isoform_unigene" : os.path.join(self.assemble.output_dir,"Trinity.filter_t2g2u"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.set_output, "cds_predict")
        self.cds_predict.on("end", self.run_pfam_filter)
        self.cds_predict.on('start', self.set_step, {'start': self.step.cds_predict})
        self.cds_predict.on('end', self.set_step, {'end': self.step.cds_predict})
        self.cds_predict.run()

    def run_diamond(self):
        self.logger.info("开始运行blast注释")
        if self.option("annot_group") in ["REFRNA_GROUP_202110"]:
            diamond_version = "v2.0.13"
        else:
            diamond_version = "v0.9.24.125"
        blast_opts = {
            'query': self.assemble_filter.option("filter_fa").prop['path'],
            'method': 'diamond',
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
            "diamond_version" : diamond_version
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            # self.diamond = self.add_module("denovo_rna_v2.annot_mapdb")
            if self.option("nr_sub_database") in ["Spermatophyta", "Reptilia", "OtherPlant", "Mammalia", "Invertebrate",
                                                  "Fishes", "Aves", "Amphibia", "Algae"]:
                blast_opts.update(
                    {
                        'nr_db': self.option("nr_sub_database"),
                    }
                )
            else:
                blast_opts.update(
                    {
                        'nr_db': self.option("nr_database"),
                    }
                )
        self.diamond.set_options(blast_opts)
        self.diamond.on('end', self.set_output, 'diamond')
        self.diamond.on('end', self.run_nr_blast_filter)
        if 'cog' in self.option('database'):
            self.diamond.on('end', self.run_cog_blast_filter)
        if 'kegg' in self.option('database'):
            self.diamond.on('end', self.run_kegg_blast_filter)
        if 'swissprot' in self.option('database'):
            self.diamond.on('end', self.run_swissprot_blast_filter)
        self.diamond.run()

    def run_nr_blast_filter(self):
        # self.nrxml=glob.glob(os.path.join(self.diamond.output_dir,"xml")+"/"+"*nr.xml")[0]
        self.nrxml = os.path.join(self.diamond.output_dir, "nr", "blast.xml")
        options = {
            'xml': self.nrxml,
            'types': "xml",
            'evalue': self.option('nr_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.nr_filter.set_options(options)
        self.nr_filter.run()

    def run_swissprot_blast_filter(self):
        # self.swissxml = glob.glob(os.path.join(self.diamond.output_dir, "xml") + "/" + "*swissprot.xml")[0]
        self.swissxml = os.path.join(self.diamond.output_dir, "swissprot", "blast.xml")
        options = {
            'xml': self.swissxml,
            'types': "xml",
            'evalue': self.option('swissprot_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.swissprot_filter.set_options(options)
        self.swissprot_filter.run()

    def run_cog_blast_filter(self):
        # self.cogxml = glob.glob(os.path.join(self.diamond.output_dir, "xml") + "/" + "*eggnog.xml")[0]
        self.cogxml = os.path.join(self.diamond.output_dir, "eggnog", "blast.xml")
        options = {
            'xml': self.cogxml,
            'types': "xml",
            'evalue': self.option('string_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.cog_filter.set_options(options)
        self.cog_filter.run()

    def run_kegg_blast_filter(self):
        # self.keggxml = glob.glob(os.path.join(self.diamond.output_dir, "xml") + "/" + "*kegg.xml")[0]
        self.keggxml = os.path.join(self.diamond.output_dir, "kegg", "blast.xml")
        options = {
            'xml': self.keggxml,
            'types': "xml",
            'evalue': self.option('kegg_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.kegg_filter.set_options(options)
        self.kegg_filter.run()

    def run_pfam_filter(self):
        pfamdomain = glob.glob(self.cds_predict.output_dir + "/*pfam_domain")[0]
        options = {
            'hmm': pfamdomain,
            # 'hmm': self.cds_predict.output_dir + "/pfam_domain",
            'types': "hmm",
            'evalue': self.option('pfam_evalue'),
        }
        self.pfam_filter.set_options(options)
        self.pfam_filter.run()

    def run_annotation(self):
        self.logger.info("开始运行注释统计")
        self.newt2g = glob.glob(self.cds_predict.output_dir + "/*all_tran2gen.txt")[0]
        self.go_table = glob.glob(self.diamond.output_dir + "/GO/*blast2go_merge.xls")[0]
        anno_opts = {
            'db': 'nr,swissprot,kegg,eggnog,pfam,go',
            'type': 'new',
            'gene2trans': self.newt2g,
            # "gene2trans" : os.path.join(self.assemble.output_dir,"Trinity.filter.gene_trans_map"),
            "taxonomy": self.option("kegg_database"),
            "blast_nr_xml": self.nr_filter.option('outxml').prop['path'],
            "blast_kegg_xml": self.kegg_filter.option('outxml').prop['path'],
            "blast_eggnog_xml": self.cog_filter.option('outxml').prop['path'],
            "blast_swissprot_xml": self.swissprot_filter.option('outxml').prop['path'],
            "pfam_domain": self.pfam_filter.option('outtable').prop['path'],
            'blast2go_annot': self.go_table,
            'link_bgcolor': 'green',
            'png_bgcolor': '#00CD00',
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "ncbi_taxonmy_version": self.annot_config_dict['ncbi_taxonomy']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
        }
        self.annot_stat.set_options(anno_opts)
        self.annot_stat.on('end', self.set_output, 'annotation')
        self.annot_stat.on('start', self.set_step, {'start': self.step.annotation})
        self.annot_stat.on('end', self.set_step, {'end': self.step.annotation})
        self.annot_stat.run()

    def run_snp(self):
        # cdsbed = glob.glob(self.cds_predict.output_dir + "/" + "*all_predicted.bed")[0]
        # newt2g = glob.glob(self.cds_predict.output_dir + "/*all_tran2gen.txt")[0]
        if self.option('datatype') == 'rawdata':
            fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        else:
            fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            "ref_fasta": self.assemble_filter.option("unigene_filter_fa").prop["path"],
            # "fq_list" : self.qc.option("fq_list").prop["path"],
            "fq_list": fq_list,
            "fq_type": self.option("fq_type"),
            # "cds_bed":cdsbed,
            # "allt2g" :newt2g,
            'call_type': self.option("snp_method").lower(),
            # "allt2g":os.path.join(self.assemble.output_dir,"Trinity.filter.gene_trans_map"),
            # "anno":os.path.join(self.annot_stat.output_dir,"all_annot.xls"),
        }
        self.snp.set_options(opts)
        self.snp.on('end', self.set_output, 'snp')
        self.snp.on('start', self.set_step, {'start': self.step.snp_analysis})
        self.snp.on('end', self.set_step, {'end': self.step.snp_analysis})
        self.snp.run()

    def run_snpfinal(self):
        cdsbed = glob.glob(self.cds_predict.output_dir + "/" + "*all_predicted.bed")[0]
        newt2g = glob.glob(self.cds_predict.output_dir + "/*all_tran2gen.txt")[0]
        # if self.option('datatype') == 'rawdata':
        #     fq_list = os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt")
        # else:
        #     fq_list = os.path.join(self.option("qc_dir").prop['path'], "fq_list.txt")
        opts = {
            # "ref_fasta": self.assemble_filter.option("unigene_filter_fa").prop["path"],
            # "fq_list" : self.qc.option("fq_list").prop["path"],
            # "fq_list": fq_list,
            # "fq_type": self.option("fq_type"),
            "cds_bed": cdsbed,
            "allt2g": newt2g,
            'method': self.option("snp_method").lower(),
            "anno": os.path.join(self.annot_stat.output_dir, "all_annot.xls"),
        }
        if self.option("snp_method").lower() == "samtools":
            with open(self.snp.work_dir + "/bamlist") as len_bam:
                sample_num = len(len_bam.readlines())
            opts.update(
                {
                    'bamlist': sample_num,
                    'call_vcf': os.path.join(self.snp.work_dir, "VcfFilterSamtools/output/final.vcf")
                }
            )
        else:
            with open(os.path.join(self.snp.work_dir, "bam_list")) as len_bam:
                sample_num = len(len_bam.readlines())
            if os.path.exists(os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")):
                filted_snp_vcf = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")
            else:
                filted_snp_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.snp.filter.recode.vcf")
            if os.path.exists(os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")):
                filted_indel_vcf = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
            else:
                filted_indel_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.indel.filter.recode.vcf")
            opts.update(
                {
                    'bamlist': sample_num,
                    'filted_snp_vcf': filted_snp_vcf,
                    'filted_indel_vcf': filted_indel_vcf
                    # 'filted_snp_vcf': os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf"),
                    # 'filted_indel_vcf': os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
                }
            )

        self.snpfinal.set_options(opts)
        self.snpfinal.on('end', self.set_output, 'snp')
        # self.snp.on('start', self.set_step, {'start': self.step.snp_analysis})
        # self.snp.on('end', self.set_step, {'end': self.step.snp_analysis})
        self.snpfinal.run()

    def run_ssr(self):
        self.logger.info("开始运行SSR分析")
        self.cdsbed = glob.glob(self.cds_predict.output_dir + "/" + "*all_predicted.bed")[0]
        opts = {
            "unigene_fa": self.assemble_filter.option("unigene_filter_fa").prop['path'],
            "bed": self.cdsbed,
        }
        self.ssr.set_options(opts)
        self.ssr.on('end', self.set_output, 'ssr')
        self.ssr.on('start', self.set_step, {'start': self.step.ssr_analysis})
        self.ssr.on('end', self.set_step, {'end': self.step.ssr_analysis})
        self.ssr.run()

    def run_diff_geneset_analysis(self):
        if self.option('exp_way').lower() == 'fpkm':
            exp = self.align.option('gene_fpkm').path
        else:
            exp = self.align.option('gene_tpm').path
        opts = {
            'diff_path': self.diffexpress_gene.output_dir,
            'annot_result': self.annot_stat.output_dir,
            'gene_count_file': self.align.option('gene_count').path,
            'diff_method': self.option('diff_method'),
            'gene_exp_file': exp,
            'group': self.option('group').path,
            'level': "G",
            "kegg_version":self.annot_config_dict['kegg']['version']
        }
        self.diff_geneset_analysis.set_options(opts)
        self.diff_geneset_analysis.on('start', self.set_step, {'start': self.step.diff_geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_step, {'end': self.step.diff_geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_output, 'diff_geneset_analysis')
        self.diff_geneset_analysis.run()

    def run_seq_detail(self):
        cds = os.path.join(self.cds_predict.output_dir, 'all_predicted.cds.fa')
        pep = os.path.join(self.cds_predict.output_dir, 'all_predicted.pep.fa')
        fasta = os.path.join(self.assemble_filter.output_dir, 'Trinity.filter.fasta')
        trans2unigene = os.path.join(self.cds_predict.output_dir, 'all_tran2gen.txt')
        opts = {
            'cds_seq': cds,
            'pep_seq': pep,
            'txpt_seq': fasta,
            'trans2gene': trans2unigene,
        }
        self.sequence_detail.set_options(opts)
        self.sequence_detail.on('end', self.set_db)
        self.sequence_detail.run()

    def stop(self, reason):
        assert reason in ['rrna', 'mapping']
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.denovo_task_info')
        task_info.add_task_info()
        # self.export_genome_info()
        api_qc = self.api.api("denovo_rna_v2.denovo_rna_qc")
        fq_type = self.option("fq_type").lower()
        if self.option('datatype') == 'rawdata':
            api_qc = self.api.api("denovo_rna_v2.denovo_rna_qc")
            qc_stat = self.hiseq_reads_stat_raw.output_dir
            fq_type = self.option("fq_type").lower()
            api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before")
            # self.api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="before")
            quality_stat_after = self.hiseq_reads_stat_clean.output_dir + "/qualityStat"  # 质控数据结果统计
            quality_stat_before = self.hiseq_reads_stat_raw.output_dir + "/qualityStat"  # 原始数据结果统计
            api_qc.add_gragh_info(quality_stat_before, "before")
        quality_stat_after = self.hiseq_reads_stat_clean.output_dir + "/qualityStat"  # 质控数据结果统计
        qc_stat = self.hiseq_reads_stat_clean.output_dir
        api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="after")
        # self.api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="after")
        api_qc.add_gragh_info(quality_stat_after, "after")
        if self.option("group").is_set:
            self.group_id, self.group_detail, self.group_category = api_qc.add_specimen_group(
                self.option("group").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
            if self.option('productive_table').is_set:
                api_qc.add_productive_name(samples=self.option('group').prop['sample'],
                                   productive_table=self.option('productive_table').path)
        if self.option("control").is_set:
            self.control_id, compare_detail = api_qc.add_control_group(self.option("control").prop["path"],
                                                                       self.group_id)
            self.compare_detail = compare_detail
        if reason == 'rrna':
            msg = 'Workflow will stop because the rRNA ratio is not up to par'
            self.logger.warn(msg)
        elif reason == 'mapping':
            self.export_denovo_align()
            msg = 'Workflow will stop because the alignment rate is not up to par'
            self.logger.warn(msg)
        receiver = ['caiping.shi@majorbio.com', 'rna_bioinfor@majorbio.com']
        a = SendEmail("897236887@qq.com", "smtp.qq.com", "fhwuvcclstjqbfga", "897236887@qq.com", ','.join(receiver),
                      "WARNING - project_sn ({}), task_id ({})".format(self._sheet.project_sn, self._sheet.id), 465)
        a.send_msg(msg)
        a.send_email()
        # db = Config().get_mongo_client(mtype='project', dydb_forbid=True)[Config().get_mongo_dbname(mtype='project', dydb_forbid=True)]
        # email = db['sg_task_email']
        # email.update({'task_id': self.task_id}, {'$set': {'status': '5'}}, upsert=True)
        super(DenovornaWorkflow, self).end()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="12000413")
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
        if event['data'] == 'hiseq_reads_stat_raw':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
            self.logger.info("开始设置qc的输出目录")
        if event['data'] == 'hiseq_reads_stat_clean':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'annotation')
        if event['data'] == 'diamond':
            self.move2outputdir(obj.output_dir, 'diamond')
        if event['data'] == 'assemble':
            self.move2outputdir(obj.output_dir, 'assemble')
        if event['data'] == 'evaluation':
            self.move2outputdir(obj.output_dir, 'assemble_evaluation')
        if event['data'] == 'filter_evaluation':
            self.move2outputdir(obj.output_dir, 'filter_evaluation')
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'snp')
        if event['data'] == 'align':
            self.move2outputdir(obj.output_dir, 'align')
        if event['data'] == 'express':
            self.move2outputdir(obj.output_dir, 'express')
        if event['data'] == 'cds_predict':
            self.move2outputdir(obj.output_dir, 'cds_predict')
        if event['data'] == 'exp_pca':
            self.move2outputdir(obj.output_dir, 'exp_pca')
        if event['data'] == 'exp_corr':
            self.move2outputdir(obj.output_dir, 'exp_corr')
        if event['data'] == 'diffexpress':
            self.move2outputdir(obj.output_dir, 'diffexpress')
        if event['data'] == 'diff_geneset_analysis':
            self.move2outputdir(obj.output_dir, 'diff_geneset_analysis')

    def run(self):
        """
        denovo-rna workflow run方法
        """
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("sample_num") != "multiple":
            if self.option('datatype') != 'rawdata':
                fq_path = self.option('qc_dir').path
                self.logger.info("处理处理qc_dir{}".format(fq_path))
                with open(os.path.join(fq_path, 'list.txt'), "r") as fl, open(os.path.join(self.work_dir, "group"),
                                                                              "w") as g:
                    items = fl.readline().strip().split("\t")
                    sample_name = items[1]
                    g.write("#sample\tgroup\n")
                    g.write(sample_name + "\t" + sample_name + "\n")
                self.set_options(dict(group=os.path.join(self.work_dir, "group")))
            else:
                fq_path = self.option('fastq_dir').path
                self.logger.info("处理处理qc_dir{}".format(fq_path))
                with open(os.path.join(fq_path, 'list.txt'), "r") as fl, open(os.path.join(self.work_dir, "group"),
                                                                              "w") as g:
                    items = fl.readline().strip().split("\t")
                    sample_name = items[1]
                    g.write("#sample\tgroup\n")
                    g.write(sample_name + "\t" + sample_name + "\n")
                self.set_options(dict(group=os.path.join(self.work_dir, "group")))
        if self.option('datatype') == 'rawdata':
            self.on_rely([self.hiseq_reads_stat_raw, self.hiseq_reads_stat_clean],
                         self.check_rrna)  # 当选择原始数据时，质控后数据结果统计进行rRNA检查
            self.filecheck.on('end', self.run_hiseq_reads_stat_raw)  # 质控前统计
            self.filecheck.on('end', self.run_qc)  # 质控
            self.qc.on('end', self.run_hiseq_reads_stat_clean)  # 质控后统计
        else:
            dct = dict()
            fq_path = self.option('qc_dir').path
            for line in open(os.path.join(fq_path, 'list.txt')):
                items = line.strip().split('\t')
                if items[1] in dct:
                    dct[items[1]][items[2]] = os.path.join(fq_path, items[0])
                    flag = 1
                else:
                    dct[items[1]] = {items[2]: os.path.join(fq_path, items[0])}
                    flag = 0
            else:
                self.fq_list = os.path.join(fq_path, 'fq_list.txt')
                open(self.fq_list, 'w').writelines(['{}\t{l}\t{r}\n'.format(k, **v) for k, v in dct.items()] if flag \
                                                       else ['{}\t{s}\n'.format(k, **v) for k, v in dct.items()])
            self.filecheck.on('end', self.run_hiseq_reads_stat_clean)  # 当客户选择上传质控数据时，直接进行质控后统计
        if self.option("assemble") == True:
            if self.option("rrna_stop"):
                self.hiseq_reads_stat_clean.on('end', self.run_assemble)  # 当客户选择rRNA截断时，质控统计结束后进行组装
            else:
                if self.option('datatype') == 'rawdata':# 当客户选择不需要rRNA截断时,如果选择raw_data,质控后直接组装
                    self.qc.on('end', self.run_assemble)
                else:# 当客户选择不需要rRNA截断时,如果选择clean_data,质控后直接组装
                    self.hiseq_reads_stat_clean.on('end', self.run_assemble)
                # self.qc.on('end', self.run_assemble)  # 当客户选择不需要rRNA截断时，质控统计结束后进行组装
            self.assemble.on('end', self.run_assemble_filter)  # 无论客户选择是否优化，均进行这一步，区别在于参数不同
        else:
            if self.option("rrna_stop"):
                self.hiseq_reads_stat_clean.on('end', self.run_assemble_filter)  # 当客户选择rRNA截断时，质控统计结束后进行组装
            else:
                if self.option('datatype') == 'rawdata':
                    self.qc.on('end', self.run_assemble_filter)  # 当客户选择不需要rRNA截断时，质控统计结束后进行组装
                else:
                    self.hiseq_reads_stat_clean.on('end', self.run_assemble_filter)
        if self.option("optimize") == True:
            self.assemble_filter.on('end', self.run_filter_evaluation)  # 当客户选择组装优化时进行优化结果评估
        if self.option("level").lower() == "transcript":
            self.assemble_filter.on('end', self.run_filter_unigene_evaluation)  # 当客户选择分析水平为转录本时，进行以下两步
            self.assemble_filter.on('end', self.run_unigene_evaluation)
        if self.option("express_method") == "RSEM":  #
            self.assemble_filter.on('end', self.run_align)
        else:
            self.assemble_filter.on('end', self.run_align)
            self.assemble_filter.on('end', self.run_express)
        self.assemble_filter.on('end', self.run_diamond)
        self.diamond.on('end', self.run_cds_predict)
        annot_filter_tools = []
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            annot_filter_tools.append(self.nr_filter)
        if 'cog' in self.option('database'):
            annot_filter_tools.append(self.cog_filter)
        if 'kegg' in self.option('database'):
            annot_filter_tools.append(self.kegg_filter)
        if 'swissprot' in self.option('database'):
            annot_filter_tools.append(self.swissprot_filter)
        if 'pfam' in self.option('database'):
            annot_filter_tools.append(self.pfam_filter)
        self.on_rely(annot_filter_tools, self.run_annotation)
        self.cds_predict.on('end', self.run_ssr)
        if self.option("sample_num") == "multiple":
            self.on_rely([self.annot_stat, self.diffexpress_gene], self.run_diff_geneset_analysis)
            if self.option("express_method") == "RSEM":
                if self.option("level").lower() == "transcript":
                    self.align.on('end', self.run_diffexpress_gene)
                    self.align.on('end', self.run_diffexpress_trans)
                else:
                    self.align.on('end', self.run_diffexpress_gene)
                if self.option("group").prop["sample_number"] > 2:
                    self.align.on('end', self.run_exp_pca)
                if len(self.option('group').prop['group_dict']) > 1:
                    self.align.on('end', self.run_exp_venn)
                self.align.on('end', self.run_exp_corr)
                if self.option("is_snp") == "True":
                    self.on_rely([self.hiseq_reads_stat_clean, self.assemble_filter], self.run_snp)
                    self.on_rely([self.annot_stat, self.cds_predict, self.assemble_filter, self.snp], self.run_snpfinal)
            else:
                if self.option("level").lower() == "transcript":
                    self.express.on('end', self.run_diffexpress_gene)
                    self.express.on('end', self.run_diffexpress_trans)
                else:
                    self.express.on('end', self.run_diffexpress_gene)
                if self.option("group").prop["sample_number"] > 2:
                    self.express.on('end', self.run_exp_pca)
                if len(self.option('group').prop['group_dict']) > 1:
                    self.express.on('end', self.run_exp_venn)
                self.express.on('end', self.run_exp_corr)
                if self.option("is_snp") == "True":
                    self.on_rely([self.hiseq_reads_stat_clean, self.assemble_filter], self.run_snp)
                    self.on_rely([self.annot_stat, self.cds_predict, self.assemble_filter, self.snp], self.run_snpfinal)
        self.on_rely(self.final_tools, self.run_chart)
        # self.on_rely(self.final_tools, self.set_db)
        self.run_filecheck()
        super(DenovornaWorkflow, self).run()

    def run_chart(self):
        '''
        绘图步骤插入在导表前
        '''
        self.chart = self.add_tool('denovo_rna_v2.chart')
        samples = self.option("group").prop["sample"]
        chart_dict = {
            "type": "workflow",
            "samples": samples,
            "qc_file_use": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=self.hiseq_reads_stat_clean.output_dir, sample_name='{sample_name}'),
        }

        if self.option('datatype') == 'rawdata':
            chart_dict.update({
                "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                    table_dir=self.hiseq_reads_stat_raw.output_dir, sample_name='{sample_name}')
            })
        if self.option('fq_type') == "SE":
            chart_dict.update({
                "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.hiseq_reads_stat_raw.output_dir, sample_name='{sample_name}'),
                "qc_file_use": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.hiseq_reads_stat_clean.output_dir, sample_name='{sample_name}')
            })
        if self.option("group").is_set:
            group_dict = self.option("group").prop["group_dict"]
            chart_dict.update({
                "group_dict": group_dict
            })
        if self.option("control").is_set:
            cmp_list = self.option("control").prop["cmp_list"]
            chart_dict.update({
                "cmp_list": cmp_list
            })
        assemble_filter_dir = self.assemble_filter.output_dir
        filter_evolution = self.filter_evaluation.output_dir

        if self.option("level").lower() == "transcript":
            unigene_evaluation = self.unigene_evaluation.output_dir
            filter_unigene_evaluation = self.filter_unigene_evaluation.output_dir
        else:
            unigene_evaluation = assemble_filter_dir
            filter_unigene_evaluation = filter_evolution
        if self.option("optimize") == "True" or self.option("optimize") == True:
            chart_dict.update({
                "filter_assemble_length": "{table_dir}/trinity_stat/trans_count_stat_500.txt".format(
                    table_dir=filter_evolution),
                "filter_assemble_length_gene": "{table_dir}/trinity_stat/unigene_count_stat_500.txt".format(
                    table_dir=filter_unigene_evaluation),
            })
        else:
            chart_dict.update({
                "assemble_length": "{table_dir}/trinity_stat/trans_count_stat_500.txt".format(
                    table_dir=assemble_filter_dir),
                "assemble_length_gene": "{table_dir}/trinity_stat/unigene_count_stat_500.txt".format(
                    table_dir=unigene_evaluation),
            })


        #annotation
        chart_dict.update({
            "annot_stat": "{table_dir}/all_stat.xls".format(
                table_dir=self.annot_stat.output_dir)
        })
        # annotation nr pie
        chart_dict.update(
            {"annot_nr_species_pie_gene": os.path.join(self.annot_stat.output_dir, 'nr', 'gene_nr_species_stat.xls'),
             "annot_nr_evalue_pie_gene": os.path.join(self.annot_stat.output_dir, 'nr', 'gene_nr_evalue.xls'),
             "annot_nr_similar_pie_gene": os.path.join(self.annot_stat.output_dir, 'nr', "gene_nr_similar.xls")
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update(
                {"annot_nr_species_pie_trans": os.path.join(self.annot_stat.output_dir, 'nr',
                                                           'tran_nr_species_stat.xls'),
                 "annot_nr_evalue_pie_trans": os.path.join(self.annot_stat.output_dir, 'nr', 'trans_nr_evalue.xls'),
                 "annot_nr_similar_pie_trans": os.path.join(self.annot_stat.output_dir, 'nr', "trans_nr_similar.xls")
                 })
        # annotation swissprot pie
        chart_dict.update(
            {"annot_swissprot_evalue_pie_gene": os.path.join(self.annot_stat.output_dir, 'swissprot', 'gene_swissprot_evalue.xls'),
             "annot_swissprot_similar_pie_gene": os.path.join(self.annot_stat.output_dir, 'swissprot', "gene_swissprot_similar.xls")
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update(
                {"annot_swissprot_evalue_pie_trans": os.path.join(self.annot_stat.output_dir, 'swissprot', 'trans_swissprot_evalue.xls'),
                 "annot_swissprot_similar_pie_trans": os.path.join(self.annot_stat.output_dir, 'swissprot', "trans_swissprot_similar.xls")
                 })
        # annotation pfam bar
        chart_dict.update({
            "annot_pfam_bar_gene": os.path.join(self.annot_stat.output_dir, 'pfam', 'pfam_domain_gene.xls')
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update({
                "annot_pfam_bar_trans": os.path.join(self.annot_stat.output_dir, 'pfam', 'pfam_domain_tran.xls')
            })
        # annotation cog bar
        chart_dict.update({
            "annot_cog_bar_gene": os.path.join(self.annot_stat.output_dir, 'cog', 'summary.G.tsv')
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update({
                "annot_cog_bar_trans": os.path.join(self.annot_stat.output_dir, 'cog', 'summary.T.tsv')
            })
        # annotation go pie level
        chart_dict.update({
            "annot_go2_pie_gene": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev2_gene.stat.xls'),
            "annot_go3_pie_gene": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev3_gene.stat.xls'),
            "annot_go4_pie_gene": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev4_gene.stat.xls'),
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update({
            "annot_go2_pie_trans": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev2_tran.stat.xls'),
            "annot_go3_pie_trans": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev3_tran.stat.xls'),
            "annot_go4_pie_trans": os.path.join(self.annot_stat.output_dir, 'go', 'go_lev4_tran.stat.xls'),
            })
        # annotation kegg layer
        chart_dict.update({
            "annot_kegg_layer_gene": os.path.join(self.annot_stat.output_dir, 'kegg', 'kegg_layer_gene.xls')
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update({
                "annot_kegg_layer_trans": os.path.join(self.annot_stat.output_dir, 'kegg', 'kegg_layer_tran.xls')
            })
        # tf
        if self.option("tf_database").lower() == "animal":
            tf_unigene_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_unigene_animal.stat')
            tf_transcript_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_transcript_animal.stat')
        elif self.option("tf_database").lower() == "plant":
            tf_unigene_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_unigene_plant.stat')
            tf_transcript_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_transcript_plant.stat')
        else:
            tf_unigene_path = ""
            tf_transcript_path = ""
        if tf_unigene_path:
            chart_dict.update({
                "tf_stat_gene": tf_unigene_path
            })
        if self.option("level").lower() == "transcript" and tf_transcript_path:
            chart_dict.update({
            "tf_stat_trans": tf_transcript_path
            })

        #quant
        if self.option("express_method") == "RSEM":
            exp_output = self.align.output_dir
        else:
            exp_output = self.express.output_dir
        gene_exp = os.path.join(exp_output, 'gene.{}.matrix'.format(self.option('exp_way').lower()))
        trans_exp = os.path.join(exp_output, 'transcript.{}.matrix'.format(self.option('exp_way').lower()))
        chart_dict.update({
            "gene_exp_all": gene_exp,
            "tran_exp_all": trans_exp
        })
        # exp venn pca corr
        if self.option('sample_num') == 'multiple':
            if self.option('group').prop['sample_number'] > 2:
                chart_dict.update({
                    "exp_pca_file": "{table_dir}/PCA.xls".format(
                        table_dir=self.exp_pca.output_dir),
                    "exp_pca_var_file": "{table_dir}/Explained_variance_ratio.xls".format(
                        table_dir=self.exp_pca.output_dir),
                    "exp_corr_file": "{table_dir}/sample_correlation.xls".format(
                        table_dir=self.exp_corr.work_dir),
                    "exp_corr_tree_file": "{table_dir}/sample.cluster_tree.txt".format(
                        table_dir=self.exp_corr.work_dir)
                })
            if len(self.option("group").prop['group_dict']) > 1:
                chart_dict.update({
                    "exp_venn": self.exp_venn.output_dir + "/venn_graph.xls"
                })
        if hasattr(self, 'ellipse'):
            chart_dict.update({
                "exp_pca_ellipse": "{table_dir}/ellipse_out.xls".format(table_dir=self.ellipse.output_dir)
            })

        #CDS
        chart_dict.update({
            'gene_cds_summary': os.path.join(self.cds_predict.output_dir, "cds_len_unigene.txt")
        })
        if self.option("level").lower() == "transcript":
            chart_dict.update({
                'trans_cds_summary': os.path.join(self.cds_predict.output_dir, 'cds_len_transcript.txt')
            })
        #SSR
        chart_dict.update({
            'ssr_stat': os.path.join(self.ssr.work_dir, 'ssr_type.txt')
        })
        # Diff
        if self.option('sample_num') == 'multiple':
            chart_dict.update({
                'gene_diff_summary': os.path.join(self.diffexpress_gene.uniform.output_dir, 'json'),
                'gene_diff_volcano': os.path.join(self.diffexpress_gene.uniform.output_dir, 'all_volcano.txt'),
                'gene_diff_scatter': os.path.join(self.diffexpress_gene.uniform.output_dir, 'all_scatter.txt')
            })
            if self.option("level").lower() == "transcript":
                chart_dict.update({
                    'trans_diff_summary': os.path.join(self.diffexpress_trans.uniform.output_dir, 'json'),
                    'trans_diff_volcano': os.path.join(self.diffexpress_trans.uniform.output_dir, 'all_volcano.txt'),
                    'trans_diff_scatter': os.path.join(self.diffexpress_trans.uniform.output_dir, 'all_scatter.txt')
                })
        # geneset
        if self.option('sample_num') == 'multiple':
            genesets = sorted([os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.diff_geneset_analysis.output_dir))])
            chart_dict.update({"genesets":genesets})
            chart_dict.update({"diff_geneset_venn":self.diff_geneset_analysis.output_dir})
            if len(genesets) > 0:
                cluster_geneset_name = os.listdir(os.path.join(self.diff_geneset_analysis.output_dir,"cluster"))[0]
                cluster_dir = os.path.join(self.diff_geneset_analysis.output_dir,"cluster",cluster_geneset_name)
                chart_dict.update({
                    "cluster_geneset_name":cluster_geneset_name,
                    "cluster_exp" : os.path.join(cluster_dir,"expression_matrix.xls"),
                    "cluster_tree" : os.path.join(cluster_dir,"seq.cluster_tree.txt"),
                    "sample_tree" : os.path.join(cluster_dir, "sample.cluster_tree.txt"),
                    "subcluster_list" : glob.glob(cluster_dir + "/*subcluster_*.xls"),
                    "gene_annot_file" :  os.path.join(self.annot_stat.output_dir, "all_annot.xls")
                })
                file_json_path = os.path.join(self.diff_geneset_analysis.file_prepare.output_dir, "prepare_json")
                with open(file_json_path, "r") as j:
                    file_dict = json.load(j)
                kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
                chart_dict.update({
                    "cog_class":"{table_dir}/{geneset_name}/diff_cog_class/cog_class_table.xls".format(table_dir = self.diff_geneset_analysis.output_dir,geneset_name = "{geneset_name}"),
                    "go_class":"{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(table_dir = self.diff_geneset_analysis.output_dir,geneset_name = "{geneset_name}"),
                    "go_enrich": "{table_dir}/{geneset_name}/diff_go_enrich/go_enrich_geneset_list_gene.xls".format(
                        table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                    "kegg_class": "{table_dir}/{geneset_name}/diff_kegg_class/kegg_stat.xls".format(
                            table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                    "kegg_enrich": "{table_dir}/{geneset_name}/diff_kegg_enrich/enrich/{geneset_name}_gene.list.DE.list.check.kegg_enrichment.xls".format(
                                table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                    "kegg_level":kegg_level_path

                })
        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
            json.dump(chart_dict, json_f, sort_keys=True, indent=4)
        self.chart.set_options({
            "file_json": self.work_dir + "/chart_workflow.json"
        })
        self.chart.on('end', self.run_seq_detail)
        self.chart.run()

    def get_options(self):
        options = {}
        for name in self._options:
            if "file" in self._options[name]._type :
                try:
                    options[name] = self.option(name).prop["path"]
                except:
                    pass
            else:
                options[name] = self._options[name].value
        return options

    def get_function_infos(self):
        function_infos = defaultdict(dict)
        #build_seq_database
        function_infos["build_seq_database"]["cds"] = os.path.join(self.cds_predict.output_dir, 'all_predicted.cds.fa')
        function_infos["build_seq_database"]["pep"] = os.path.join(self.cds_predict.output_dir, 'all_predicted.pep.fa')
        function_infos["build_seq_database"]["fasta"] = os.path.join(self.assemble_filter.output_dir, 'Trinity.filter.fasta')
        function_infos["build_seq_database"]["trans2unigene"] =  os.path.join(self.cds_predict.output_dir, 'all_tran2gen.txt')
        function_infos["build_seq_database"]["seq_db"] = os.path.join(self.work_dir, 'seq_db.sqlite3')
        function_infos["build_seq_database"]["gene_stat"] = os.path.join(self.sequence_detail.work_dir, 'gene_stat')
        function_infos["build_seq_database"]["trans_stat"] = os.path.join(self.sequence_detail.work_dir, 'tran_stat')

        #export_qc
        function_infos["export_qc"]["qc_stat_before"] = self.hiseq_reads_stat_raw.output_dir
        function_infos["export_qc"]["quality_stat_after"] = self.hiseq_reads_stat_clean.output_dir + "/qualityStat"  # 质控数据结果统计
        function_infos["export_qc"]["quality_stat_before"] = self.hiseq_reads_stat_raw.output_dir + "/qualityStat"  # 原始数据结果统计
        function_infos["export_qc"]["qc_stat_after"] = self.hiseq_reads_stat_clean.output_dir

        #export_denovo_assembly
        function_infos["export_denovo_assembly"]["assemble_filter_dir"] = self.assemble_filter.output_dir
        function_infos["export_denovo_assembly"]["filter_evolution"] = self.filter_evaluation.output_dir
        if self.option("level").lower() == "transcript":
            function_infos["export_denovo_assembly"]["unigene_evaluation"] = self.unigene_evaluation.output_dir
            function_infos["export_denovo_assembly"]["filter_unigene_evaluation"] = self.filter_unigene_evaluation.output_dir
        else:
            function_infos["export_denovo_assembly"]["unigene_evaluation"] = self.assemble_filter.output_dir
            function_infos["export_denovo_assembly"]["filter_unigene_evaluation"] = self.filter_evaluation.output_dir

        #export_denovo_align
        function_infos["export_denovo_align"]["result_dir"] = self.align.output_dir

        #export_denovo_annotation
        function_infos["export_denovo_annotation"]["result_dir"] = self.annot_stat.output_dir
        function_infos["export_denovo_annotation"]["trans2gene"] = os.path.join(self.annot_stat.output_dir, "all_tran2gene.txt")
        if self.option("express_method") == "RSEM":
            exp_output = self.align.output_dir
            function_infos["export_denovo_annotation"]["exp_output"] = self.align.output_dir
        else:
            exp_output = self.express.output_dir
            function_infos["export_denovo_annotation"]["exp_output"] = self.express.output_dir
        function_infos["export_denovo_annotation"]["gene_exp"] = os.path.join(exp_output, 'gene.tpm.matrix')
        function_infos["export_denovo_annotation"]["trans_exp"] = os.path.join(exp_output, 'transcript.tpm.matrix')

        #export_snp
        function_infos["export_snp"]["snpfinal_work_dir"] = self.snpfinal.work_dir

        #export_ssr
        function_infos["export_ssr"]["ssr_work_dir"] = self.ssr.work_dir

        #export_tf
        function_infos["export_tf"]["cds_predict_work_dir"] = self.cds_predict.work_dir

        #export_cds
        function_infos["export_cds"]["cds_predict_output_dir"] = self.cds_predict.output_dir

        #export_expression
        if self.option("strand_specific") == True:
            if self.option('strand_dir') == 'forward':
                if self.option("fq_type") == "PE":
                    libtype = "rf"
                else:
                    libtype = "r"
            else:
                if self.option("fq_type") == "SE":
                    libtype = "fr"
                else:
                    libtype = "f"
        else:
            libtype = None
        function_infos["export_expression"]["libtype"] = libtype
        if self.option("sample_num") == "multiple":
            function_infos["export_expression"]["corr_output"] = self.exp_corr.work_dir
            function_infos["export_expression"]["graph_table"] = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
            function_infos["export_expression"]["pca_output"] = self.exp_pca.output_dir
            if hasattr(self, 'ellipse'):
                function_infos["export_expression"]["ellipse"] = os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')
            if self.option("level").lower() == "transcript":
                function_infos["export_expression"]["trans_diff_output"] = self.diffexpress_trans.tool.output_dir
                function_infos["export_expression"]["trans_uniform_output"] = self.diffexpress_trans.uniform.output_dir
            function_infos["export_expression"]["gene_diff_output"] = self.diffexpress_gene.tool.output_dir
            function_infos["export_expression"]["gene_uniform_output"] = self.diffexpress_gene.uniform.output_dir
        return function_infos


    def set_db(self):
        if self.task_id == "8503_a6g31n9navf21qhl3lu7ml":
            self.IMPORT_REPORT_DATA = True
            self.IMPORT_REPORT_AFTER_END = False
            import datetime
            class DateEncoder(json.JSONEncoder):
                def default(self, obj):
                    if isinstance(obj, datetime.datetime):
                        return obj.strftime("%Y-%m-%d %H:%M:%S")
                    elif isinstance(obj, datetime.date):
                        return obj.strftime("%Y-%m-%d")
                    else:
                        return json.JSONEncoder.default(self, obj)
            with open(os.path.join(self.work_dir, "workflow_sheet_data.json"), 'w') as json_f:
                json.dump(self.sheet.data, json_f, sort_keys=True, indent=4, cls=DateEncoder)
            with open(os.path.join(self.work_dir, "option_data.json"), 'w') as json_f:
                options = self.get_options()
                json.dump(options, json_f, sort_keys=True, indent=4, cls=DateEncoder)
            # self.stop_timeout_check()
            with open(os.path.join(self.work_dir, "function_infos.json"), 'w') as json_f:
                function_infos = self.get_function_infos()
                json.dump(function_infos, json_f, sort_keys=True, indent=4, cls=DateEncoder)
            if self._sheet.rerun:
                self.logger.info("该项目是重运行的,SetDb_rerun.pk文件!")
                tool_dir = os.path.join(self.work_dir, "SetDb")
                if os.path.exists(tool_dir):
                    with open(os.path.join(tool_dir, "rerun.pk"), "w") as w:
                        pass
            self.set_db_api = self.add_tool("denovo_rna_v3.large.set_db")
            options = {
                "sheet_data_json": os.path.join(self.work_dir, "workflow_sheet_data.json"),
                "option_data_json": os.path.join(self.work_dir, "option_data.json"),
                "function_json": os.path.join(self.work_dir, "function_infos.json"),
                "task_id": self.task_id,
                # "analysis_content": str(self.analysis_content),
                # "is_assemble": self.option("is_assemble"),
                "group": self.option("group")
            }
            if self.option('productive_table').is_set:
                options.update({'productive_table': self.option("productive_table").prop["path"]})
            if self.option('datatype') == 'cleandata':
                options.update({'qc_dir': self.option("qc_dir")})
            else:
                options.update({'fastq_dir': self.option("fastq_dir")})
            if self.option('sample_num') == 'multiple':
                options.update({'control': self.option("control")})

            # eval(self.set_db_api_name).set_options(options)
            # eval(self.set_db_api_name).on('end', self.set_upload)
            # eval(self.set_db_api_name).run()
            # )
            self.set_db_api.set_options(options)
            self.set_db_api.on('end', self.set_upload)
            self.set_db_api.run()


        else:
            self.stop_timeout_check()
            self.build_seq_database()  # 创建序列数据库
            # self.merge_annotation_exp_matrix() # 表达量表增加注释信息
            # if self.option("sample_num") == "multiple":
            #     self.merge_annotation_diffexp_matrix() # 差异表达量表增加注释信息
            self.run_api()  # 运行导表函数
            self.set_upload()

    def set_upload(self):
        self.merge_annotation_exp_matrix()  # 表达量表增加注释信息
        if self.option('sample_num') == 'multiple':
            self.merge_annotation_diffexp_matrix()  # 差异表达量表增加注释信息
        # 页面交互运行所需中间结果文件，不在结果目录中呈现
        self.intermediate_result = os.path.join(self.work_dir, 'intermediate_results')
        if os.path.isdir(self.intermediate_result):
            shutil.rmtree(self.intermediate_result)
        os.mkdir(self.intermediate_result)
        ## QC
        # 文件上传到页面后，需要修改文件的内容，才能用于交互分析
        # old_name = self.work_dir + "/HiseqQc/output/sickle_dir/fq_list.txt"
        if self.option('datatype') == 'rawdata':
            old_name = self.qc.option('fq_list').prop['path']
            a = open(old_name, "r").read()
            start_dir = self.work_dir + "/HiseqQc/output/sickle_dir/"
            end_dir = self.workflow_output + "/QC/cleandata/"
            b = a.replace(start_dir, end_dir)
            new_name = self.work_dir + "/HiseqQc/output/sickle_dir/tmp.txt"
            if self.option('qc_soft') == "fastp":
                new_name = self.work_dir + "/FastpRna/output/fastq/tmp.txt"
            with open(new_name, "w") as f:
                f.write(b)
            os.remove(old_name)
            os.rename(new_name, old_name)
            os.mkdir(os.path.join(self.intermediate_result, 'QC'))
            os.link(old_name, os.path.join(self.intermediate_result, 'QC', os.path.basename(old_name)))
        ## SequenceDatabase
        os.mkdir(os.path.join(self.intermediate_result, 'SequenceDatabase'))
        seq_db = os.path.join(self.work_dir, 'seq_db.sqlite3')
        os.link(seq_db, os.path.join(self.intermediate_result, 'SequenceDatabase', os.path.basename(seq_db)))

        os.mkdir(os.path.join(self.intermediate_result, 'SequenceDetail'))
        cds_seq = os.path.join(self.sequence_detail.work_dir, 'cds_seq')
        pep_seq = os.path.join(self.sequence_detail.work_dir, 'pep_seq')
        txpt_seq = os.path.join(self.sequence_detail.work_dir, 'txpt_seq')
        gene_detail = os.path.join(self.sequence_detail.work_dir, 'gene_detail')
        trans_detail = os.path.join(self.sequence_detail.work_dir, 'tran_detail')
        gene_stat = os.path.join(self.sequence_detail.work_dir, 'gene_stat')
        trans_stat = os.path.join(self.sequence_detail.work_dir, 'tran_stat')
        os.link(cds_seq, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(cds_seq)))
        os.link(pep_seq, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(pep_seq)))
        os.link(txpt_seq, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(txpt_seq)))
        os.link(gene_detail, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(gene_detail)))
        os.link(trans_detail, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(trans_detail)))
        os.link(gene_stat, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(gene_stat)))
        os.link(trans_stat, os.path.join(self.intermediate_result, 'SequenceDetail', os.path.basename(trans_stat)))
        ##Assemble
        os.mkdir(os.path.join(self.intermediate_result, 'Assemble'))
        assemble_fasta = self.assemble_filter.option("filter_fa").prop["path"]
        os.link(assemble_fasta, os.path.join(self.intermediate_result, 'Assemble', os.path.basename(assemble_fasta)))
        final_t2g = os.path.join(self.cds_predict.output_dir, "all_tran2gen.txt")
        os.link(final_t2g, os.path.join(self.intermediate_result, 'Assemble', os.path.basename(final_t2g)))
        assemble_unigene_fasta = self.assemble_filter.option("unigene_filter_fa").prop['path']
        os.link(assemble_unigene_fasta,
                os.path.join(self.intermediate_result, 'Assemble', os.path.basename(assemble_unigene_fasta)))
        ## Annotation
        os.makedirs(os.path.join(self.intermediate_result, 'Annotation'))
        CopyFile().linkdir(os.path.join(self.output_dir, 'annotation'),
                           os.path.join(self.intermediate_result, 'Annotation'))
        rm_files = glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/*/*.html.mark')) + \
                   glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/*/*.KOs.txt')) + \
                   glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/*/*.pdf')) + \
                   glob.glob(os.path.join(self.intermediate_result, 'Annotation/cog/summary.*.tsv')) + \
                   glob.glob(os.path.join(self.intermediate_result, 'Annotation/kegg/kegg_layer_*.xls'))
        for rm_file in rm_files:
            os.remove(rm_file)
        os.makedirs(os.path.join(self.intermediate_result, 'Annotation', "blast_xml"))
        blast_nr = os.path.join(self.diamond.output_dir, "nr", "blast.xml")
        os.link(blast_nr, self.intermediate_result + "/Annotation/blast_xml/Trinity_vs_nr.xml")
        blast_swiss = os.path.join(self.diamond.output_dir, "swissprot", "blast.xml")
        os.link(blast_swiss, self.intermediate_result + "/Annotation/blast_xml/Trinity_vs_swissprot.xml")
        blast_string = os.path.join(self.diamond.output_dir, "eggnog", "blast.xml")
        os.link(blast_string, self.intermediate_result + "/Annotation/blast_xml/Trinity_vs_string.xml")
        blast_kegg = os.path.join(self.diamond.output_dir, "kegg", "blast.xml")
        os.link(blast_kegg, self.intermediate_result + "/Annotation/blast_xml/Trinity_vs_kegg.xml")
        blast_pfam = glob.glob(self.cds_predict.output_dir + "/*pfam_domain")[0]
        os.link(blast_pfam, self.intermediate_result + "/Annotation/blast_xml/pfam_domain")
        blast_go = glob.glob(self.diamond.output_dir + "/GO/*blast2go_merge.xls")[0]
        os.link(blast_go, self.intermediate_result + "/Annotation/blast_xml/blast2go_merge.xls")

        ## Express
        if self.option('express_method') == 'RSEM':
            os.makedirs(os.path.join(self.intermediate_result, 'Express/ExpAnnalysis'))
            os.link(os.path.join(self.align.output_dir, 'gene.tpm.matrix'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.xls'))
            os.link(os.path.join(self.align.output_dir, 'gene.count.matrix'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.xls'))
            os.link(os.path.join(self.align.output_dir, 'gene.tpm.matrix.annot.xls'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.annot.xls'))
            os.link(os.path.join(self.align.output_dir, 'gene.count.matrix.annot.xls'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.annot.xls'))
            os.link(os.path.join(self.align.output_dir, 'gene.fpkm.matrix'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.fpkm.matrix.xls'))
            os.link(os.path.join(self.align.output_dir, 'gene.fpkm.matrix.annot.xls'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.fpkm.matrix.annot.xls'))
            if self.option('level').lower() == 'transcript':
                os.link(os.path.join(self.align.output_dir, 'transcript.tpm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.tpm.matrix.xls'))
                os.link(os.path.join(self.align.output_dir, 'transcript.count.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.count.matrix.xls'))
                os.link(os.path.join(self.align.output_dir, 'transcript.tpm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                os.link(os.path.join(self.align.output_dir, 'transcript.count.matrix.annot.xls'),
                        os.path.join(self.intermediate_result,
                                     'Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))
                os.link(os.path.join(self.align.output_dir, 'transcript.fpkm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.fpkm.matrix.xls'))
                os.link(os.path.join(self.align.output_dir, 'transcript.fpkm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls'))
        else:
            os.makedirs(os.path.join(self.intermediate_result, 'Express/ExpAnnalysis'))
            os.link(os.path.join(self.express.output_dir, 'gene.tpm.matrix'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.xls'))
            os.link(os.path.join(self.express.output_dir, 'gene.count.matrix'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.xls'))
            os.link(os.path.join(self.express.output_dir, 'gene.tpm.matrix.annot.xls'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.annot.xls'))
            os.link(os.path.join(self.express.output_dir, 'gene.count.matrix.annot.xls'),
                    os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.annot.xls'))
            if self.option('level').lower() == 'transcript':
                os.link(os.path.join(self.express.output_dir, 'transcript.tpm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.tpm.matrix.xls'))
                os.link(os.path.join(self.express.output_dir, 'transcript.count.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.count.matrix.xls'))
                os.link(os.path.join(self.express.output_dir, 'transcript.tpm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                os.link(os.path.join(self.express.output_dir, 'transcript.count.matrix.annot.xls'),
                        os.path.join(self.intermediate_result,
                                     'Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))
        if self.option('sample_num') == 'multiple':
            ## Diffexpress
            os.makedirs(os.path.join(self.intermediate_result, 'Diffexpress', 'unigene'))
            diff_output = self.diffexpress_gene.output_dir
            files = glob.glob(diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(
                diff_output + '/' + '*_vs_*.sizeFactor.xls')
            for file in files:
                os.link(file, os.path.join(self.intermediate_result, 'Diffexpress', 'unigene', os.path.basename(file)))
            if self.option("level") == "transcript":
                os.makedirs(os.path.join(self.intermediate_result, 'Diffexpress', 'transcript'))
                diff_output = self.diffexpress_trans.output_dir
                files = glob.glob(diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(
                    diff_output + '/' + '*_vs_*.sizeFactor.xls')
                for file in files:
                    os.link(file,
                            os.path.join(self.intermediate_result, 'Diffexpress', 'transcript', os.path.basename(file)))

        ## SNP
        if self.option('sample_num') == 'multiple':
            os.makedirs(os.path.join(self.intermediate_result, 'SNP'))
            if self.option("is_snp") == "True":
                if self.option('snp_method').lower() == 'samtools':
                    final_vcf = os.path.join(self.snp.work_dir, "VcfFilterSamtools/output/final.vcf")
                    os.link(final_vcf, os.path.join(self.intermediate_result, 'SNP', os.path.basename(final_vcf)))
                elif self.option('snp_method').lower() == 'gatk':
                    if os.path.exists(
                            os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")):
                        filted_snp_vcf = os.path.join(self.snp.work_dir,
                                                      "VcfFilterGatk/output/pop.snp.filter.recode.vcf")
                    else:
                        filted_snp_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.snp.filter.recode.vcf")
                    if os.path.exists(
                            os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")):
                        filted_indel_vcf = os.path.join(self.snp.work_dir,
                                                        "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
                    else:
                        filted_indel_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.indel.filter.recode.vcf")
                    os.link(filted_snp_vcf,
                            os.path.join(self.intermediate_result, 'SNP', os.path.basename(filted_snp_vcf)))
                    os.link(filted_indel_vcf,
                            os.path.join(self.intermediate_result, 'SNP', os.path.basename(filted_indel_vcf)))
                else:
                    if os.path.exists(
                            os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")):
                        filted_snp_vcf = os.path.join(self.snp.work_dir,
                                                      "VcfFilterGatk/output/pop.snp.filter.recode.vcf")
                    else:
                        filted_snp_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.snp.filter.recode.vcf")
                    if os.path.exists(
                            os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")):
                        filted_indel_vcf = os.path.join(self.snp.work_dir,
                                                        "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
                    else:
                        filted_indel_vcf = os.path.join(self.snp.vcffilter.output_dir, "pop.indel.filter.recode.vcf")
                    os.link(filted_snp_vcf,
                            os.path.join(self.intermediate_result, 'SNP', os.path.basename(filted_snp_vcf)))
                    os.link(filted_indel_vcf,
                            os.path.join(self.intermediate_result, 'SNP', os.path.basename(filted_indel_vcf)))
                os.makedirs(os.path.join(os.path.join(self.intermediate_result, 'Bam')))
                files = glob.glob(self.snp.work_dir + "/bwa/*")
                for file in files:
                    os.link(file, os.path.join(self.intermediate_result, 'Bam', os.path.basename(file)))

        # Cds
        os.makedirs(os.path.join(self.intermediate_result, 'CDS'))
        bed_path = os.path.join(self.cds_predict.output_dir, "all_predicted.bed")
        os.link(bed_path, os.path.join(self.intermediate_result, 'CDS', os.path.basename(bed_path)))
        transdecoder_pep = os.path.join(self.output_dir, "cds_predict", "all_predicted.pep.fa")
        os.link(transdecoder_pep, os.path.join(self.intermediate_result, 'CDS', "all_predicted.pep"))

        # 比对需要
        os.makedirs(os.path.join(self.intermediate_result, "Assembly_dir"))
        if os.path.exists(os.path.join(self.assemble.output_dir, "assemble_raw.fasta")):
            raw_fasta = os.path.join(self.assemble.output_dir, "assemble_raw.fasta")
            os.link(raw_fasta, os.path.join(self.intermediate_result, "Assembly_dir", "assemble_raw.fasta"))
        assemble_fasta = self.assemble_filter.option("filter_fa").prop["path"]
        os.link(assemble_fasta, os.path.join(self.intermediate_result, 'Assembly_dir', "filter.fasta"))
        final_t2g = os.path.join(self.cds_predict.output_dir, "all_tran2gen.txt")
        os.link(final_t2g, os.path.join(self.intermediate_result, 'Assembly_dir', "all_tran2gen.txt"))
        transdecoder_pep = os.path.join(self.output_dir, "cds_predict", "all_predicted.pep.fa")
        os.link(transdecoder_pep, os.path.join(self.intermediate_result, 'Assembly_dir', "all_predicted.pep"))
        assemble_unigene_fasta = self.assemble_filter.option("unigene_filter_fa").prop['path']
        os.link(assemble_unigene_fasta, os.path.join(self.intermediate_result, 'Assembly_dir', "filter.unigene.fasta"))

        # 项目结果目录中呈现文件
        self.output_dir = self.output_dir
        self.target_dir = os.path.join(self.work_dir, 'upload')
        if os.path.isdir(self.target_dir):
            shutil.rmtree(self.target_dir)
        os.mkdir(self.target_dir)
        # 01Background
        os.mkdir(os.path.join(self.target_dir, '01Background'))
        ## software_info
        software_info = os.path.join(self.target_dir, '01Background', "software_info.xls")
        self.db = Config().get_mongo_client(mtype='denovo_rna_v2')[Config().get_mongo_dbname('denovo_rna_v2')]
        my_collection = self.db['sg_software_database']
        my_results = my_collection.find({})
        with open(software_info, "w") as w:
            w.write("\t".join(["Soft/Database", "Version", "Analysis", "Source"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["software_database"]), str(collection["version"]), str(collection["usage"]),
                     str(collection["source"])]) + "\n")
        ## get run parameter
        if self.option("get_run_log"):
            self.run_parameter(os.path.join(self.target_dir, '01Background'))
        ## sample_info
        sample_info = os.path.join(self.target_dir, '01Background', "sample_info.xls")
        productive_names = dict()
        mj_names = dict()
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
        with open(sample_info, "w") as w, open(self.option('group').prop["path"], "r") as f:
            if self.option("productive_table").is_set:
                if mj_names:
                    w.write("\t".join(["Sample Productive Name", "MJ_name", "Sample Name", "Group Name"]) + "\n")
                else:
                    w.write("\t".join(["Sample Productive Name", "Sample Name", "Group Name"]) + "\n")
            else:
                w.write("\t".join(["Sample Name", "Group Name"]) + "\n")
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    if line.strip().split("\t")[0] in productive_names:
                        if line.strip().split("\t")[0] in mj_names:
                            w.write("\t".join([productive_names[line.strip().split("\t")[0]], mj_names[line.strip().split("\t")[0]],
                                               line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")
                        else:
                            w.write("\t".join([productive_names[line.strip().split("\t")[0]],
                                           line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")
                    else:
                        w.write("\t".join([line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")

        # 02QC
        os.mkdir(os.path.join(self.target_dir, '02QC'))
        if self.option('datatype') == 'rawdata':
            raw_stat = os.path.join(self.output_dir, 'QC_stat/before_qc/fastq_stat.xls')
            os.link(raw_stat, os.path.join(self.target_dir, '02QC/rawdata_statistics.xls'))
        clean_stat = os.path.join(self.output_dir, 'QC_stat/after_qc/fastq_stat.xls')
        os.link(clean_stat, os.path.join(self.target_dir, '02QC/cleandata_statistics.xls'))

        # 03Assemble
        os.mkdir(self.target_dir + "/03Assemble")
        os.mkdir(self.target_dir + "/03Assemble/01Assemble_Sequence")
        if self.option("assemble") == False:
            orgin_trinity = self.option("assembly_file").prop['path']
        else:
            orgin_trinity = os.path.join(self.assemble.output_dir + "/assemble_raw.fasta")
        if os.path.exists(orgin_trinity):
            os.link(orgin_trinity,
                    self.target_dir + "/03Assemble/01Assemble_Sequence/{}.fasta".format(self.option("assemble_soft")))
        if self.option("assemble") == False:
            orgin_trinity_g2t = self.option("gene_to_trans").prop['path']
        else:
            orgin_trinity_g2t = os.path.join(self.assemble.output_dir + "/assemble_raw.gene_trans_map")
        if os.path.exists(orgin_trinity_g2t):
            os.link(orgin_trinity_g2t, self.target_dir + "/03Assemble/01Assemble_Sequence/{}.gene_trans_map".format(
                self.option("assemble_soft")))
        orgin_trinity_filter = os.path.join(self.output_dir + "/assemble/Trinity.filter.fasta")
        os.link(orgin_trinity_filter, self.target_dir + "/03Assemble/01Assemble_Sequence/{}.filter.fasta".format(
            self.option("assemble_soft")))
        if os.path.exists(self.output_dir + "/assemble/Trinity.filter.g2t"):
            orgin_trinity_filter_g2t = self.output_dir + "/assemble/Trinity.filter.g2t"
        else:
            orgin_trinity_filter_g2t = os.path.join(self.output_dir + "/assemble/Trinity.filter.gene_trans_map")
        os.link(orgin_trinity_filter_g2t,
                self.target_dir + "/03Assemble/01Assemble_Sequence/{}.filter.gene_trans_map".format(
                    self.option("assemble_soft")))
        orgin_unigene_filter = os.path.join(self.output_dir + "/assemble/Trinity.filter.unigene.fasta")
        os.link(orgin_unigene_filter,
                self.target_dir + "/03Assemble/01Assemble_Sequence/{}.filter.unigene.fasta".format(
                    self.option("assemble_soft")))
        os.mkdir(self.target_dir + "/03Assemble/02Assemble_Evaluation")
        before_t = os.path.join(self.output_dir + "/assemble/trinity_stat/Trinity_stat.xls")
        # os.link(before_t, self.target_dir + "/03Assemble/AssembleEvaluation/Optimize_before/Trinity_stat.xls")
        if self.option("level").lower() == "transcript":
            before_g = os.path.join(self.unigene_evaluation.output_dir, "trinity_stat/Trinity_stat.xls")
            bg = pd.read_table(before_g, index_col=0).rename(columns={"Resource": "unigene"})
            bt = pd.read_table(before_t, index_col=0).rename(columns={"Resource": "transcript"})
            bs = pd.concat([bg, bt], axis=1)
            bs.to_csv(
                self.target_dir + "/03Assemble/02Assemble_Evaluation/{}_stat.xls".format(self.option("assemble_soft")),
                sep="\t")
        else:
            os.link(before_t, self.target_dir + "/03Assemble/02Assemble_Evaluation/{}_stat.xls".format(
                self.option("assemble_soft")))
        if self.option("optimize") == True:
            after_t = os.path.join(self.output_dir + "/filter_evaluation/trinity_stat/Trinity_stat.xls")
            if self.option("level").lower() == "transcript":
                after_g = os.path.join(self.filter_unigene_evaluation.output_dir, "trinity_stat/Trinity_stat.xls")
                ag = pd.read_table(after_g, index_col=0).rename(columns={"Resource": "unigene"})
                at = pd.read_table(after_t, index_col=0).rename(columns={"Resource": "transcript"})
                ast = pd.concat([ag, at], axis=1)
                ast.to_csv(self.target_dir + "/03Assemble/02Assemble_Evaluation/{}_filter_stat.xls".format(
                    self.option("assemble_soft")), sep="\t")
            else:
                os.link(after_t, self.target_dir + "/03Assemble/02Assemble_Evaluation/{}_filter_stat.xls".format(
                    self.option("assemble_soft")))
        os.mkdir(self.target_dir + "/03Assemble/03Length_Distribution")
        if self.option("optimize") == False:
            gene_len_dis = os.path.join(self.output_dir + "/assemble/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, self.target_dir + "/03Assemble/03Length_Distribution/unigene_count_stat_500.xls")
            if self.option("level").lower() == "transcript":
                trans_len_dis = os.path.join(self.output_dir + "/assemble/trinity_stat/trans_count_stat_500.txt")
                os.link(trans_len_dis,
                        self.target_dir + "/03Assemble/03Length_Distribution/transcript_count_stat_500.xls")
        else:
            gene_len_dis = os.path.join(self.output_dir + "/filter_evaluation/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, self.target_dir + "/03Assemble/03Length_Distribution/unigene_count_stat_500.xls")
            if self.option("level").lower() == "transcript":
                trans_len_dis = os.path.join(
                    self.output_dir + "/filter_evaluation/trinity_stat/trans_count_stat_500.txt")
                os.link(trans_len_dis,
                        self.target_dir + "/03Assemble/03Length_Distribution/transcript_count_stat_500.xls")
        # os.mkdir(self.target_dir + "/03Assemble/Mapping")
        # align_result = os.path.join(self.align.output_dir + "/alignment_rate.txt")
        # os.link(align_result, self.target_dir + "/03Assemble/Mapping/alignment_rate.txt")
        align_result = os.path.join(self.align.output_dir + "/alignment_rate.txt")
        os.link(align_result, self.target_dir + "/03Assemble/alignment_rate.xls")

        ##04Annotation
        os.mkdir(self.target_dir + "/04Annotation")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/GO")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/NR")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/COG")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/KEGG")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/Pfam")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/Swiss_Prot")
        # all
        all_anno_stat = os.path.join(self.annot_stat.output_dir, "all_stat_detail.xls")
        os.link(all_anno_stat, self.target_dir + "/04Annotation/all_annotation_statistics.xls")
        # Unigene_Anno
        # GO
        go_list = os.path.join(self.output_dir, "annotation", "go", "go_list_gene.xls")
        os.link(go_list, self.target_dir + "/04Annotation/Unigene_Anno/GO/unigene_gos.list")
        go_stat1 = os.path.join(self.output_dir, "annotation", "go", "go_lev2_gene.stat.xls")
        os.link(go_stat1, self.target_dir + "/04Annotation/Unigene_Anno/GO/unigene_go12level_statistics.xls")
        go_stat2 = os.path.join(self.output_dir, "annotation", "go", "go_lev3_gene.stat.xls")
        os.link(go_stat2, self.target_dir + "/04Annotation/Unigene_Anno/GO/unigene_go123level_statistics.xls")
        go_stat3 = os.path.join(self.output_dir, "annotation", "go", "go_lev4_gene.stat.xls")
        os.link(go_stat3, self.target_dir + "/04Annotation/Unigene_Anno/GO/unigene_go1234level_statistics.xls")
        # NR
        NR_detail = os.path.join(self.output_dir, "annotation", "nr", "nr_blast_gene.xls")
        anno_nr = pd.read_table(NR_detail)
        anno_nr.rename(columns={"Query-Name": "Gene_id", "Hit-Description": "Description", "HSP-Len": "HSP-Length",
                                "Identity-%": "Identity(%)", "Similarity-%": "Similarity(%)"}, inplace=True)
        anno_nr = anno_nr[
            ["Gene_id", "Hit-Name", "Description", "HSP-Length", "E-Value", "Score", "Identity(%)", "Similarity(%)"]]
        anno_nr.to_csv(self.target_dir + "/04Annotation/Unigene_Anno/NR/unigene_nr_anno_detail.xls", sep="\t", index=0)
        # os.link(NR_detail, self.target_dir + "/04Annotation/Unigene_Anno/NR/unigene_nr_anno_detail.xls")
        # COG
        COG_detail = os.path.join(self.output_dir, "annotation", "cog", "summary.G.tsv")
        os.link(COG_detail, self.target_dir + "/04Annotation/Unigene_Anno/COG/unigene_cog_summary.xls")
        # KEGG
        kegg_path = os.path.join(self.output_dir, "annotation", "kegg", "kegg_pathway_gene.xls")
        os.link(kegg_path, self.target_dir + "/04Annotation/Unigene_Anno/KEGG/unigene_pathways_table.xls")
        pathway_level = os.path.join(self.output_dir, "annotation", "kegg", "kegg_layer_gene.xls")
        os.link(pathway_level, self.target_dir + "/04Annotation/Unigene_Anno/KEGG/unigene_pathway_level.xls")
        with open(self.target_dir + "/04Annotation/Unigene_Anno/KEGG/unigene_pathway_level.xls", 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('First_Category\tSecond_Category\tnumber\tseq_list\n' + content)
        kegg_table = os.path.join(self.output_dir, "annotation", "kegg", "kegg_gene_gene.xls")
        os.link(kegg_table, self.target_dir + "/04Annotation/Unigene_Anno/KEGG/unigene_kegg_table.xls")
        os.mkdir(self.target_dir + "/04Annotation/Unigene_Anno/KEGG/unigene_pathways")
        pathway_files = glob.glob(os.path.join(self.output_dir, 'annotation/kegg/kegg_pathway_gene_dir/*.html')) + \
                        glob.glob(os.path.join(self.output_dir, 'annotation/kegg/kegg_pathway_gene_dir/*.png'))
        for file in pathway_files:
            os.link(file, os.path.join(self.target_dir, "04Annotation/Unigene_Anno/KEGG/unigene_pathways",
                                       os.path.basename(file)))
        # Pfam
        pfam_path = os.path.join(self.output_dir, "annotation", "pfam", "pfam_domain_gene.xls")
        os.link(pfam_path, self.target_dir + "/04Annotation/Unigene_Anno/Pfam/unigene_pfam_anno_detail.xls")
        # swissprot
        swiss_path = os.path.join(self.output_dir, "annotation", "swissprot", "swissprot_blast_gene.xls")
        os.link(swiss_path, self.target_dir + "/04Annotation/Unigene_Anno/Swiss_Prot/unigene_swissprot_anno_detail.xls")

        # Transcript_Anno
        # os.mkdir(self.target_dir + "/04Annotation")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/GO")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/NR")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/COG")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/KEGG")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/Pfam")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/Swiss_Prot")
        # GO
        go_list = os.path.join(self.output_dir, "annotation", "go", "go_list_tran.xls")
        os.link(go_list, self.target_dir + "/04Annotation/Transcript_Anno/GO/transcript_gos.list")
        go_stat1 = os.path.join(self.output_dir, "annotation", "go", "go_lev2_tran.stat.xls")
        os.link(go_stat1, self.target_dir + "/04Annotation/Transcript_Anno/GO/transcript_go12level_statistics.xls")
        go_stat2 = os.path.join(self.output_dir, "annotation", "go", "go_lev3_tran.stat.xls")
        os.link(go_stat2, self.target_dir + "/04Annotation/Transcript_Anno/GO/transcript_go123level_statistics.xls")
        go_stat3 = os.path.join(self.output_dir, "annotation", "go", "go_lev4_tran.stat.xls")
        os.link(go_stat3, self.target_dir + "/04Annotation/Transcript_Anno/GO/transcript_go1234level_statistics.xls")
        # NR
        NR_detail = os.path.join(self.output_dir, "annotation", "nr", "nr_blast_tran.xls")
        annot_output = self.annot_stat.output_dir
        all_annot = os.path.join(annot_output, 'all_annot.xls')
        all_annot = pd.read_table(all_annot, header=0)
        pattern = '.*?[(](.*?)[)]'
        all_annot["nr_description"] = all_annot["nr"].str.extract(pattern)
        trans_info_pd = all_annot[['transcript', 'gene_id', 'nr_description']].set_index('transcript')
        anno_nrt = pd.read_table(NR_detail)
        anno_nrt.rename(
            columns={"Query-Name": "Transcript_id", "Hit-Description": "Description", "HSP-Len": "HSP-Length",
                     "Identity-%": "Identity(%)", "Similarity-%": "Similarity(%)"}, inplace=True)
        anno_nrt = anno_nrt[
            ["Transcript_id", "Hit-Name", "Description", "HSP-Length", "E-Value", "Score", "Identity(%)",
             "Similarity(%)"]].set_index("Transcript_id")
        nrt_final = pd.concat([trans_info_pd, anno_nrt], axis=1, join_axes=[anno_nrt.index])
        nrt_final.to_csv(self.target_dir + "/04Annotation/Transcript_Anno/NR/transcript_nr_anno_detail.xls", sep="\t")

        # COG
        COG_detail = os.path.join(self.output_dir, "annotation", "cog", "summary.T.tsv")
        os.link(COG_detail, self.target_dir + "/04Annotation/Transcript_Anno/COG/transcript_cog_summary.xls")
        # KEGG
        kegg_path = os.path.join(self.output_dir, "annotation", "kegg", "kegg_pathway_tran.xls")
        os.link(kegg_path, self.target_dir + "/04Annotation/Transcript_Anno/KEGG/transcript_pathways_table.xls")
        pathway_level = os.path.join(self.output_dir, "annotation", "kegg", "kegg_layer_tran.xls")
        os.link(pathway_level, self.target_dir + "/04Annotation/Transcript_Anno/KEGG/transcript_pathway_level.xls")
        with open(self.target_dir + "/04Annotation/Transcript_Anno/KEGG/transcript_pathway_level.xls", 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('First_Category\tSecond_Category\tnumber\tseq_list\n' + content)
        kegg_table = os.path.join(self.output_dir, "annotation", "kegg", "kegg_gene_tran.xls")
        os.link(kegg_table, self.target_dir + "/04Annotation/Transcript_Anno/KEGG/transcript_kegg_table.xls")
        os.mkdir(self.target_dir + "/04Annotation/Transcript_Anno/KEGG/transcript_pathways")
        pathway_files = glob.glob(os.path.join(self.output_dir, 'annotation/kegg/kegg_pathway_tran_dir/*.html')) + \
                        glob.glob(os.path.join(self.output_dir, 'annotation/kegg/kegg_pathway_tran_dir/*.png'))
        for file in pathway_files:
            os.link(file, os.path.join(self.target_dir, "04Annotation/Transcript_Anno/KEGG/transcript_pathways",
                                       os.path.basename(file)))
        # Pfam
        pfam_path = os.path.join(self.output_dir, "annotation", "pfam", "pfam_domain_tran.xls")
        os.link(pfam_path, self.target_dir + "/04Annotation/Transcript_Anno/Pfam/transcript_pfam_anno_detail.xls")
        # swissprot
        swiss_path = os.path.join(self.output_dir, "annotation", "swissprot", "swissprot_blast_tran.xls")
        os.link(swiss_path,
                self.target_dir + "/04Annotation/Transcript_Anno/Swiss_Prot/transcript_swissprot_anno_detail.xls")

        # #Xml文件
        # os.mkdir(self.target_dir + "/04Annotation/blast_xml")
        # blast_nr = os.path.join(self.diamond.output_dir, "nr", "blast.xml")
        # os.link(blast_nr, self.target_dir + "/04Annotation/blast_xml/Trinity_vs_nr.xml")
        # blast_swiss = os.path.join(self.diamond.output_dir, "swissprot", "blast.xml")
        # os.link(blast_swiss, self.target_dir + "/04Annotation/blast_xml/Trinity_vs_swissprot.xml")
        # blast_string = os.path.join(self.diamond.output_dir, "eggnog", "blast.xml")
        # os.link(blast_string, self.target_dir + "/04Annotation/blast_xml/Trinity_vs_string.xml")
        # blast_kegg = os.path.join(self.diamond.output_dir, "kegg", "blast.xml")
        # os.link(blast_kegg, self.target_dir + "/04Annotation/blast_xml/Trinity_vs_kegg.xml")
        # blast_pfam = glob.glob(self.cds_predict.output_dir + "/*pfam_domain")[0]
        # os.link(blast_pfam, self.target_dir + "/04Annotation/blast_xml/pfam_domain")
        # blast_go = glob.glob(self.diamond.output_dir + "/GO/*blast2go_merge.xls")[0]
        # os.link(blast_go, self.target_dir + "/04Annotation/blast_xml/blast2go_merge.xls")

        ##05AnnoQuery
        os.mkdir(self.target_dir + "/05AnnoQuery")
        anno_query = os.path.join(self.output_dir, "annotation", "all_annot.xls")
        annot_file = RefAnnotation()
        out_gene = self.target_dir + "/05AnnoQuery/unigene_anno_detail.xls"
        out_trans = self.target_dir + "/05AnnoQuery/transcript_anno_detail.xls"
        annot_file.change_annot_result(anno_query, out_gene, out_trans)
        if self.option("level").lower() != "transcript":
            os.remove(out_trans)

        ##06TF
        if self.option("tf_database").lower() != "other":
            os.mkdir(self.target_dir + "/06TF")
            if self.option('level').lower() == 'transcript':
                trans_tf = os.path.join(self.output_dir, "cds_predict", "transcript_tf_detail.xls")
                tran_df = pd.read_table(trans_tf)
                tran_df.rename(columns={"DNA_binding domain": "DNA_binding_domain"}, inplace=True)
                tran_df.rename(
                    columns={"DNA_binding_domain": "DNA domain", "query_name": "Transcript_id", "accession": "PF ID",
                             "description_of_target": "Description"}, inplace=True)
                trans_f = tran_df[["Transcript_id", "PF ID", "DNA domain", "Description", "Family", "E_value", "score"]]
                trans_f.to_csv(self.target_dir + "/06TF/transcript_tf_predict.xls", sep="\t", index=0)
            gene_tf = os.path.join(self.output_dir, "cds_predict", "unigene_tf_detail.xls")
            gene_df = pd.read_table(gene_tf)
            gene_df.rename(columns={"DNA_binding domain": "DNA_binding_domain"}, inplace=True)
            gene_df.rename(columns={"DNA_binding_domain": "DNA domain", "unigene": "Gene_id", "accession": "PF ID",
                                    "description_of_target": "Description"}, inplace=True)
            gene_f = gene_df[["Gene_id", "PF ID", "DNA domain", "Description", "Family", "E_value", "score"]]
            gene_f.to_csv(self.target_dir + "/06TF/unigene_tf_predict.xls", sep="\t", index=0)

        ##07Express
        os.makedirs(os.path.join(self.target_dir, '07Express/ExpAnnalysis'))
        if self.option('express_method') == 'RSEM':
            os.link(os.path.join(self.align.output_dir, 'gene.count.matrix.annot.xls'),
                    os.path.join(self.target_dir, '07Express/ExpAnnalysis/unigene.count.matrix.annot.xls'))
            if self.option('exp_way').lower() == 'tpm':
                os.link(os.path.join(self.align.output_dir, 'gene.tpm.matrix.annot.xls'),
                        os.path.join(self.target_dir, '07Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls'))
            if self.option('exp_way').lower() == 'fpkm':
                os.link(os.path.join(self.align.output_dir, 'gene.fpkm.matrix.annot.xls'),
                        os.path.join(self.target_dir, '07Express/ExpAnnalysis/unigene.fpkm.matrix.annot.xls'))
            if self.option('level').lower() == 'transcript':
                os.link(os.path.join(self.align.output_dir, 'transcript.count.matrix.annot.xls'),
                        os.path.join(self.target_dir, '07Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))
                if self.option('exp_way').lower() == 'tpm':
                    os.link(os.path.join(self.align.output_dir, 'transcript.tpm.matrix.annot.xls'),
                            os.path.join(self.target_dir, '07Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                if self.option('exp_way').lower() == 'fpkm':
                    os.link(os.path.join(self.align.output_dir, 'transcript.fpkm.matrix.annot.xls'),
                            os.path.join(self.target_dir, '07Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls'))
        else:
            os.link(os.path.join(self.express.output_dir, 'gene.tpm.matrix.annot.xls'),
                    os.path.join(self.target_dir, '07Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls'))
            os.link(os.path.join(self.express.output_dir, 'gene.count.matrix.annot.xls'),
                    os.path.join(self.target_dir, '07Express/ExpAnnalysis/unigene.count.matrix.annot.xls'))
            if self.option('level').lower() == 'transcript':
                os.link(os.path.join(self.express.output_dir, 'transcript.tpm.matrix.annot.xls'),
                        os.path.join(self.target_dir, '07Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                os.link(os.path.join(self.express.output_dir, 'transcript.count.matrix.annot.xls'),
                        os.path.join(self.target_dir, '07Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))
        if self.option('sample_num') == 'multiple':
            os.makedirs(os.path.join(self.target_dir, '07Express/ExpCorr'))
            os.link(os.path.join(self.exp_corr.output_dir, 'sample_correlation.xls'),
                    os.path.join(self.target_dir, '07Express/ExpCorr/sample_correlation.xls'))
            if self.option('group').prop['sample_number'] > 2:
                os.makedirs(os.path.join(self.target_dir, '07Express/ExpPCA'))
                os.link(os.path.join(self.exp_pca.output_dir, 'PCA.xls'),
                        os.path.join(self.target_dir, '07Express/ExpPCA/PCA.xls'))
                exppca = pd.read_table(os.path.join(self.exp_pca.output_dir, 'Explained_variance_ratio.xls'),
                                       header=None)
                exppca.columns = ["", "Proportion of Variance"]
                exppca.to_csv(os.path.join(self.target_dir, '07Express/ExpPCA/Explained_variance_ratio.xls'), sep="\t",
                              index=0)
                # os.link(os.path.join(self.exp_pca.output_dir, 'Explained_variance_ratio.xls'),
                #         os.path.join(self.target_dir, '07Express/ExpPCA/Explained_variance_ratio.xls'))

        ##08DiffExpress
        if self.option('sample_num') == 'multiple':
            os.makedirs(os.path.join(self.target_dir, '08DiffExpress/unigene'))
            files = glob.glob(self.diffexpress_gene.output_dir + "/*.annot.xls") + glob.glob(
                self.diffexpress_gene.output_dir + "/diff_summary*")
            for b in files:
                if not os.path.exists(os.path.join(self.target_dir, '08DiffExpress/unigene', os.path.basename(b))):
                    os.link(b, os.path.join(self.target_dir, '08DiffExpress/unigene', os.path.basename(b)))

            if self.option("level").lower() == "transcript":
                os.makedirs(os.path.join(self.target_dir, '08DiffExpress/transcript'))
                files = glob.glob(self.diffexpress_trans.output_dir + "/*.annot.xls") + glob.glob(
                    self.diffexpress_trans.output_dir + "/diff_summary*")
                for b in files:
                    if not os.path.exists(
                            os.path.join(self.target_dir, '08DiffExpress/transcript', os.path.basename(b))):
                        os.link(b, os.path.join(self.target_dir, '08DiffExpress/transcript', os.path.basename(b)))

        ##09CDS
        os.makedirs(os.path.join(self.target_dir, '09CDS'))
        # pfam_domain=os.path.join(self.output_dir,"cds_predict","pfam_domain")
        # os.link(pfam_domain, os.path.join(self.target_dir, '09CDS', os.path.basename(pfam_domain)))
        if self.option("level").lower() == "transcript":
            tran_cds = os.path.join(self.output_dir, "cds_predict", "cds_len_transcript.txt")
            os.link(tran_cds, os.path.join(self.target_dir, '09CDS', "transcript_cds_len.xls"))
        gene_cds = os.path.join(self.output_dir, "cds_predict", "cds_len_unigene.txt")
        os.link(gene_cds, os.path.join(self.target_dir, '09CDS', "unigene_cds_len.xls"))
        all_cds = os.path.join(self.output_dir, "cds_predict", "all_predicted.xls")
        os.link(all_cds, os.path.join(self.target_dir, '09CDS', os.path.basename(all_cds)))
        transdecoder_bed = os.path.join(self.output_dir, "cds_predict", "all_predicted.bed")
        os.link(transdecoder_bed, os.path.join(self.target_dir, '09CDS', os.path.basename(transdecoder_bed)))
        transdecoder_cds = os.path.join(self.output_dir, "cds_predict", "all_predicted.cds.fa")
        os.link(transdecoder_cds, os.path.join(self.target_dir, '09CDS', "all_predicted.cds"))
        transdecoder_pep = os.path.join(self.output_dir, "cds_predict", "all_predicted.pep.fa")
        os.link(transdecoder_pep, os.path.join(self.target_dir, '09CDS', "all_predicted.pep"))

        ##10SNP
        if self.option('sample_num') == 'multiple':
            os.makedirs(os.path.join(self.target_dir, '10SNP'))
            os.makedirs(os.path.join(self.target_dir, 'others'))
            if self.option("is_snp") == "True":
                if self.option('snp_method').lower() == 'samtools':
                    final_vcf = os.path.join(self.snp.work_dir, "VcfFilterSamtools/output/final.vcf")
                    os.link(final_vcf, os.path.join(self.target_dir, 'others', "final.vcf"))
                elif self.option('snp_method').lower() == 'gatk':
                    final_vcf_snp = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")
                    final_vcf_indel = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
                    os.link(final_vcf_snp, os.path.join(self.target_dir, 'others', "filter_snp_vcf"))
                    os.link(final_vcf_indel,os.path.join(self.target_dir, 'others', "filter_indel_vcf"))
                else:
                    if os.path.exists(os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")):
                        final_vcf_snp = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.snp.filter.recode.vcf")
                    else:
                        final_vcf_snp = os.path.join(self.snp.vcffilter.output_dir,
                                                     "pop.snp.filter.recode.vcf")
                    if os.path.exists(os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")):
                        final_vcf_indel = os.path.join(self.snp.work_dir, "VcfFilterGatk/output/pop.indel.filter.recode.vcf")
                    else:
                        final_vcf_indel = os.path.join(self.snp.vcffilter.output_dir,
                                                     "pop.indel.filter.recode.vcf")
                    os.link(final_vcf_snp, os.path.join(self.target_dir, 'others',"filter_snp_vcf"))
                    os.link(final_vcf_indel,os.path.join(self.target_dir, 'others',"filter_indel_vcf"))

                if self.snpfinal.work_dir + '/new_snp_rewrite' is not None:
                    new_snp_rewrite = self.snpfinal.work_dir + '/new_snp_rewrite'
                    depth_path = self.snpfinal.work_dir + '/depth_new_per'
                    hh_path = self.snpfinal.work_dir + '/statis_hh'
                    tt_new_per_path = self.snpfinal.work_dir + '/transition_tranversion'
                    cds_path = self.snpfinal.work_dir + "/statis_cds"
                    anno_path = self.snpfinal.work_dir + "/snp_anno_stat"
                    anno_snp = pd.read_table(anno_path)
                    anno_snp.rename(
                        columns={"sample": "Sample", "genes": "unigene contain SNPs", "totalre": "total_annot SNPs",
                                 "gore": "GO_annot SNPs", "keggre": "KEGG_annot SNPs", "cogre": "COG_annot SNPs",
                                 "nrre": "NR_annot SNPs", "swissre": "Swiss-Prot_annot SNPs"}, inplace=True)
                    anno_snp.to_csv(os.path.join(self.target_dir, '10SNP', "snp_anno_statistics.xls"), sep="\t",
                                    index=0)
                    type_snp = pd.read_table(tt_new_per_path, index_col=0, header=0)
                    type_snp.loc['total'] = type_snp.apply(lambda x: x.sum())
                    type_snp.to_csv(os.path.join(self.target_dir, '10SNP', "snp_transition_tranversion_statistics.xls"),
                                    sep="\t")
                    os.link(new_snp_rewrite, os.path.join(self.target_dir, '10SNP', "snp_detail.xls"))
                    os.link(depth_path, os.path.join(self.target_dir, '10SNP', "snp_depth_statistics.xls"))
                    os.link(hh_path, os.path.join(self.target_dir, '10SNP', "snp_homo_hete_statistics.xls"))
                    # os.link(tt_new_per_path, os.path.join(self.target_dir, '10SNP', "snp_transition_tranversion_statistics.xls"))
                    os.link(cds_path, os.path.join(self.target_dir, '10SNP', "snp_cds_statistics.xls"))

        ##11SSR
        os.makedirs(os.path.join(self.target_dir, '11SSR'))
        ssr_files = glob.glob(self.ssr.output_dir + "/*")
        for file in ssr_files:
            os.link(file, os.path.join(self.target_dir, '11SSR', os.path.basename(file) + ".xls"))
        ## GENESET
        os.makedirs(os.path.join(self.target_dir, '12GeneSet'))
        if self.option('sample_num') == 'multiple' and os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, "cluster")):
            geneset_upload_dirs = {
                "cluster": os.path.join(self.target_dir, '12GeneSet', '01_Cluster_Analysis'),
                "go_class": os.path.join(self.target_dir, '12GeneSet', '03_GO_Annotation'),
                "cog_class": os.path.join(self.target_dir, '12GeneSet', '02_COG_Annotation'),
                "kegg_class": os.path.join(self.target_dir, '12GeneSet', '04_KEGG_Annotation'),
                "go_enrich": os.path.join(self.target_dir, '12GeneSet', '05_GO_Enrich'),
                "kegg_enrich": os.path.join(self.target_dir, '12GeneSet', '06_KEGG_Enrich'),
                "venn": os.path.join(self.target_dir, '12GeneSet', '07_GenesetVenn')
            }
            for paths in geneset_upload_dirs.values():
                if os.path.exists(paths):
                    shutil.rmtree(paths)
                os.makedirs(paths)
            diff_genesets = []
            for gt in os.listdir(os.path.join(self.diff_geneset_analysis.output_dir)):
                if gt not in ["cluster", "prepare_json"]:
                    diff_genesets.append(gt)


            select_geneset = os.listdir(os.path.join(self.diff_geneset_analysis.output_dir,"cluster"))[0]

            shutil.copytree(os.path.join(self.diff_geneset_analysis.output_dir,"cluster",select_geneset),
                            os.path.join(geneset_upload_dirs["cluster"], select_geneset))

            #GO_Annotion
            for geneset in diff_genesets:
                if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset,"diff_cog_class")):
                    os.makedirs(os.path.join(
                        geneset_upload_dirs["cog_class"]
                        ,geneset))
                    raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_cog_class','cog_class_table.xls')
                    file_path = os.path.join(geneset_upload_dirs["cog_class"], geneset,"cog_class_stat.xls")
                    os.link(raw_path, file_path)
                if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset,"diff_go_class")):
                    os.makedirs(os.path.join(
                        geneset_upload_dirs["go_class"]
                        ,geneset))
                    raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_go_class','go_class_table.xls')
                    file_path = os.path.join(geneset_upload_dirs["go_class"], geneset,"go_class_stat.xls")
                    os.link(raw_path, file_path)

                #02KEGG_Annotion
                if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_kegg_class")):
                    os.makedirs(os.path.join(geneset_upload_dirs["kegg_class"],geneset))

                    raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_kegg_class','kegg_stat.xls')
                    file_path = os.path.join(geneset_upload_dirs["kegg_class"], geneset,"kegg_class_stat.xls")
                    os.link(raw_path, file_path)

                    os.link(
                        os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_class', 'pathways.tar.gz'),
                        os.path.join(geneset_upload_dirs["kegg_class"], geneset,
                                        "pathways.tar.gz"))

                # GO_Enrich
                if os.path.exists(
                        os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_go_enrich")):
                    os.makedirs(os.path.join(geneset_upload_dirs["go_enrich"], geneset ))
                    raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_go_enrich',
                                            'go_enrich_geneset_list_gene.xls')
                    file_path = os.path.join(geneset_upload_dirs["go_enrich"], geneset,  "go_enrich_stat.xls")
                    os.link(raw_path,file_path)

                # 02KEGG_Enrich
                if os.path.exists(
                        os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_kegg_enrich")):
                    os.makedirs(os.path.join(geneset_upload_dirs["kegg_enrich"], geneset))
                    raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich','enrich',
                                            '{}_gene.list.DE.list.check.kegg_enrichment.xls'.format(geneset))
                    file_path = os.path.join(geneset_upload_dirs["kegg_enrich"], geneset , "kegg_enrich_stat.xls")
                    os.link(raw_path,file_path)

                    os.link(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich', 'class',
                                         'pathways.tar.gz'),
                            os.path.join(geneset_upload_dirs["kegg_enrich"], geneset,
                                         "pathways.tar.gz"))



        # seq_db = self.work_dir + "/seq_db.sqlite3"
        # os.makedirs(os.path.join(self.target_dir, 'Sequence_database'))
        # os.link(seq_db,self.target_dir + "/Sequence_database/seq_db.sqlite3")
        self.move_chart_file()
        self.def_upload()

    def run_parameter(self, dir_path):
        document = self.db['sg_table_relation'].find_one({})
        targets = document['target']
        parsed_tables = list()
        for target in targets:
            if target[0] in ['sg_status', 'sg_task', 'sg_software_para', 'sg_geneset', 'task_workflow_params',
                             'sg_one_touch', 'sg_file_dowload', 'sg_transcripts', 'sg_result_table_deposit',
                             'sg_splicing_rmats_count', 'sg_annotation_query'] or target[0] in parsed_tables:
                continue
            parsed_tables.append(target[0])
            if target[0] in ['sg_exp']:
                records = self.db[target[0]].find_one({'task_id': self.task_id})
            else:
                records = self.db[target[0]].find({'task_id': self.task_id})
            if type(records) == dict:
                record = records
                if 'params' in record and record['params'] and record['params'] != 'none':
                    if type(record['params']) == dict:
                        params = json.loads(record['params'])
                    else:
                        params = record['params']
                    if 'submit_location' in params:
                        get_run_log = GetRunLog("denovo_rna_v2", table=target[0], main_id=record['main_id'],
                                                dir_path=dir_path, append=True)
                        get_run_log.run()
            else:
                for record in records:
                    if 'params' in record and record['params'] and record['params'] != 'none':
                        if type(record['params']) == dict:
                            params = json.loads(record['params'])
                        else:
                            params = record['params']
                        if 'submit_location' in params:
                            get_run_log = GetRunLog("denovo_rna_v2", table=target[0], main_id=record['main_id'],
                                                    dir_path=dir_path, append=True)
                            get_run_log.run()

    def move_chart_file(self):
        os.makedirs(os.path.join(self.target_dir, '07Express/ExpVenn/'))
        file2uploads = [
            ("*raw_qc_qual.box*.pdf", '02QC/'),
            ("*raw_qc_error.line*.pdf", '02QC/'),
            ("*raw_qc_base.line*.pdf", '02QC/'),
            ("*clean_qc_qual.box*.pdf", '02QC/'),
            ("*clean_qc_error.line*.pdf", '02QC/'),
            ("*clean_qc_base.line*.pdf", '02QC/'),
            ("*length_distribution.*.pdf", '03Assemble/'),
            ("*annot_gene_stat*.pdf", '04Annotation/Unigene_Anno/'),
            ("*annot_tran_stat*.pdf", '04Annotation/Transcript_Anno/'),
            ("*Gene.nr*.pdf", '04Annotation/Unigene_Anno/NR/'),
            ("*Trans.nr*.pdf", '04Annotation/Transcript_Anno/NR/'),
            ("*Gene*swissprot*.pdf", '04Annotation/Unigene_Anno/Swiss_Prot/'),
            ("Trans*swissprot*.pdf", '04Annotation/Transcript_Anno/Swiss_Prot/'),
            ("Pfam*Gene*.pdf", '04Annotation/Unigene_Anno/Pfam/'),
            ("Pfam*Trans*.pdf", '04Annotation/Transcript_Anno/Pfam/'),
            ("COG*Gene*.pdf", '04Annotation/Unigene_Anno/COG/'),
            ("COG*Trans*.pdf", '04Annotation/Transcript_Anno/COG/'),
            ("Gene_GO*.pdf", '04Annotation/Unigene_Anno/GO/'),
            ("Trans_GO*.pdf", '04Annotation/Transcript_Anno/GO/'),
            ("Gene_Histogram_of_KEGG*.pdf", '04Annotation/Unigene_Anno/KEGG/'),
            ("Trans_Histogram_of_KEGG*.pdf", '04Annotation/Transcript_Anno/KEGG/'),
            ("Unigene_TF*.pdf", '06TF/'),
            ("Transcript_TF*.pdf", '06TF/'),
            ("Genes*exp.distribution*.pdf", '07Express/ExpAnnalysis/'),
            ("Trans*exp.distribution*.pdf", '07Express/ExpAnnalysis/'),
            ("Genes*exp_distribution.violin.pdf", '07Express/ExpAnnalysis/'),
            ("Trans*exp_distribution.violin.pdf", '07Express/ExpAnnalysis/'),
            ("all.exp.venn.pdf", '07Express/ExpVenn/'),
            ("all.exp.heat_corr.pdf", '07Express/ExpCorr/'),
            ("all.exp_relation_pca*.pdf", '07Express/ExpPCA/'),
            ("Gene*diff*.pdf", '08DiffExpress/unigene/'),
            ("Trans*.diff*.pdf", '08DiffExpress/transcript/'),
            ("*length_distribution_of_CDS*.pdf", '09CDS/'),
            ("SSR.stat.pie.pdf", '11SSR/')
        ]
        if self.option('sample_num') == 'multiple':
            file2uploads_geneset = [
                ("geneset.cluster.heat_corr.pdf".format("{geneset_name}"),
                 "12GeneSet/01_Cluster_Analysis/{}".format("{geneset_name}")),
                ("{}.cog_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "12GeneSet/02_COG_Annotation/{}".format("{geneset_name}")),
                ("{}.go_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "12GeneSet/03_GO_Annotation/{}".format("{geneset_name}")),
                ("{}*.kegg_annot.*.column.pdf".format("{geneset_name}"),
                 "12GeneSet/04_KEGG_Annotation/{}".format("{geneset_name}")),
                ("{}.go_enrich.*.pdf".format("{geneset_name}"),
                 "12GeneSet/05_GO_Enrich/{}".format("{geneset_name}")),
                ("{}*kegg_enrich.gene_set*.pdf".format("{geneset_name}"),
                 "12GeneSet/06_KEGG_Enrich/{}".format("{geneset_name}")),
                ("{}.*venn*.pdf".format("diff_genesets"),
                 "12GeneSet/07_GenesetVenn/")

            ]
        for filefrom, fileto in file2uploads:
            pdf_file = glob.glob(self.chart.work_dir + "/" + filefrom)
            for p in pdf_file:
                if os.path.exists(self.target_dir + "/" + fileto + "/" + os.path.basename(p)):
                    os.remove(self.target_dir + "/" + fileto + "/" + os.path.basename(p))
                os.link(p, self.target_dir + "/" + fileto + "/" + os.path.basename(p))

        if self.option('sample_num') == 'multiple':
            genesets = [os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.diff_geneset_analysis.output_dir))]
            #以基因集为单位,专门为差异一键化分析准备
            for geneset in genesets :
                for filefrom, fileto in file2uploads_geneset:
                    pdf_file = glob.glob((self.chart.work_dir + "/" + filefrom).format(geneset_name = geneset))
                    for p in pdf_file:
                        if os.path.exists(self.target_dir + "/" + fileto.format(geneset_name = geneset) + "/" + os.path.basename(p)):
                            os.remove(self.target_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
                        try:
                            os.link(p, self.target_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
                        except:
                            self.logger.info('cuolacuolaraw:{} new:{}'.format(p,
                                                                              self.target_dir + "/" + fileto.format(
                                                                                  geneset_name=geneset) + "/" + os.path.basename(
                                                                                  p)))
        pass

    def def_upload(self):
        if self._sheet.output.endswith("/"):
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        else:
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results/')
        if self.intermediate_result.endswith("/"):
            self.upload_to_s3(self.intermediate_result, intermediate_dir)
        else:
            self.upload_to_s3(self.intermediate_result + "/", intermediate_dir)
        if self.option("report_img"):
            s3 = self._sheet.output.split(":")[0]
            report_img_dir = self.chart.work_dir + '/png/'
            report_img_s3 = s3 + "://commonbucket/files/report_img/drna/" + self.task_id + "/"
            self.upload_to_s3(report_img_dir, report_img_s3)

        # sdir = self.add_upload_dir(self.intermediate_result)
        repaths = [
            [".", "", "流程分析结果目录", 0, "201204"],
            ["Sequence_database", "", "序列文件数据库", 1, "201205"],
            ["Sequence_database/seq_db.sqlite3", "", "序列文件", 1, "201206"],
            ['01Background', '', '项目背景结果目录', 0, "201207"],
            ['01Background/sample_info.xls', 'xls', '样本信息表', 0, "201208"],
            ['01Background/software_info.xls', 'xls', '软件信息表', 0, "201209"],
            ['01Background/run_parameter.txt', 'txt', '分析参数日志', 0],
            ["02QC", "", "测序数据统计与质控结果目录", 0, "201210"],
            ["02QC/rawdata_statistics.xls", "", "原始数据统计表", 0, "201211"],
            ["02QC/cleandata_statistics.xls", "", "质控数据统计表", 0, "201212"],
            # ["QC/cleandata", "", "质控结果目录", 1, '201007'],
            ["03Assemble", "", "从头组装结果目录", 0, "201213"],
            ["03Assemble/01Assemble_Sequence", "", "组装序列文件", 0, "201214"],
            ["03Assemble/01Assemble_Sequence/Trinity.fasta", "", "Transcript（原始组装）序列文件", 0, "201215"],
            ["03Assemble/01Assemble_Sequence/Trinity.gene_trans_map", "", "原始组装结果基因和转录组对应关系表", 0, "201216"],
            ["03Assemble/01Assemble_Sequence/Trinity.filter.fasta", "", "Transcript（优化组装）序列文件", 0, "201217"],
            ["03Assemble/01Assemble_Sequence/Trinity.filter.gene_trans_map", "", "优化组装结果基因和转录组对应关系表", 0, "201218"],
            ["03Assemble/01Assemble_Sequence/Trinity.filter.unigene.fasta", "", "Unigene（优化组装）序列文件", 0, "201219"],
            ["03Assemble/01Assemble_Sequence/Spadas.fasta", "", "Transcript（原始组装）序列文件", 0, "201215"],
            ["03Assemble/01Assemble_Sequence/Spadas.filter.fasta", "", "Transcript（优化组装）序列文件", 0, "201217"],
            ["03Assemble/01Assemble_Sequence/Spadas.filter.unigene.fasta", "", "Unigene（优化组装）序列文件", 0, "201219"],
            ["03Assemble/02Assemble_Evaluation", "", "组装结果评估", 0, "201220"],
            # ["03Assemble/02Assemble_Evaluation/Optimize_before", "", "原始组装结果评估", 0],
            ["03Assemble/02Assemble_Evaluation/Trinity_stat.xls", "", "原始组装结果评估表", 0, "201221"],
            # ["03Assemble/02Assemble_Evaluation/Optimize_after", "", "优化组装结果评估", 0],
            ["03Assemble/02Assemble_Evaluation/Trinity_filter_stat.xls", "", "优化组装结果评估表", 0, "201222"],
            ["03Assemble/03Length_Distribution", "", "转录本长度分布", 0, "201223"],
            ["03Assemble/03Length_Distribution/unigene_count_stat_500.xls", "", "Unigene长度分布统计表", 0, "201224"],
            ["03Assemble/03Length_Distribution/transcript_count_stat_500.xls", "", "Transcript长度分布统计表", 0, "201225"],
            ["03Assemble/Mapping", "", "测序数据与组装结果比对", 0, "201226"],
            ["03Assemble/alignment_rate.xls", "", "测序数据与组装结果比对结果表", 0, "201227"],
            ["04Annotation", "", "功能注释结果目录", 0, "201228"],
            ["04Annotation/all_annotation_statistics.xls", "", "注释概况统计表", 0, "201229"],
            # ["04Annotation/blast_xml/all_annotation_statistics.xls", "", "注释概况统计表", 0, '201027'],
            # ["04Annotation/blast_xml/Trinity_vs_kegg.xml", "", "与KEGG数据库比对xml文件", 0, '201028'],
            # ["04Annotation/blast_xml/blast2go_merge.xls", "", "与GO数据库比对获得注释文件", 0, '201029'],
            # ["04Annotation/blast_xml/Trinity_vs_nr.xml", "", "与NR数据库比对xml文件", 0, '201030'],
            # ["04Annotation/blast_xml/Trinity_vs_string.xml", "", "与COG数据库比对xml文件", 0, '201031'],
            # ["04Annotation/blast_xml/pfam_domain", "", "与Pfam数据库比对获得注释文件", 0, '201032'],
            # ["04Annotation/blast_xml/Trinity_vs_swissprot.xml", "", "与Swiss-Prot数据库比对xml文件", 0, '201033'],
            # ["04Annotation/blast_xml", "", "blast结果文件", 0, '201034'],
            ["04Annotation/Unigene_Anno", "", "Unigene注释结果", 0, "201230"],
            ["04Annotation/Unigene_Anno/GO", "", "GO分类统计结果", 0, "201231"],
            ["04Annotation/Unigene_Anno/GO/unigene_gos.list", "", "序列对应GO分类列表", 0, "201232"],
            ["04Annotation/Unigene_Anno/GO/unigene_go12level_statistics.xls", "", "GO二级分类统计表", 0, "201233"],
            ["04Annotation/Unigene_Anno/GO/unigene_go123level_statistics.xls", "", "GO三级分类统计表", 0, "201234"],
            ["04Annotation/Unigene_Anno/GO/unigene_go1234level_statistics.xls", "", "GO四级分类统计表", 0, "201235"],
            ["04Annotation/Unigene_Anno/NR", "", "NR注释结果", 0, "201236"],
            ["04Annotation/Unigene_Anno/NR/unigene_nr_anno_detail.xls", "", "unigene NR注释详情表", 0, 0, "201237"],
            ["04Annotation/Unigene_Anno/Swiss_Prot", "", "Swiss-Prot注释", 0, "201238"],
            ["04Annotation/Unigene_Anno/Swiss_Prot/unigene_swissprot_anno_detail.xls", "", "Swissprot注释详情表", 0,
             "201239"],
            ["04Annotation/Unigene_Anno/COG", "", "COG注释结果", 0, "201240"],
            ["04Annotation/Unigene_Anno/COG/unigene_cog_summary.xls", "", "COG分类统计表", 0, "201241"],
            ["04Annotation/Unigene_Anno/KEGG", "", "KEGG注释结果", 0, "201242"],
            ["04Annotation/Unigene_Anno/KEGG/unigene_pathway_level.xls", "", "Pathway分类统计表", 0, "201243"],
            ["04Annotation/Unigene_Anno/KEGG/unigene_kegg_table.xls", "", "序列对应的Ko编号统计表", 0, "201244"],
            ["04Annotation/Unigene_Anno/KEGG/unigene_pathways_table.xls", "", "Pathway对应的序列统计表", 0, "201245"],
            ["04Annotation/Unigene_Anno/KEGG/unigene_pathways", "", "KEGG分析结果通路图", 0, "201246"],
            ["04Annotation/Unigene_Anno/Pfam", "", "Pfam注释结果", 0, "201247"],
            ["04Annotation/Unigene_Anno/Pfam/unigene_pfam_anno_detail.xls", "", "Pfam注释详情表", 0, "201248"],
            # ["04Annotation/Unigene_Anno/GO", "", "GO分类统计结果", 0, '201048'],
            # ["04Annotation/Unigene_Anno/GO/unigene_gos.list", "", "序列对应GO分类列表", 0, '201052'],
            # ["04Annotation/Unigene_Anno/GO/unigene_go12level_statistics.xls", "", "GO二级分类统计表", 0, '201052'],
            # ["04Annotation/Unigene_Anno/GO/unigene_go123level_statistics.xls", "", "GO三级分类统计表", 0, '201052'],
            # ["04Annotation/Unigene_Anno/GO/unigene_go1234level_statistics.xls", "", "GO四级分类统计表", 0, '201052'],
            # ["04Annotation/Unigene_Anno/NR", "", "NR注释结果", 0, '201035'],
            # ["04Annotation/Unigene_Anno/NR/unigene_nr_anno_detail.xls", "", "unigene NR注释详情表", 0, '201036'],
            # ["04Annotation/Unigene_Anno/Swiss_Prot", "", "Swiss-Prot注释", 0, '201038'],
            # ["04Annotation/Unigene_Anno/Swiss_Prot/unigene_swissprot_anno_detail.xls", "", "Swissprot注释详情表", 0,
            #  '201038'],
            # ["04Annotation/Unigene_Anno/COG", "", "COG注释结果", 0, '201038'],
            # ["04Annotation/Unigene_Anno/COG/unigene_cog_summary.xls", "", "COG分类统计表", 0, '201038'],
            # ["04Annotation/Unigene_Anno/KEGG", "", "KEGG注释结果", 0, '201038'],
            # ["04Annotation/Unigene_Anno/KEGG/unigene_pathway_level.xls", "", "Pathway分类统计表", 0, '201038'],
            # ["04Annotation/Unigene_Anno/KEGG/unigene_kegg_table.xls", "", "序列对应的Ko编号统计表", 0, '201038'],
            # ["04Annotation/Unigene_Anno/KEGG/unigene_pathways_table.xls", "", "Pathway对应的序列统计表", 0, '201038'],
            # ["04Annotation/Unigene_Anno/KEGG/unigene_pathways", "", "KEGG分析结果通路图", 0, '201038'],
            # ["04Annotation/Unigene_Anno/Pfam", "", "Pfam注释结果", 0, '201038'],
            # ["04Annotation/Unigene_Anno/Pfam/unigene_pfam_anno_detail.xls", "", "Pfam注释详情表", 0, '201038'],
            ["04Annotation/Transcript_Anno", "", "Transcript注释结果", 0, "201249"],
            ["04Annotation/Transcript_Anno/GO", "", "GO分类统计结果", 0, "201250"],
            ["04Annotation/Transcript_Anno/GO/transcript_gos.list", "", "序列对应GO分类列表", 0, "201251"],
            ["04Annotation/Transcript_Anno/GO/transcript_go12level_statistics.xls", "", "GO二级分类统计表", 0, "201252"],
            ["04Annotation/Transcript_Anno/GO/transcript_go123level_statistics.xls", "", "GO三级分类统计表", 0, "201253"],
            ["04Annotation/Transcript_Anno/GO/transcript_go1234level_statistics.xls", "", "GO四级分类统计表", 0, "201254"],
            ["04Annotation/Transcript_Anno/NR", "", "NR注释结果", 0, "201255"],
            ["04Annotation/Transcript_Anno/NR/transcript_nr_anno_detail.xls", "", "unigene NR注释详情表", 0, "201256"],
            ["04Annotation/Transcript_Anno/Swiss_Prot", "", "Swiss-Prot注释", 0, "201257"],
            ["04Annotation/Transcript_Anno/Swiss_Prot/transcript_swissprot_anno_detail.xls", "", "Swissprot注释详情表", 0,
             "201258"],
            ["04Annotation/Transcript_Anno/COG", "", "COG注释结果", 0, "201259"],
            ["04Annotation/Transcript_Anno/COG/transcript_cog_summary.xls", "", "COG分类统计表", 0, "201260"],
            ["04Annotation/Transcript_Anno/KEGG", "", "KEGG注释结果", 0, "201261"],
            ["04Annotation/Transcript_Anno/KEGG/transcript_pathway_level.xls", "", "Pathway分类统计表", 0, "201262"],
            ["04Annotation/Transcript_Anno/KEGG/transcript_kegg_table.xls", "", "序列对应的Ko编号统计表", 0, "201263"],
            ["04Annotation/Transcript_Anno/KEGG/transcript_pathways_table.xls", "", "Pathway对应的序列统计表", 0, "201264"],
            ["04Annotation/Transcript_Anno/KEGG/transcript_pathways", "", "KEGG分析结果通路图", 0, "201265"],
            ["04Annotation/Transcript_Anno/Pfam", "", "Pfam注释结果", 0, "201266"],
            ["04Annotation/Transcript_Anno/Pfam/transcript_pfam_anno_detail.xls", "", "Pfam注释详情表", 0, "201267"],
            ["05AnnoQuery", "", "功能查询结果目录", 0, "201268"],
            ["05AnnoQuery/unigene_anno_detail.xls", "", "Unigene功能查询详情表", 0, "201269"],
            ["05AnnoQuery/transcript_anno_detail.xls", "", "Transcript功能查询详情表", 0, "201270"],
            ["06TF", "", "转录因子分析结果目录", 0, "201271"],
            ["06TF/transcript_tf_predict.xls", "", "Transcript转录因子预测结果表", 0, "201272"],
            ["06TF/unigene_tf_predict.xls", "", "Unigene转录因子预测结果表", 0, "201273"],
            ["07Express", "", "表达量分析结果目录", 0, "201274"],
            ["07Express/ExpAnnalysis", "", "表达量统计", 0, "201275"],
            ["07Express/ExpAnnalysis/unigene.count.matrix.annot.xls", "", "Unigene count表达定量结果表", 0, "201276"],
            ["07Express/ExpAnnalysis/transcript.count.matrix.annot.xls", "", "Transcript_count表达定量注释结果表", 0, "201277"],
            ["07Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls", "", "Unigene tpm表达定量注释结果表", 0, "201278"],
            ["07Express/ExpAnnalysis/unigene.fpkm.matrix.annot.xls", "", "Unigene fpkm表达定量注释结果表", 0, "201279"],
            ["07Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls", "", "Transcript_tmp表达定量注释结果表", 0, "201280"],
            ["07Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls", "", "Transcript_fpkm表达定量注释结果表", 0, "201281"],
            ["07Express/ExpCorr", "", "样本间相关性分析", 0, "201282"],
            ["07Express/ExpCorr/sample_correlation.xls", "", "样本间相关性分析矩阵表", 0, "201283"],
            ["07Express/ExpPCA", "", "样本间PCA分析", 0, "201284"],
            ["07Express/ExpPCA/PCA.xls", "", "样本间PCA分析结果表", 0, "201285"],
            ["07Express/ExpPCA/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表", 0, "201286"],
            ["08DiffExpress", "", "表达量差异分析结果目录", 0, "201287"],
            ["08DiffExpress/unigene", "", "Unigene差异表达分析", 0, "201288"],
            ["08DiffExpress/transcript", "", "Transcript差异表达分析", 0, "201289"],
            ["09CDS", "", "CDS预测结果目录", 0, "201290"],
            # ["09CDS/pfam_domain", "", "转录本序列pfam注释结果", 0, '201091'],
            ["09CDS/all_predicted.bed", "", "CDS序列在组装结果中的相对位置文件", 0, "201291"],
            ["09CDS/all_predicted.cds", "", "CDS的核苷酸序列文件", 0, "201292"],
            ["09CDS/all_predicted.pep", "", "CDS的蛋白序列文件", 0, "201293"],
            ["09CDS/unigene_cds_len.xls", "", "Unigene_cds长度分布", 0, "201294"],
            ["09CDS/transcript_cds_len.xls", "", "Transcript_cds长度分布", 0, "201295"],
            ["09CDS/all_predicted.xls", "", "CDS预测详情表", 0, "201296"],
            ["10SNP", "", "SNP分析结果目录", 0, "201297"],
            # ["10SNP/call.vcf", "", "SNP原始结果文件", 0, '201099'],
            # ["10SNP/filter_snp_vcf", "", "SNP原始结果文件", 0, '201101'],
            # ["10SNP/filter_indel_vcf", "", "SNP原始结果文件", 0, '201102'],
            ["10SNP/snp_detail.xls", "", "SNP结果详情表", 0, "201298"],
            ["10SNP/snp_depth_statistics.xls", "", "SNP测序深度统计表", 0, "201299"],
            ["10SNP/snp_homo_hete_statistics.xls", "", "SNP类型统计表", 0, "201300"],
            ["10SNP/snp_transition_tranversion_statistics.xls", "", "SNP位点统计表", 0, "201301"],
            ["10SNP/snp_cds_statistics.xls", "", "SNP功能区域统计表", 0, "201302"],
            ["10SNP/snp_anno_statistics.xls", "", "SNP功能注释统计表", 0, "201303"],
            ["11SSR", "", "SSR分析结果目录", 0, "201304"],
            ["11SSR/ssr_type.xls", "", "SSR类型统计表", 0, "201305"],
            ["11SSR/ssr_repeats.xls", "", "SSR类型分布统计表", 0, "201306"],
            ["11SSR/ssr_analysis_details.xls", "", "SSR分析详情表", 0, "201307"],
            ['12GeneSet', '', '基因集分析结果目录', 0],
            ['12GeneSet/07_GenesetVenn', '', '差异基因集venn分析结果目录', 0],
            ['12GeneSet/07_GenesetVenn/*.pdf', '', '差异基因集venn图', 0],
            ['12GeneSet/01_Cluster_Analysis', '', '基因集聚类分析结果目录', 0],
            ['12GeneSet/01_Cluster_Analysis/*', '', '基因集聚类分析结果目录', 0],
            ["12GeneSet/01_Cluster_Analysis/*/seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件",0,"211531"],
            ["12GeneSet/01_Cluster_Analysis/*/seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
            ["12GeneSet/01_Cluster_Analysis/*/sample.cluster_tree.txt", "txt", "样本聚类树文件",0,"211532"],
            ["12GeneSet/01_Cluster_Analysis/*/expression_matrix.xls", "xls", "聚类热图分析表",0,"211533"],
            ["12GeneSet/01_Cluster_Analysis/*/seq.subcluster_*.xls", "xls", "子聚类分析表",0],
            ["12GeneSet/01_Cluster_Analysis/*/heatmap.pdf", 'pdf', "聚类热图",0],
            ["12GeneSet/01_Cluster_Analysis/*/*heat_corr.pdf", 'pdf', "聚类热图", 0],
            ["12GeneSet/01_Cluster_Analysis/*/*.line.pdf", 'pdf', "子聚类折线图",0],
            ['12GeneSet/02_COG_Annotation', '', '基因集聚COG分类结果目录', 0],
            ['12GeneSet/02_COG_Annotation/*', '', '基因集聚COG分类结果目录', 0],
            ["12GeneSet/02_COG_Annotation/*/cog_class_stat.xls", "", "COG分类统计表",0,"211529"],
            ["12GeneSet/02_COG_Annotation/*/*.pdf", "pdf", "COG分类统计图",0],
            ['12GeneSet/03_GO_Annotation', '', '基因集GO分类结果目录', 0],
            ['12GeneSet/03_GO_Annotation/*', '', '基因集GO分类结果目录', 0],
            ["12GeneSet/03_GO_Annotation/*/go_class_stat.xls", "", "GO分类统计表",0,"211527"],
            ["12GeneSet/03_GO_Annotation/*/*.pdf", "pdf", "GO分类统计图",0],
            ['12GeneSet/04_KEGG_Annotation', '', '基因集KEGG分类结果目录', 0],
            ['12GeneSet/04_KEGG_Annotation/*', '', '基因集KEGG分类结果目录', 0],
            ["12GeneSet/04_KEGG_Annotation/*/pathways.tar.gz", "", "KEGG通路图",0,"211555"],

            ["12GeneSet/04_KEGG_Annotation/*/kegg_class_stat.xls", "", "Pathway分类统计表",0,"211557"],
            ["12GeneSet/04_KEGG_Annotation/*/*.pdf", "pdf", "KEGG分类统计图",0],
            ['12GeneSet/05_GO_Enrich', '', '基因集GO富集分析结果目录', 0],
            ['12GeneSet/05_GO_Enrich/*', '', '基因集GO富集分析结果目录', 0],
            ['12GeneSet/05_GO_Enrich/*/go_enrich_stat.xls', 'xls', 'GO富集分析统计表',0,"211536"],
            ["12GeneSet/05_GO_Enrich/*/*bar.pdf", "pdf", "GO富集分析柱形图",0],
            ["12GeneSet/05_GO_Enrich/*/*bar_line.pdf", "pdf", "GO富集分析柱形图(带折线)", 0],
            ["12GeneSet/05_GO_Enrich/*/*buble1.pdf", "pdf", "GO富集分析气泡图(单基因集)", 0],
            ["12GeneSet/05_GO_Enrich/*/*buble2.pdf", "pdf", "GO富集分析气泡图(分散型)", 0],
            ['12GeneSet/06_KEGG_Enrich', '', '基因集KEGG富集分析结果目录', 0],
            ['12GeneSet/06_KEGG_Enrich/*', '', '基因集KEGG富集分析结果目录', 0],
            ['12GeneSet/06_KEGG_Enrich/*/kegg_enrich_stat.xls', '', 'KEGG富集分析统计表',0,"211538"],
            ['12GeneSet/06_KEGG_Enrich/*/pathways.tar.gz', ' ', 'KEGG富集通路图',0,"211539"],
            ["12GeneSet/06_KEGG_Enrich/*/*bar.pdf", "pdf", "KEGG富集分析柱形图", 0],
            ["12GeneSet/06_KEGG_Enrich/*/*bar_line.pdf", "pdf", "KEGG富集分析柱形图(带折线)", 0],
            ["12GeneSet/06_KEGG_Enrich/*/*buble1.pdf", "pdf", "KEGG富集分析气泡图(单基因集)", 0],
            ["12GeneSet/06_KEGG_Enrich/*/*buble2.pdf", "pdf", "KEGG富集分析气泡图(分散型)", 0],
        ]
        sdir = self.add_upload_dir(self.target_dir)
        sdir.add_regexp_rules([
            [r'04Annotation/Transcript_Anno/KEGG/transcript_pathways/.*\.png', 'png', 'KEGG通路png图片', 0, "201308"],
            [r'04Annotation/Transcript_Anno/KEGG/transcript_pathways/.*\.html', 'html', 'KEGG通路html图片', 0, "201309"],
            [r'04Annotation/Unigene_Anno/KEGG/unigene_pathways/.*\.png', 'png', '差异表达unigene详情表', 0, "201310"],
            [r'04Annotation/Unigene_Anno/KEGG/unigene_pathways/.*\.html', 'html', '差异表达unigene详情表', 0, "201311"],
            [r'08DiffExpress/unigene/.*_vs_.*\.xls', 'xls', '差异表达unigene详情表', 0, "201312"],
            [r'08DiffExpress/unigene/.*\.DE\.list', 'xls', '差异基因列表', 0, "201313"],
            [r'08DiffExpress/unigene/.*summary.*\.xls', 'xls', '差异表达unigene基因统计表', 0, "201314"],
            [r'08DiffExpress/unigene/.*summary.*\.annot\.xls', 'xls', '差异表达unigene统计详情表', 0, "201315"],
            [r'08DiffExpress/unigene/.*total_diff_stat.*\.xls', 'xls', '差异表达unigene详情总表', 0, "201316"],
            [r'08DiffExpress/transcript/.*_vs_.*\.xls', 'xls', '差异表达transcript详情表', 0, "201317"],
            [r'08DiffExpress/transcript/.*\.DE\.list', 'xls', '差异基因列表', 0, "201318"],
            [r'08DiffExpress/transcript/.*summary.*\.xls', 'xls', '差异表达transcript统计表', 0, "201319"],
            [r'08DiffExpress/transcript/.*summary.*\.annot\.xls', 'xls', '差异表达transcript统计详情表', 0, "201320"],
            [r'08DiffExpress/transcript/.*total_diff_stat.*\.xls', 'xls', '差异表达transcript详情总表', 0, "201321"],
            ## 图片描述
            [r"02QC/.*raw_qc_qual\.box\.pdf", 'pdf', '原始数据碱基质量分布图', 0],
            [r"02QC/.*raw_qc_error\.line\.pdf", 'pdf', '原始数据碱基错误率分布图', 0],
            [r"02QC/.*raw_qc_base\.line\.pdf", 'pdf', '原始数据碱基含量分布图', 0],
            [r"02QC/.*clean_qc_qual\.box\.pdf", 'pdf', '质控数据碱基质量分布图', 0],
            [r"02QC/.*clean_qc_error\.line\.pdf", 'pdf', '质控数据碱基错误率分布图', 0],
            [r"02QC/.*clean_qc_base\.line\.pdf", 'pdf', '质控数据碱基含量分布图', 0],
            [r"03Assemble/.*length_distribution.*\.pdf", 'pdf', '转录本长度分布图', 0],
            [r"04Annotation/Unigene_Anno/.*annot_gene_stat.*\.pdf", 'pdf', '注释结果图', 0],
            [r"04Annotation/Transcript_Anno/.*annot_tran_stat.*\.pdf", 'pdf', '注释结果图', 0],
            [r"04Annotation/Unigene_Anno/NR/.*nr_species\.pie\.pdf", 'pdf', 'NR注释物种分布饼图', 0],
            [r"04Annotation/Unigene_Anno/NR/.*nr_evalue\.pie.\pdf", 'pdf', 'E-value分布饼图', 0],
            [r"04Annotation/Unigene_Anno/NR/.*nr_similary\.pie.\pdf", 'pdf', '相似度分布饼图', 0],
            [r"04Annotation/Transcript_Anno/NR/.*nr_species\.pie\.pdf", 'pdf', 'NR注释物种分布饼图', 0],
            [r"04Annotation/Transcript_Anno/NR/.*nr_evalue\.pie.\pdf", 'pdf', 'E-value分布饼图', 0],
            [r"04Annotation/Transcript_Anno/NR/.*nr_similary\.pie.\pdf", 'pdf', '相似度分布饼图', 0],
            [r"04Annotation/Unigene_Anno/Swiss_Prot/.*nr_evalue\.pie.\pdf", 'pdf', 'E-value分布饼图', 0],
            [r"04Annotation/Unigene_Anno/Swiss_Prot/.*nr_similary\.pie.\pdf", 'pdf', '相似度分布饼图', 0],
            [r"04Annotation/Transcript_Anno/Swiss_Prot/.*nr_evalue\.pie.\pdf", 'pdf', 'E-value分布饼图', 0],
            [r"04Annotation/Transcript_Anno/Swiss_Prot/.*nr_similary\.pie.\pdf", 'pdf', '相似度分布饼图', 0],
            [r"04Annotation/Unigene_Anno/Pfam/Pfam.*Gene.*\.pdf", 'pdf', 'Pfam注释柱状图', 0],
            [r"04Annotation/Transcript_Anno/Pfam/Pfam.*Trans*\.pdf", 'pdf', 'Pfam注释柱状图', 0],
            [r"04Annotation/Unigene_Anno/COG/.*COG.*Gene.*\.pdf", 'pdf', 'COG分类统计柱状图', 0],
            [r"04Annotation/Transcript_Anno/COG/.*COG.*Trans.*\.pdf", 'pdf', 'COG分类统计柱状图', 0],
            [r"04Annotation/Unigene_Anno/GO/.*Gene_GO.*\.pdf", 'pdf', 'GO分类统计图', 0],
            [r"04Annotation/Transcript_Anno/GO/.*Trans_GO.*\.pdf", 'pdf', 'GO分类统计图', 0],
            [r"04Annotation/Unigene_Anno/KEGG/.*Gene_Histogram_of_KEGG.*\.pdf", 'pdf', 'Pathway分类统计柱状图', 0],
            [r"04Annotation/Transcript_Anno/KEGG/.*Trans_Histogram_of_KEGG.*\.pdf", 'pdf', 'Pathway分类统计柱状图', 0],
            [r"06TF/.*Unigene_TF.*\.pdf", 'pdf', '转录因子家族统计图', 0],
            [r"06TF/.*Transcript_TF.*\.pdf", 'pdf', '转录因子家族统计图', 0],
            [r"07Express/ExpAnnalysis/.*exp.*distribution.*\.pdf", 'pdf', '表达量分布图', 0],
            [r"07Express/ExpVenn/all.exp.venn.pdf", 'pdf', '样本间Venn图', 0],
            [r"07Express/ExpCorr/all.exp.heat_corr.pdf", 'pdf', '样本间相关性热图', 0],
            [r"07Express/ExpPCA/.*all.exp_relation_pca.*\.pdf", 'pdf', 'PCA图', 0],
            [r"08DiffExpress/unigene/Gene\.differential\.summary\.bar.*pdf", 'pdf', '表达量差异统计图', 0],
            [r"08DiffExpress/transcript/Trans\.differential\.summary\.bar.*pdf", 'pdf', '表达量差异统计图', 0],
            [r"08DiffExpress/unigene/.*diff\.volcano.*pdf", 'pdf', '火山图', 0],
            [r"08DiffExpress/transcript/.*diff\.volcano.*pdf", 'pdf', '火山图', 0],
            [r"08DiffExpress/unigene/.*diff\.scatter.*pdf", 'pdf', '散点图', 0],
            [r"08DiffExpress/transcript/.*diff\.scatter.*pdf", 'pdf', '散点图', 0],
            [r"09CDS/.*length_distribution_of_CDS.*\.pdf", 'pdf', 'CDS长度分布图', 0],
            [r"11SSR/SSR.stat.pie.pdf", 'pdf', 'SSR统计图'],

        ])
        sdir.add_relpath_rules(repaths)
        self.update_collections()

    def update_collections(self):
        if self._sheet.output.endswith("/"):
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        else:
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results/')
        db = Config().get_mongo_client(mtype='denovo_rna_v2')[Config().get_mongo_dbname('denovo_rna_v2')]
        col0 = db['sg_task']
        col0.update({'task_id': self.task_id},
                    {'$set': {'seq_db': os.path.join(intermediate_dir, 'SequenceDatabase/seq_db.sqlite3')}},
                    upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'seqdetail': os.path.join(intermediate_dir, 'SequenceDetail/')}},
                    upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {"fastq": os.path.join(intermediate_dir, "QC/fq_list.txt")}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'assemble_fa': os.path.join(intermediate_dir, 'Assemble',
                                                                                     os.path.basename(
                                                                                         self.assemble_filter.option(
                                                                                             "filter_fa").prop[
                                                                                             "path"]))}}, upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'assemble_t2g': os.path.join(intermediate_dir, 'Assemble', "all_tran2gen.txt")}},
                    upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'bedpath': os.path.join(intermediate_dir, 'CDS', "all_predicted.bed")}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'unigene_fa': os.path.join(intermediate_dir, 'Assemble',
                                                                                    os.path.basename(
                                                                                        self.assemble_filter.option(
                                                                                            "unigene_filter_fa").prop[
                                                                                            'path']))}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'fq_type': self.option('fq_type')}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'strand_specific': self.option('strand_specific')}},
                    upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'strand_dir': self.option('strand_dir')}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'version': 'v2'}}, upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'assembly_dir': os.path.join(intermediate_dir, 'Assembly_dir')}}, upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'assembly_object': {"raw": "assemble_raw.fasta", "filter": "filter.fasta",
                                                  "unigene": "filter.unigene.fasta", "peptide": "all_predicted.pep",
                                                  "map": "all_tran2gen.txt"}}}, upsert=True)
        annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        if self.annot_config_dict['kegg']['version'] > "2020":
            pass
            # if self.option("kegg_org") not in [None, ""]:
            #     annot_version_dict['kegg'] += "_spe"
        else:
            del annot_version_dict['kegg']
        col0.update({'task_id': self.task_id},
                    {'$set': {'database_version': annot_version_dict,
                              'annot_group': self.option("annot_group")}}, upsert=True)
        col1 = db['sg_annotation_stat']
        col1.update({'task_id': self.task_id},
                    {'$set': {'result_dir': os.path.join(intermediate_dir, 'Annotation')}}, upsert=True)
        col2 = db['sg_exp']
        col2.update({'task_id': self.task_id, 'exp_level': 'G'}, {
            '$set': {'count_file': os.path.join(intermediate_dir,
                                                'Express/ExpAnnalysis/gene.count.matrix.xls')}}, upsert=True)
        col2.update({'task_id': self.task_id, 'exp_level': 'T'}, {
            '$set': {'count_file': os.path.join(intermediate_dir,
                                                'Express/ExpAnnalysis/transcript.count.matrix.xls')}}, upsert=True)
        super(DenovornaWorkflow, self).end()

    def modify_output(self):
        # 文件上传到页面后，需要修改文件的内容，才能用于交互分析
        # old_name = self.work_dir + "/HiseqQc/output/sickle_dir/fq_list.txt"
        old_name = self.qc.option('fq_list').prop['path']
        a = open(old_name, "r").read()
        start_dir = self.work_dir + "/HiseqQc/output/sickle_dir/"
        end_dir = self.workflow_output + "/QC/cleandata/"
        b = a.replace(start_dir, end_dir)
        new_name = self.work_dir + "/HiseqQc/output/sickle_dir/tmp.txt"
        if self.option('fastp'):
            new_name = self.work_dir + "/FastpRna/output/fastq/tmp.txt"
        with open(new_name, "w") as f:
            f.write(b)
        os.remove(old_name)
        os.rename(new_name, old_name)

        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"
        # Sequence_database
        os.mkdir(target_dir + "/Sequence_database")
        seq_db = self.work_dir + "/seq_db.sqlite3"
        os.link(seq_db, target_dir + "/Sequence_database/seq_db.sqlite3")
        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        fq_list = origin_dir + "/QC_stat/sickle_dir/fq_list.txt"
        if self.option('fastp'):
            fq_list = origin_dir + "/QC_stat/fastq/fq_list.txt"
        os.mkdir(target_dir + "/QC")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
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
        # Assemble
        os.mkdir(target_dir + "/Assemble")
        os.mkdir(target_dir + "/Assemble/AssembleResult")
        if self.option("assemble") != False:
            orgin_trinity = os.path.join(self.assemble.output_dir + "/Trinity.fasta")
            if os.path.exists(orgin_trinity):
                os.link(orgin_trinity, target_dir + "/Assemble/AssembleResult/Trinity.fasta")
            orgin_trinity_g2t = os.path.join(self.assemble.output_dir + "/Trinity.fasta.gene_trans_map")
            if os.path.exists(orgin_trinity_g2t):
                os.link(orgin_trinity_g2t, target_dir + "/Assemble/AssembleResult/Trinity.gene_trans_map")
        orgin_trinity_filter = os.path.join(origin_dir + "/assemble/Trinity.filter.fasta")
        os.link(orgin_trinity_filter, target_dir + "/Assemble/AssembleResult/Trinity.filter.fasta")
        orgin_trinity_filter_g2t = os.path.join(origin_dir + "/assemble/Trinity.filter.gene_trans_map")
        os.link(orgin_trinity_filter_g2t, target_dir + "/Assemble/AssembleResult/Trinity.filter.gene_trans_map")
        orgin_unigene_filter = os.path.join(origin_dir + "/assemble/Trinity.filter.unigene.fasta")
        os.link(orgin_unigene_filter, target_dir + "/Assemble/AssembleResult/Trinity.filter.unigene.fasta")
        os.mkdir(target_dir + "/Assemble/AssembleEvaluation")
        os.mkdir(target_dir + "/Assemble/AssembleEvaluation/Optimize_before")
        optimize_before_stat = os.path.join(origin_dir + "/assemble/trinity_stat/Trinity_stat.xls")
        os.link(optimize_before_stat, target_dir + "/Assemble/AssembleEvaluation/Optimize_before/Trinity_stat.xls")
        if self.option("optimize") == True:
            os.mkdir(target_dir + "/Assemble/AssembleEvaluation/Optimize_after")
            optimize_after_stat = os.path.join(origin_dir + "/assemble_evaluation/trinity_stat/Trinity_stat.xls")
            os.link(optimize_after_stat, target_dir + "/Assemble/AssembleEvaluation/Optimize_after/Trinity_stat.xls")
        os.mkdir(target_dir + "/Assemble/LengthDistribution")
        if self.option("optimize") == False:
            gene_len_dis = os.path.join(origin_dir + "/assemble/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, target_dir + "/Assemble/LengthDistribution/unigene_count_stat_500.txt")
            trans_len_dis = os.path.join(origin_dir + "/assemble/trinity_stat/trans_count_stat_500.txt")
            os.link(trans_len_dis, target_dir + "/Assemble/LengthDistribution/transcript_count_stat_500.txt")
        else:
            gene_len_dis = os.path.join(origin_dir + "/assemble_evaluation/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, target_dir + "/Assemble/LengthDistribution/unigene_count_stat_500.txt")
            trans_len_dis = os.path.join(origin_dir + "/assemble_evaluation/trinity_stat/trans_count_stat_500.txt")
            os.link(trans_len_dis, target_dir + "/Assemble/LengthDistribution/transcript_count_stat_500.txt")
        os.mkdir(target_dir + "/Assemble/Mapping")
        align_result = os.path.join(self.align.output_dir + "/alignment_rate.txt")
        os.link(align_result, target_dir + "/Assemble/Mapping/alignment_rate.txt")
        # Annotation
        os.mkdir(target_dir + "/Annotation")
        os.mkdir(target_dir + "/Annotation/blast_xml")
        blast_nr = self.diamond_nr.option('outxml').prop['path']
        os.link(blast_nr, target_dir + "/Annotation/blast_xml/Trinity_vs_nr.xml")
        blast_swiss = self.blast_swissprot.option('outxml').prop['path']
        os.link(blast_swiss, target_dir + "/Annotation/blast_xml/Trinity_vs_swissprot.xml")
        blast_string = self.diamond_string.option('outxml').prop['path']
        os.link(blast_string, target_dir + "/Annotation/blast_xml/Trinity_vs_string.xml")
        blast_kegg = self.diamond_kegg.option('outxml').prop['path']
        os.link(blast_kegg, target_dir + "/Annotation/blast_xml/Trinity_vs_kegg.xml")
        blast_pfam = self.cds_predict.output_dir + "/pfam_domain"
        os.link(blast_pfam, target_dir + "/Annotation/blast_xml/pfam_domain")
        os.mkdir(target_dir + "/Annotation/AnnoStat")
        gene_anno_stat = os.path.join(origin_dir + "/annotation/anno_stat/all_annotation_statistics.xls")
        os.link(gene_anno_stat, target_dir + "/Annotation/AnnoStat/all_annotation_statistics.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/NR")
        gene_nr = os.path.join(origin_dir + "/annotation/anno_stat/blast/gene_nr.xls")
        os.link(gene_nr, target_dir + "/Annotation/AnnoStat/NR/unigene_nr_anno_detail.xls")
        trans_nr = os.path.join(origin_dir + "/annotation/anno_stat/blast/nr.xls")
        os.link(trans_nr, target_dir + "/Annotation/AnnoStat/NR/transcript_nr_anno_detail.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/Swiss-Prot")
        gene_swiss = os.path.join(origin_dir + "/annotation/anno_stat/blast/gene_swissprot.xls")
        os.link(gene_swiss, target_dir + "/Annotation/AnnoStat/Swiss-Prot/unigene_swissprot_anno_detail.xls")
        trans_swiss = os.path.join(origin_dir + "/annotation/anno_stat/blast/swissprot.xls")
        os.link(trans_swiss, target_dir + "/Annotation/AnnoStat/Swiss-Prot/transcript_swissprot_anno_detail.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/Pfam")
        gene_pfam = os.path.join(origin_dir + "/annotation/anno_stat/pfam_stat/gene_pfam_domain")
        os.link(gene_pfam, target_dir + "/Annotation/AnnoStat/Pfam/unigene_pfam_anno_detail.xls")
        trans_pfam = os.path.join(origin_dir + "/annotation/blast_xml/pfam_domain")
        os.link(trans_pfam, target_dir + "/Annotation/AnnoStat/Pfam/transcript_pfam_anno_detail.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/COG")
        gene_coglist = os.path.join(origin_dir + "/annotation/anno_stat/cog_stat/gene_cog_list.xls")
        os.link(gene_coglist, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_list.xls")
        gene_coganno = os.path.join(origin_dir + "/annotation/anno_stat/cog_stat/gene_cog_table.xls")
        os.link(gene_coganno, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_anno_detail.xls")
        gene_cogsum = os.path.join(origin_dir + "/annotation/anno_stat/cog_stat/gene_cog_summary.xls")
        os.link(gene_cogsum, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_summary.xls")
        trans_coglist = os.path.join(origin_dir + "/annotation/cog/cog_list.xls")
        os.link(trans_coglist, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_list.xls")
        trans_coganno = os.path.join(origin_dir + "/annotation/cog/cog_table.xls")
        os.link(trans_coganno, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_anno_detail.xls")
        trans_cogsum = os.path.join(origin_dir + "/annotation/cog/cog_summary.xls")
        os.link(trans_cogsum, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_summary.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/GO")
        gene_golist = os.path.join(origin_dir + "/annotation/anno_stat/go_stat/gene_gos.list")
        os.link(gene_golist, target_dir + "/Annotation/AnnoStat/GO/unigene_gos.list")
        gene_goanno = os.path.join(origin_dir + "/annotation/anno_stat/go_stat/gene_blast2go.annot")
        os.link(gene_goanno, target_dir + "/Annotation/AnnoStat/GO/unigene_blast2go.annot")
        gene_golevel12 = os.path.join(origin_dir + "/annotation/anno_stat/go_stat/gene_go12level_statistics.xls")
        os.link(gene_golevel12, target_dir + "/Annotation/AnnoStat/GO/unigene_go12level_statistics.xls")
        gene_golevel123 = os.path.join(origin_dir + "/annotation/anno_stat/go_stat/gene_go123level_statistics.xls")
        os.link(gene_golevel123, target_dir + "/Annotation/AnnoStat/GO/unigene_go123level_statistics.xls")
        gene_golevel1234 = os.path.join(origin_dir + "/annotation/anno_stat/go_stat/gene_go1234level_statistics.xls")
        os.link(gene_golevel1234, target_dir + "/Annotation/AnnoStat/GO/unigene_go1234level_statistics.xls")
        trans_golist = os.path.join(origin_dir + "/annotation/go/query_gos.list")
        os.link(trans_golist, target_dir + "/Annotation/AnnoStat/GO/transcript_gos.list")
        trans_goanno = os.path.join(origin_dir + "/annotation/go/blast2go.annot")
        os.link(trans_goanno, target_dir + "/Annotation/AnnoStat/GO/transcript_blast2go.annot")
        trans_golevel12 = os.path.join(origin_dir + "/annotation/go/go12level_statistics.xls")
        os.link(trans_golevel12, target_dir + "/Annotation/AnnoStat/GO/transcript_go12level_statistics.xls")
        trans_golevel123 = os.path.join(origin_dir + "/annotation/go/go123level_statistics.xls")
        os.link(trans_golevel123, target_dir + "/Annotation/AnnoStat/GO/transcript_go123level_statistics.xls")
        trans_golevel1234 = os.path.join(origin_dir + "/annotation/go/go1234level_statistics.xls")
        os.link(trans_golevel1234, target_dir + "/Annotation/AnnoStat/GO/transcript_go1234level_statistics.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathways")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathways")
        gene_kegg_table = os.path.join(origin_dir + "/annotation/anno_stat/kegg_stat/gene_kegg_table.xls")
        os.link(gene_kegg_table, target_dir + "/Annotation/AnnoStat/KEGG/unigene_kegg_table.xls")
        gene_pathway_table = os.path.join(origin_dir + "/annotation/anno_stat/kegg_stat/gene_pathway_table.xls")
        os.link(gene_pathway_table, target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathway_table.xls")
        gene_kegg_layer = os.path.join(origin_dir + "/annotation/anno_stat/kegg_stat/gene_kegg_layer.xls")
        os.link(gene_kegg_layer, target_dir + "/Annotation/AnnoStat/KEGG/unigene_kegg_layer.xls")
        trans_kegg_table = os.path.join(origin_dir + "/annotation/kegg/kegg_table.xls")
        os.link(trans_kegg_table, target_dir + "/Annotation/AnnoStat/KEGG/transcript_kegg_table.xls")
        trans_pathway_table = os.path.join(origin_dir + "/annotation/kegg/pathway_table.xls")
        os.link(trans_pathway_table, target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathway_table.xls")
        trans_kegg_layer = os.path.join(origin_dir + "/annotation/kegg/kegg_layer.xls")
        os.link(trans_kegg_layer, target_dir + "/Annotation/AnnoStat/KEGG/transcript_kegg_layer.xls")
        for file in os.listdir(origin_dir + "/annotation/anno_stat/kegg_stat/gene_pathway"):
            gene_pathway_path = os.path.join(origin_dir + "/annotation/anno_stat/kegg_stat/gene_pathway", file)
            os.link(gene_pathway_path, target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathways/" + file)
        for file in os.listdir(origin_dir + "/annotation/kegg/pathways"):
            trans_pathway_path = os.path.join(origin_dir + "/annotation/kegg/pathways", file)
            os.link(trans_pathway_path, target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathways/" + file)
        # AnnoQuery
        os.mkdir(target_dir + "/AnnoQuery")
        gene_anno_detail = os.path.join(origin_dir + "/annotation/anno_stat/gene_anno_detail.xls")
        os.link(gene_anno_detail, target_dir + "/AnnoQuery/unigene_anno_detail.xls")
        trans_anno_detail = os.path.join(origin_dir + "/annotation/anno_stat/trans_anno_detail.xls")
        os.link(trans_anno_detail, target_dir + "/AnnoQuery/transcript_anno_detail.xls")
        # Express
        os.mkdir(target_dir + "/Express")
        os.mkdir(target_dir + "/Express/ExpAnnalysis")
        if self.option("express_method") == "RSEM":
            gene_count = os.path.join(self.align.output_dir + "/gene.count.matrix")
            os.link(gene_count, target_dir + "/Express/ExpAnnalysis/unigene.count.matrix.xls")
            gene_tpm = os.path.join(self.align.output_dir + "/gene.tpm.matrix")
            os.link(gene_tpm, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.xls")
            gene_fpkm = os.path.join(self.align.output_dir + "/gene.fpkm.matrix")
            os.link(gene_fpkm, target_dir + "/Express/ExpAnnalysis/unigene.fpkm.matrix.xls")
            trans_count = os.path.join(self.align.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
            trans_tpm = os.path.join(self.align.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
            trans_fpkm = os.path.join(self.align.output_dir + "/transcript.fpkm.matrix")
            os.link(trans_fpkm, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.xls")
            gene_tpm_anno = os.path.join(self.align.output_dir + "/gene.tpm.matrix.annot.xls")
            os.link(gene_tpm_anno, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls")
            gene_fpkm_anno = os.path.join(self.align.output_dir + "/gene.fpkm.matrix.annot.xls")
            os.link(gene_fpkm_anno, target_dir + "/Express/ExpAnnalysis/unigene.fpkm.matrix.annot.xls")
            trans_tpm_anno = os.path.join(self.align.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
            trans_fpkm_anno = os.path.join(self.align.output_dir + "/transcript.fpkm.matrix.annot.xls")
            os.link(trans_fpkm_anno, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls")
        else:
            gene_count = os.path.join(self.express.output_dir + "/gene.count.matrix")
            os.link(gene_count, target_dir + "/Express/ExpAnnalysis/unigene.count.matrix.xls")
            gene_tpm = os.path.join(self.express.output_dir + "/gene.tpm.matrix")
            os.link(gene_tpm, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.xls")
            trans_count = os.path.join(self.express.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
            trans_tpm = os.path.join(self.express.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
            gene_tpm_anno = os.path.join(self.express.output_dir + "/gene.tpm.matrix.annot.xls")
            os.link(gene_tpm_anno, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls")
            trans_tpm_anno = os.path.join(self.express.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
        if self.option("sample_num") == "multiple":
            os.mkdir(target_dir + "/Express/ExpCorr")
            gene_corr = os.path.join(self.exp_corr.output_dir + "/sample_correlation.xls")
            os.link(gene_corr, target_dir + "/Express/ExpCorr/sample_correlation.xls ")
            if self.option("group").prop["sample_number"] > 2:
                os.mkdir(target_dir + "/Express/ExpPCA")
                gene_pca = os.path.join(self.exp_pca.output_dir + "/PCA.xls")
                os.link(gene_pca, target_dir + "/Express/ExpPCA/PCA.xls")
                gene_ratio = os.path.join(self.exp_pca.output_dir + "/Explained_variance_ratio.xls")
                os.link(gene_ratio, target_dir + "/Express/ExpPCA/Explained_variance_ratio.xls")
            # DiffExpress
            os.mkdir(target_dir + "/DiffExpress")
            os.mkdir(target_dir + "/DiffExpress/DiffExpress_unigene")
            for file in os.listdir(self.diffexpress_gene.output_dir):
                diffexpress_unigene = os.path.join(self.diffexpress_gene.output_dir, file)
                os.link(diffexpress_unigene, target_dir + "/DiffExpress/DiffExpress_unigene/" + file)
            os.mkdir(target_dir + "/DiffExpress/DiffExpress_transcript")
            for file in os.listdir(self.diffexpress_trans.output_dir):
                diffexpress_unigene = os.path.join(self.diffexpress_trans.output_dir, file)
                os.link(diffexpress_unigene, target_dir + "/DiffExpress/DiffExpress_transcript/" + file)
            # SNP
            if self.option("is_snp") == "True":
                os.mkdir(target_dir + "/SNP")
                call = os.path.join(self.snp.output_dir + "/pop.variant.vcf")
                if os.path.exists(call):
                    os.link(call, target_dir + "/SNP/call.vcf")
                dp = os.path.join(self.snp.output_dir + "/snp_depth_statistics")
                if os.path.exists(dp):
                    os.link(dp, target_dir + "/SNP/snp_depth_statistics.xls")
                hh = os.path.join(self.snp.output_dir + "/snp_homo_hete_statistics")
                if os.path.exists(hh):
                    os.link(hh, target_dir + "/SNP/snp_homo_hete_statistics.xls")
                tt = os.path.join(self.snp.output_dir + "/snp_transition_tranversion_statistics")
                if os.path.exists(tt):
                    os.link(tt, target_dir + "/SNP/snp_transition_tranversion_statistics.xls")
                detail = os.path.join(self.snp.output_dir + "/snp_detail")
                if os.path.exists(detail):
                    os.link(detail, target_dir + "/SNP/snp_detail.xls")
        # SSR
        os.mkdir(target_dir + "/SSR")
        type = os.path.join(self.ssr.output_dir + "/ssr_type")
        os.link(type, target_dir + "/SSR/ssr_type.xls")
        rep = os.path.join(self.ssr.output_dir + "/ssr_repeats")
        os.link(rep, target_dir + "/SSR/ssr_repeats.xls")
        ssr = os.path.join(self.ssr.output_dir + "/ssr_analysis_details")
        os.link(ssr, target_dir + "/SSR/ssr_analysis_details.xls")
        # CDS
        os.mkdir(target_dir + "/CDS")
        pfam_domain = os.path.join(self.cds_predict.output_dir + "/pfam_domain")
        os.link(pfam_domain, target_dir + "/CDS/pfam_domain")
        bed = os.path.join(self.cds_predict.output_dir + "/Trinity.filter.fasta.transdecoder.bed")
        os.link(bed, target_dir + "/CDS/Trinity.filter.fasta.transdecoder.bed")
        cds = os.path.join(self.cds_predict.output_dir + "/Trinity.filter.fasta.transdecoder.cds")
        os.link(cds, target_dir + "/CDS/Trinity.filter.fasta.transdecoder.cds")
        pep = os.path.join(self.cds_predict.output_dir + "/Trinity.filter.fasta.transdecoder.pep")
        os.link(pep, target_dir + "/CDS/Trinity.filter.fasta.transdecoder.pep")
        gene_len = os.path.join(self.cds_predict.output_dir + "/cds_len_unigene.txt")
        os.link(gene_len, target_dir + "/CDS/unigene_cds_len.xls")
        trans_len = os.path.join(self.cds_predict.output_dir + "/cds_len_transcript.txt")
        os.link(trans_len, target_dir + "/CDS/transcript_cds_len.xls")
        # TF
        if self.option("tf_database") == "Animal":
            os.mkdir(target_dir + "/TF")
            gene_tf = os.path.join(self.cds_predict.output_dir + "/merge_only_unigene_animal")
            os.link(gene_tf, target_dir + "/TF/unigene_tf_predict.xls")
            trans_tf = os.path.join(self.cds_predict.output_dir + "/merge_only_transcript_animal")
            os.link(trans_tf, target_dir + "/TF/transcript_tf_predict.xls")
        if self.option("tf_database") == "Plant":
            os.mkdir(target_dir + "/TF")
            gene_tf = os.path.join(self.cds_predict.output_dir + "/merge_only_unigene_plant")
            os.link(gene_tf, target_dir + "/TF/unigene_tf_predict.xls")
            trans_tf = os.path.join(self.cds_predict.output_dir + "/merge_only_transcript_plant")
            os.link(trans_tf, target_dir + "/TF/transcript_tf_predict.xls")
        repaths = [
            [".", "", "流程分析结果目录", 0, '201001'],
            ["Sequence_database", "", "序列文件数据库", 1, '201002'],
            ["Sequence_database/seq_db.sqlite3", "", "序列文件", 1, '201003'],
            ["QC", "", "测序数据统计与质控结果", 0, '201004'],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表", 0, '201005'],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表", 0, '201006'],
            ["QC/cleandata", "", "质控结果目录", 1, '201007'],
            ["Assemble", "", "从头组装", 0, '201008'],
            ["Assemble/AssembleResult", "", "组装结果", 0, '201009'],
            ["Assemble/AssembleResult/Trinity.fasta", "", "原始组装结果文件", 0, '201010'],
            ["Assemble/AssembleResult/Trinity.gene_trans_map", "", "原始组装结果基因和转录组对应关系表", 0, '201011'],
            ["Assemble/AssembleResult/Trinity.filter.fasta", "", "优化组装结果文件", 0, '201012'],
            ["Assemble/AssembleResult/Trinity.filter.gene_trans_map", "", "优化组装结果基因和转录组对应关系表", 0, '201013'],
            ["Assemble/AssembleResult/Trinity.filter.unigene.fasta", "", "优化组装unigene序列文件", 0, '201014'],
            ["Assemble/AssembleEvaluation", "", "组装结果评估", 0, '201015'],
            ["Assemble/AssembleEvaluation/Optimize_before", "", "原始组装结果评估", 0, '201016'],
            ["Assemble/AssembleEvaluation/Optimize_before/Trinity_stat.xls", "", "原始组装结果评估表", 0, '201017'],
            ["Assemble/AssembleEvaluation/Optimize_after", "", "优化组装结果评估", 0, '201018'],
            ["Assemble/AssembleEvaluation/Optimize_after/Trinity_stat.xls", "", "优化组装结果评估表", 0, '201019'],
            ["Assemble/LengthDistribution", "", "转录本长度分布", 0, '201020'],
            ["Assemble/LengthDistribution/unigene_count_stat_500.txt", "", "unigene长度分布", 0, '201021'],
            ["Assemble/LengthDistribution/transcript_count_stat_500.txt", "", "transcript长度分布", 0, '201022'],
            ["Assemble/Mapping", "", "测序数据与组装结果比对", 0, '201023'],
            ["Assemble/Mapping/alignment_rate.txt", "", "测序数据与组装结果比对结果表", 0, '201024'],
            ["Annotation", "", "功能注释", 0, '201025'],
            ["Annotation/AnnoStat", "", "功能注释统计", 0, '201026'],
            ["Annotation/AnnoStat/all_annotation_statistics.xls", "", "注释概况统计表", 0, '201027'],
            ["Annotation/AnnoStat/KEGG", "", "与KEGG数据库比对获得注释文件", 0, '201028'],
            ["Annotation/AnnoStat/GO", "", "与GO数据库比对获得注释文件", 0, '201029'],
            ["Annotation/AnnoStat/NR", "", "与NR数据库比对获得注释文件", 0, '201030'],
            ["Annotation/AnnoStat/COG", "", "与COG数据库比对获得注释文件", 0, '201031'],
            ["Annotation/AnnoStat/Pfam", "", "与Pfam数据库比对获得注释文件", 0, '201032'],
            ["Annotation/AnnoStat/Swiss-Prot", "", "与Swiss-Prot数据库比对获得注释文件", 0, '201033'],
            ["Annotation/blast_xml", "", "blast结果文件", 0, '201034'],
            ["Annotation/NR", "", "NR注释", 0, '201035'],
            ["Annotation/NR/unigene_nr_anno_detail.xls", "", "unigene NR注释详情表", 0, '201036'],
            ["Annotation/NR/transcript_nr_anno_detail.xls", "", "transcript NR注释详情表", 0, '201037'],
            ["Annotation/Swiss-Prot", "", "Swiss-Prot注释", 0, '201038'],
            ["Annotation/Swiss-Prot/unigene_swissprot_anno_detail.xls", "", "unigene SwissProt注释详情表", 0, '201039'],
            ["Annotation/Swiss-Prot/transcript_swissprot_anno_detail.xls", "", "transcript SwissProt注释详情表", 0,
             '201040'],
            ["Annotation/Pfam", "", "Pfam注释", 0, '201041'],
            ["Annotation/Pfam/unigene_pfam_anno_detail.xls", "", "unigene Pfam注释详情表", 0, '201042'],
            ["Annotation/Pfam/transcript_pfam_anno_detail.xls", "", "transcript Pfam注释详情表", 0, '201043'],
            ["Annotation/COG", "", "COG注释", 0, '201044'],
            ["Annotation/COG/unigene_cog_list.xls", "", "unigene COG注释列表", 0, '201045'],
            ["Annotation/COG/transcript_cog_list.xls", "", "transcript COG注释列表", 0, '201046'],
            ["Annotation/COG/unigene_cog_anno_detail.xls", "", "unigene COG注释详情表", 0, '201047'],
            ["Annotation/COG/transcript_cog_anno_detail.xls", "", "transcript COG注释详情表", 0, '201048'],
            ["Annotation/COG/unigene_cog_summary.xls", "", "unigene COG分类统计表", 0, '201049'],
            ["Annotation/COG/transcript_cog_summary.xls", "", "transcript COG分类统计表", 0, '201050'],
            ["Annotation/GO", "", "GO注释", 0, '201051'],
            ["Annotation/GO/unigene_gos.list", "", "unigene GO注释列表", 0, '201052'],
            ["Annotation/GO/transcript_gos.list", "", "transcript GO注释列表", 0, '201053'],
            ["Annotation/GO/unigene_blast2go.annot", "", "unigene GO注释详情表", 0, '201054'],
            ["Annotation/GO/transcript_blast2go.annot", "", "transcript GO注释详情表", 0, '201055'],
            ["Annotation/GO/unigene_go12level_statistics.xls", "", "unigene GO12层级注释详情表", 0, '201056'],
            ["Annotation/GO/unigene_go123level_statistics.xls", "", "unigene GO123层级注释详情表", 0, '201057'],
            ["Annotation/GO/unigene_go1234level_statistics.xls", "", "unigene GO1234层级注释详情表", 0, '201058'],
            ["Annotation/GO/transcript_go12level_statistics.xls", "", "transcript GO12层级注释详情表", 0, '201059'],
            ["Annotation/GO/transcript_go123level_statistics.xls", "", "transcript GO123层级注释详情表", 0, '201060'],
            ["Annotation/GO/transcript_go1234level_statistics.xls", "", "transcript GO1234层级注释详情表", 0, '201061'],
            ["Annotation/KEGG", "", "KEGG注释", 0, '201062'],
            ["Annotation/KEGG/unigene_kegg_table.xls", "", "unigene kegg注释详情表", 0, '201063'],
            ["Annotation/KEGG/unigene_pathway_table.xls", "", "unigene pathway注释详情表", 0, '201064'],
            ["Annotation/KEGG/unigene_kegg_layer.xls", "", "unigene kegg pathway层级", 0, '201065'],
            ["Annotation/KEGG/transcript_kegg_table.xls", "", "transcript kegg注释详情表", 0, '201066'],
            ["Annotation/KEGG/transcript_pathway_table.xls", "", "transcript pathway注释详情表", 0, '201067'],
            ["Annotation/KEGG/transcript_kegg_layer.xls", "", "transcript kegg pathway层级", 0, '201068'],
            ["Annotation/KEGG/unigene_pathways", "", "unigene pathway图片展示", 0, '201069'],
            ["Annotation/KEGG/transcript_pathways", "", "transcript pathway图片展示", 0, '201070'],
            ["AnnoQuery", "", "功能查询", 0, '201071'],
            ["AnnoQuery/unigene_anno_detail.xls", "", "unigene功能查询详情表", 0, '201072'],
            ["AnnoQuery/transcript_anno_detail.xls", "", "transcript功能查询详情表", 0, '201073'],
            ["Express", "", "表达量分析", 0, '201074'],
            ["Express/ExpAnnalysis", "", "表达量统计", 0, '201075'],
            ["Express/ExpAnnalysis/unigene.count.matrix.xls", "", "unigene count表达定量结果表", 0, '201076'],
            ["Express/ExpAnnalysis/unigene.tpm.matrix.xls", "", "unigene tpm表达定量结果表", 0, '201077'],
            ["Express/ExpAnnalysis/unigene.fpkm.matrix.xls", "", "unigene fpkm表达定量结果表", 0, '201078'],
            ["Express/ExpAnnalysis/transcript.count.matrix.xls", "", "transcript count表达定量结果表", 0, '201079'],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.xls", "", "transcript tpm表达定量结果表", 0, '201080'],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.xls", "", "transcript fpkm表达定量结果表", 0, '201081'],
            ["Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls", "", "unigene tpm表达定量注释结果表", 0, '201082'],
            ["Express/ExpAnnalysis/unigene.fpkm.matrix.annot.xls", "", "unigene fpkm表达定量注释结果表", 0, '201083'],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls", "", "transcript tpm表达定量注释结果表", 0, '201084'],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量注释结果表", 0, '201085'],
            ["Express/ExpCorr", "", "样本间相关性分析", 0, '201086'],
            ["Express/ExpCorr/sample_correlation.xls ", "", "样本间相关性分析矩阵表", 0, '201087'],
            ["Express/ExpPCA", "", "样本间PCA分析", 0, '201088'],
            ["Express/ExpPCA/PCA.xls", "", "样本间PCA分析结果表", 0, '201089'],
            ["Express/ExpPCA/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表", 0, '201090'],
            ["DiffExpress", "", "表达量差异分析", 0, '201091'],
            ["DiffExpress/DiffExpress_unigene", "", "unigene差异表达分析结果", 0, '201092'],
            ["DiffExpress/DiffExpress_transcript", "", "transcript差异表达分析结果", 0, '201093'],
            ["SNP", "", "SNP分析", 0, '201094'],
            ["SNP/call.vcf", "", "SNP原始结果文件", 0, '201095'],
            ["SNP/snp_detail.xls", "", "SNP结果详情表", 0, '201096'],
            ["SNP/snp_depth_statistics.xls", "", "SNP测序深度统计表", 0, '201097'],
            ["SNP/snp_homo_hete_statistics.xls", "", "SNP类型统计表", 0, '201098'],
            ["SNP/snp_transition_tranversion_statistics.xls", "", "SNP位点统计表", 0, '201099'],
            ["SSR", "", "SSR分析", 0, '201100'],
            ["SSR/ssr_type.xls", "", "SSR类型统计表", 0, '201101'],
            ["SSR/ssr_repeats.xls", "", "SSR类型分布统计表", 0, '201102'],
            ["SSR/ssr_analysis_details.xls", "", "SSR分析详情表", 0, '201103'],
            ["CDS", "", "CDS预测", 0, '201104'],
            ["CDS/Trinity.filter.fasta.transdecoder.bed", "", "转录本序列cds预测结果bed文件", 0, '201105'],
            ["CDS/Trinity.filter.fasta.transdecoder.cds", "", "转录本序列cds预测结果文件", 0, '201106'],
            ["CDS/Trinity.filter.fasta.transdecoder.pep", "", "转录本序列cds预测结果蛋白文件", 0, '201107'],
            ["CDS/pfam_domain", "", "转录本序列pfam注释结果", 0, '201108'],
            ["CDS/unigene_cds_len.xls", "", "unigene cds长度分布", 0, '201109'],
            ["CDS/transcript_cds_len.xls", "", "transcript cds长度分布", 0, '201110'],
            ["TF", "", "转录因子分析", 0, '201111'],
            ["TF/unigene_tf_predict.xls", "", "unigene转录因子预测结果文件", 0, '201112'],
            ["TF/transcript_tf_predict.xls", "", "transcript转录因子预测结果文件", 0, '201113'],
        ]
        sdir = self.add_upload_dir(target_dir)
        sdir.add_relpath_rules(repaths)

    def run_api(self, test=False):
        greenlets_list_first = []
        greenlets_list_sec = []
        # greenlets_list_third = []
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.denovo_task_info')
        task_info.add_task_info()
        self.logger.info("进行第一阶段导表")
        greenlets_list_first.append(gevent.spawn(self.export_qc))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_assembly))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_align))
        gevent.joinall(greenlets_list_first)
        self.logger.info("进行第二阶段导表")
        greenlets_list_sec.append(gevent.spawn(self.export_denovo_annotation))
        if self.option("sample_num") == "multiple":
            if self.option("is_snp") == "True":
                greenlets_list_sec.append(gevent.spawn(self.export_snp))
        greenlets_list_sec.append(gevent.spawn(self.export_ssr))
        if self.option("tf_database") != "Other":
            greenlets_list_sec.append(gevent.spawn(self.export_tf))
        greenlets_list_sec.append(gevent.spawn(self.export_cds))
        gevent.joinall(greenlets_list_sec)
        self.logger.info("进行第三阶段导表")
        # greenlets_list_third.append(gevent.spawn(self.export_expression))
        # gevent.joinall(greenlets_list_third)
        self.export_expression()
        if self.option("report_img"):
            self.export_report_img()
        self.logger.info("导表完成")

    def export_test(self):
        gevent.sleep()
        self.api_qc = self.api.ref_rna_qc
        from bson import ObjectId
        self.group_id = ObjectId("59ae0a75a4e1af55d523f91a")
        self.api_qc.add_control_group(self.option("control").prop["path"], self.group_id)

    @time_count
    def export_qc(self):
        gevent.sleep()
        self.api_qc = self.api.api("denovo_rna_v2.denovo_rna_qc")
        qc_stat = self.hiseq_reads_stat_raw.output_dir
        fq_type = self.option("fq_type").lower()
        if self.option("datatype") == "rawdata":
            self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before", group=self.option('group').path)
        # self.api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="before")
        quality_stat_after = self.hiseq_reads_stat_clean.output_dir + "/qualityStat"  # 质控数据结果统计
        quality_stat_before = self.hiseq_reads_stat_raw.output_dir + "/qualityStat"  # 原始数据结果统计
        if self.option("datatype") == "rawdata":
            self.api_qc.add_gragh_info(quality_stat_before, "before")
        qc_stat = self.hiseq_reads_stat_clean.output_dir
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="after", group=self.option('group').path)
        # self.api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="after")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        if self.option("group").is_set:
            self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(
                self.option("group").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
            if self.option('productive_table').is_set:
                self.api_qc.add_productive_name(samples=self.option('group').prop['sample'],
                                   productive_table=self.option('productive_table').path)
        if self.option("control").is_set:
            self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control").prop["path"],
                                                                            self.group_id)
            self.compare_detail = compare_detail
        # if self.option("is_snp") == True:
        intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results/')
        self.api_qc.add_bam_path(intermediate_dir)

    @time_count
    # def export_denovo_assembly(self):
    #     '''
    #     导入组装结果表格 liubinxu
    #     '''
    #     gevent.sleep()
    #     self.api_denovo_assembly = self.api.api("denovo_rna_v2.denovo_assemble")
    #     result_dir = self.assemble.output_dir
    #     if self.option("optimize") == True:
    #         evolution_dir = self.assemble_evaluation.output_dir
    #         self.api_denovo_assembly.run(result_dir, evolution_dir)
    #     else:
    #         self.api_denovo_assembly.run(result_dir)

    def export_denovo_assembly(self):
        '''
        导入组装结果表格 liubinxu
        '''
        gevent.sleep()
        self.api_denovo_assembly = self.api.api("denovo_rna_v2.denovoass_assemble2")
        assemble_filter_dir = self.assemble_filter.output_dir
        filter_evolution = self.filter_evaluation.output_dir

        if self.option("level").lower() == "transcript":
            unigene_evaluation = self.unigene_evaluation.output_dir
            filter_unigene_evaluation = self.filter_unigene_evaluation.output_dir
        else:
            unigene_evaluation = assemble_filter_dir
            filter_unigene_evaluation = filter_evolution

        result_dir = self.assemble.output_dir
        self.logger.info("{}".format(self.option("optimize")))
        if self.option("optimize") == "True" or self.option("optimize") == True:
            self.api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation, filter_evolution,
                                         filter_unigene_evaluation)
        else:
            self.api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation)

    @time_count
    def export_denovo_align(self):
        '''
        导入mapping率统计表格 liubinxu
        '''
        gevent.sleep()
        api_denovo_align = self.api.api("denovo_rna_v2.denovo_align")
        result_dir = self.align.output_dir
        api_denovo_align.run(result_dir, group=self.option('group').path)

    @time_count
    def export_denovo_annotation(self):
        '''
        导入注释结果表格 liubinxu
        '''
        gevent.sleep()
        self.api_denovo_annotation = self.api.api("denovo_rna_v2.denovo_annotation")
        result_dir = self.annot_stat.output_dir
        trans2gene = os.path.join(result_dir, "all_tran2gene.txt")
        if self.option("express_method") == "RSEM":
            exp_output = self.align.output_dir
        else:
            exp_output = self.express.output_dir
        gene_exp = os.path.join(exp_output, 'gene.tpm.matrix')
        trans_exp = os.path.join(exp_output, 'transcript.tpm.matrix')
        params = {
            "nr_evalue": str(self.option("nr_evalue")),
            "nr_similarity": str(0),
            "nr_identity": str(0),
            "swissprot_evalue": str(self.option("swissprot_evalue")),
            "swissprot_similarity": str(0),
            "swissprot_identity": str(0),
            "cog_evalue": str(self.option("string_evalue")),
            "cog_similarity": str(0),
            "cog_identity": str(0),
            "kegg_evalue": str(self.option("kegg_evalue")),
            "kegg_similarity": str(0),
            "kegg_identity": str(0),
            "pfam_evalue": str(self.option("pfam_evalue")),
            "submit_location": "annotationstat",
            "task_id": self.task_id,
            "task_type": 2,
        }
        self.api_denovo_annotation.anno_type = 'origin'
        self.api_denovo_annotation.run(result_dir, trans2gene, params, taxon=self.option("kegg_database"), version="v2",
                                       exp_level="T", gene_exp=gene_exp, trans_exp=trans_exp)
        # self.api_denovo_annotation.run(result_dir, trans2gene, params, taxon=self.option("kegg_database"))

    @time_count
    def export_assembly(self):
        gevent.sleep()
        self.api_assembly = self.api.api("ref_rna.ref_assembly")
        if self.option("assemble_method") == "cufflinks":
            all_gtf_path = self.assembly.output_dir + "/Cufflinks"
            merged_path = self.assembly.output_dir + "/Cuffmerge"
        else:
            all_gtf_path = self.assembly.output_dir + "/Stringtie"
            merged_path = self.assembly.output_dir + "/StringtieMerge"
        self.api_assembly.add_assemble_result(all_gtf_path=all_gtf_path, merged_path=merged_path,
                                              statistics_path=self.assembly.output_dir + "/Statistics")

    # -----------------------add by gdq-------------------------------------------------------------------------------------
    @time_count
    def export_expression(self):
        gevent.sleep()
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        if self.option("group").is_set:
            group_dict = self.option('group').prop['group_dict']
            group_id = self.group_id
        else:
            group_dict = None
            group_id = None
        if self.option("control").is_set:
            control_id = self.control_id
        else:
            control_id = None
        quant_method = self.option('express_method')
        task_id = self.task_id
        project_sn = self.project_sn

        # add exp matrix
        ## 还需确认路径信息，因为后续交互需要用到express.output_dir
        if self.option("express_method") == "RSEM":
            exp_output = self.align.output_dir
        else:
            exp_output = self.express.output_dir
        if self.option("express_method") == "RSEM" and self.option("exp_way") == "fpkm":
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
                    if self.option("fq_type") == "SE":
                        self.libtype = "fr"
                    else:
                        self.libtype = "f"
            else:
                self.libtype = None

            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            if self.option("sample_num") == "multiple":
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=self.libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=True,
                                               exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=self.libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=False,
                                               exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)

            exp_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            if self.option("sample_num") == "multiple":
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=self.libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=True,
                                              exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=self.libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=False,
                                              exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)

        else:
            params = dict(
                task_id=task_id,
                submit_location="exp_detail",
                task_type=2,
                method=quant_method,
                exp_type='TPM',
            )
            if self.option("strand_specific") == True:
                if self.option('strand_dir') == 'forward':
                    if self.option("fq_type") == "PE":
                        self.libtype = "rf"
                    else:
                        self.libtype = "r"
                else:
                    if self.option("fq_type") == "SE":
                        self.libtype = "fr"
                    else:
                        self.libtype = "f"
            else:
                self.libtype = None

            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            if self.option("sample_num") == "multiple":
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=self.libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=True,
                                               exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=self.libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=False,
                                               exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)

            exp_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
            if self.option("sample_num") == "multiple":
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=self.libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=True,
                                              exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=self.libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=False,
                                              exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)

        self.exp_ids = {
            "G": gene_exp_id
        }

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
                draw_in_groups="no"
            )
            all_exp.add_exp_corr2(corr_output, exp_level='G', quant_method=quant_method, params=params,
                                  project_sn=project_sn, task_id=task_id)

            # add gene venn
            if len(self.option('group').prop['group_dict']) > 1:
                graph_table = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
                group_dict = self.option('group').prop['group_dict']
                if len(group_dict) > 6:
                    group_dict = OrderedDict(group_dict.items()[:6])
                params = json.dumps(dict(
                    task_id=self.task_id,
                    submit_location='expvenn',
                    task_type=2,
                    exp_id=str(gene_exp_id),
                    group_id=str(group_id),
                    exp_level='G',
                    group_dict=group_dict,
                    threshold=1,
                ), sort_keys=True, separators=(',', ':'))
                import datetime
                time_now = datetime.datetime.now()
                name = 'ExpVenn_G_{}_{}_{}'.format(
                    self.option('express_method'), self.option('exp_way').upper(), time_now.strftime('%Y%m%d_%H%M%S'))
                main_info = dict(
                    project_sn=self.project_sn,
                    task_id=self.task_id,
                    version='v2',
                    name=name,
                    created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                    exp_id=str(gene_exp_id),
                    desc='Expression venn analysis main table',
                    params=params,
                    status='start'
                )
                main_id = all_exp.create_db_table('sg_exp_venn', [main_info])
                all_exp.add_exp_venn(graph_table, main_id=main_id)

            # add gene pca
            if self.option("group").prop["sample_number"] > 2:
                pca_output = self.exp_pca.output_dir
                params = dict(
                    task_id=task_id,
                    submit_location="exppca",
                    task_type=2,
                    exp_id=str(gene_exp_id),
                    group_id=str(group_id),
                    exp_level="G",
                    group_dict=self.option('group').prop['group_dict'],
                    draw_in_groups="no"
                    # quant_method=quant_method,
                )
                main_id = all_exp.add_exp_pca2(pca_output, quant_method=quant_method, exp_id=gene_exp_id, exp_level="G",
                                               params=params, project_sn=project_sn, task_id=task_id)
                if hasattr(self, 'ellipse'):
                    all_exp.insert_ellipse_table(os.path.join(self.ellipse.work_dir, 'ellipse_out.xls'), main_id)

            # add transcript diff
            if self.option("level").lower() == "transcript":
                diff_output = self.diffexpress_trans.tool.output_dir
                uniform_output = self.diffexpress_trans.uniform.output_dir
                exp_id, exp_level = trans_exp_id, 'T'
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
                    group_dict=self.option('group').prop['group_dict'],
                    fc=str(float(self.option('fc'))),
                    stat_type=stat_type,
                    stat_cutoff=self.option('diff_fdr_ci'),
                    # correct_method=self.option('padjust_way'),
                    # quant_method=quant_method,
                    # filter_method="no",
                    is_batch="False",
                    diff_method=diff_method,
                )
                if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                    params.update({"correct_method":self.option('padjust_way')})
                all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                                group_dict=group_dict, group_id=group_id,
                                                exp_level=exp_level, quant_method=quant_method,
                                                diff_method=diff_method,
                                                project_sn=project_sn, task_id=task_id, params=params,
                                                pvalue_padjust=stat_type
                                            )
                    # all_exp.add_diffexp_new(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
                    #                     exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
                    #                     project_sn=project_sn, task_id=task_id, params=params,
                    #                     pvalue_padjust=stat_type)
                # else:
                #     all_exp.add_diffexp_new_noiseq(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
                #                             exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
                #                             project_sn=project_sn, task_id=task_id, params=params,
                #                             pvalue_padjust=stat_type)

            # add gene diff
            diff_output = self.diffexpress_gene.tool.output_dir
            uniform_output = self.diffexpress_gene.uniform.output_dir
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
                group_dict=self.option('group').prop['group_dict'],
                fc=str(float(self.option('fc'))),
                # correct_method=self.option('padjust_way'),
                stat_type=stat_type,
                stat_cutoff=self.option('diff_fdr_ci'),
                # quant_method=quant_method,
                # filter_method="no",
                is_batch="False",
                diff_method=diff_method,
            )
            if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                params.update({"correct_method": self.option('padjust_way')})

            self.diff_id = all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                            group_dict=group_dict, group_id=group_id,
                                            exp_level=exp_level, quant_method=quant_method,
                                            diff_method=diff_method,
                                            project_sn=project_sn, task_id=task_id, params=params,
                                            pvalue_padjust=stat_type
                                         )
            self.export_diff_geneset_analysis()
            #     all_exp.add_diffexp_new(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
            #                         exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
            #                         project_sn=project_sn, task_id=task_id, params=params,
            #                         pvalue_padjust=stat_type)
            # else:
            #     all_exp.add_diffexp_new_noiseq(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
            #                             exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
            #                             project_sn=project_sn, task_id=task_id, params=params,
            #                             pvalue_padjust=stat_type)

        # transcript and gene mapping file
        # t2g_file = os.path.join(self.assemble.output_dir, "Trinity.filter.gene_trans_map")
        # all_exp.add_t2g_info(t2g_file, project_sn=project_sn, task_id=task_id)

    def export_diff_geneset_analysis(self):
        self.export_temporary = os.path.join(self.work_dir,"temporary")
        if os.path.exists(self.export_temporary):
            shutil.rmtree(self.export_temporary)
        os.makedirs(os.path.join(self.work_dir,"temporary"))
        diff_geneset_pipline_result = self.diff_geneset_analysis.output_dir
        diff_id = self.diff_id
        task_id = self.task_id
        analysis_names = ["kegg","go", "cog"]
        file_json_path = os.path.join(self.diff_geneset_analysis.file_prepare.output_dir, "prepare_json")
        with open(file_json_path,"r") as j:
            file_dict = json.load(j)
        kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
        all_annot_path = os.path.join(self.annot_stat.output_dir,"all_annot.xls")
        all_annot_df  = pd.read_table(all_annot_path)
        annot_df = all_annot_df[["gene_id", "KO_name", "nr"]]
        annot_df.to_csv(os.path.join(self.export_temporary,"gene_detail"),sep="\t",index=False)
        gene_detail = os.path.join(self.export_temporary,"gene_detail")
        api = self.api.api('denovo_rna_v2.diff_geneset_work_pipline')
        api.add_diff_genest_pipline_table(diff_geneset_pipline_result,diff_id = diff_id, task_id=task_id , analysis_names = analysis_names,
                                          kegg_level_path=kegg_level_path,inter_path= self.export_temporary,exp_id =self.exp_ids['G'])

    def export_report_img(self):
        report_config = os.path.join(self.chart.work_dir, 'report_config.json')
        api =  self.api.api('denovo_rna_v2.report_model')
        s3 = self._sheet.output.split(":")[0]
        report_img_s3 = s3 + ":commonbucket/files/report_img/drna/" + self.task_id
        api.add_report_image(self.task_id, report_config, report_img_s3)


    def merge_annotation_exp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        if self.option("express_method") == "RSEM":
            exp_output = self.align.output_dir
        else:
            exp_output = self.express.output_dir
        group_dict = self.option('group').prop['group_dict']
        annot_output = self.annot_stat.output_dir
        all_annot = os.path.join(annot_output, 'all_annot.xls')
        all_annot = pd.read_table(all_annot, header=0)
        pattern = '.*?[(](.*?)[)]'
        all_annot["nr_description"] = all_annot["nr"].str.extract(pattern)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index(
            'gene_id')
        trans_annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript')
        order = ["go", "KO_id", "KO_name", "paths", "cog", "cog_description", "swissprot", "pfam"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][["gene_id", 'nr_description']].set_index('gene_id')
        trans_annot_pd = pd.DataFrame(trans_annot_pd, columns=order)
        trans_info_pd = all_annot[['transcript', 'gene_id', 'nr_description']].set_index('transcript')

        # gene
        ## gene tpm
        gene_tpm_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
        gene_tpm_pd = pd.read_table(gene_tpm_matrix, header=0, index_col=0)
        gene_tpm_dicts = gene_tpm_pd.to_dict('index')
        gene_group_tpm = OrderedDict()
        for seq_id in sorted(gene_tpm_dicts):
            tmp_exp_dict = gene_tpm_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_tpm:
                    gene_group_tpm[seq_id] = dict()
                gene_group_tpm[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_tpm, columns=gene_group_tpm.keys())).T
        gene_tpm_result = pd.concat([gene_info_pd, gene_tpm_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_tpm_out = os.path.join(exp_output, 'gene.tpm.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_tpm_result.columns.tolist())
        with open(gene_tpm_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_tpm_result.to_csv(gene_tpm_out, header=False, index=True, sep='\t', mode='a')

        ## gene count
        gene_count_matrix = os.path.join(exp_output, 'gene.count.matrix')
        gene_count_pd = pd.read_table(gene_count_matrix, header=0, index_col=0)
        gene_count_dicts = gene_count_pd.to_dict('index')
        gene_group_count = OrderedDict()
        for seq_id in sorted(gene_count_dicts):
            tmp_exp_dict = gene_count_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_count:
                    gene_group_count[seq_id] = dict()
                gene_group_count[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_count, columns=gene_group_count.keys())).T
        gene_count_result = pd.concat([gene_info_pd, gene_count_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_count_out = os.path.join(exp_output, 'gene.count.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_count_result.columns.tolist())
        with open(gene_count_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_count_result.to_csv(gene_count_out, header=False, index=True, sep='\t', mode='a')

        ## gene fpkm
        if self.option('express_method') == 'RSEM':
            gene_fpkm_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            gene_fpkm_pd = pd.read_table(gene_fpkm_matrix, header=0, index_col=0)
            gene_fpkm_dicts = gene_fpkm_pd.to_dict('index')
            gene_group_fpkm = OrderedDict()
            for seq_id in sorted(gene_fpkm_dicts):
                tmp_exp_dict = gene_fpkm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in gene_group_fpkm:
                        gene_group_fpkm[seq_id] = dict()
                    gene_group_fpkm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            gene_group_pd = (pd.DataFrame(data=gene_group_fpkm, columns=gene_group_fpkm.keys())).T
            gene_fpkm_result = pd.concat([gene_info_pd, gene_fpkm_pd, gene_group_pd, gene_annot_pd], axis=1)
            gene_fpkm_out = os.path.join(exp_output, 'gene.fpkm.matrix.annot.xls')
            header = ['gene_id']
            header.extend(gene_fpkm_result.columns.tolist())
            with open(gene_fpkm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            gene_fpkm_result.to_csv(gene_fpkm_out, header=False, index=True, sep='\t', mode='a')
        if self.option('level').lower() == 'transcript':
            ## transcript tpm
            transcript_tpm_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            transcript_tpm_pd = pd.read_table(transcript_tpm_matrix, header=0, index_col=0)
            transcript_tpm_dicts = transcript_tpm_pd.to_dict('index')
            transcript_group_tpm = OrderedDict()
            for seq_id in sorted(transcript_tpm_dicts):
                tmp_exp_dict = transcript_tpm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_tpm:
                        transcript_group_tpm[seq_id] = dict()
                    transcript_group_tpm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_tpm, columns=transcript_group_tpm.keys())).T
            transcript_tpm_result = pd.concat([trans_info_pd, transcript_tpm_pd, transcript_group_pd, trans_annot_pd],
                                              axis=1)
            transcript_tpm_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_tpm_result.columns.tolist())
            with open(transcript_tpm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_tpm_result.to_csv(transcript_tpm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript fpkm
            if self.option('express_method') == 'RSEM':
                transcript_fpkm_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
                transcript_fpkm_pd = pd.read_table(transcript_fpkm_matrix, header=0, index_col=0)
                transcript_fpkm_dicts = transcript_fpkm_pd.to_dict('index')
                transcript_group_fpkm = OrderedDict()
                for seq_id in sorted(transcript_fpkm_dicts):
                    tmp_exp_dict = transcript_fpkm_dicts[seq_id]
                    for group in group_dict:
                        if seq_id not in transcript_group_fpkm:
                            transcript_group_fpkm[seq_id] = dict()
                        transcript_group_fpkm[seq_id].update({group: round(
                            sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
                transcript_group_pd = (pd.DataFrame(data=transcript_group_fpkm, columns=transcript_group_fpkm.keys())).T
                transcript_fpkm_result = pd.concat(
                    [trans_info_pd, transcript_fpkm_pd, transcript_group_pd, trans_annot_pd], axis=1)
                transcript_fpkm_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
                header = ['transcript_id']
                header.extend(transcript_fpkm_result.columns.tolist())
                with open(transcript_fpkm_out, "w") as w:
                    w.write("\t".join(header) + "\n")
                transcript_fpkm_result.to_csv(transcript_fpkm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript count
            transcript_count_matrix = os.path.join(exp_output, 'transcript.count.matrix')
            transcript_count_pd = pd.read_table(transcript_count_matrix, header=0, index_col=0)
            transcript_count_dicts = transcript_count_pd.to_dict('index')
            transcript_group_count = OrderedDict()
            for seq_id in sorted(transcript_count_dicts):
                tmp_exp_dict = transcript_count_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_count:
                        transcript_group_count[seq_id] = dict()
                    transcript_group_count[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_count, columns=transcript_group_count.keys())).T
            transcript_count_result = pd.concat(
                [trans_info_pd, transcript_count_pd, transcript_group_pd, trans_annot_pd], axis=1)
            transcript_count_out = os.path.join(exp_output, 'transcript.count.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_count_result.columns.tolist())
            with open(transcript_count_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_count_result.to_csv(transcript_count_out, header=False, index=True, sep='\t', mode='a')

        # gene_exp_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
        # trans_exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
        # annot_output = self.annot_stat.output_dir
        # all_annot = os.path.join(annot_output, 'all_annot.xls')
        # all_annot = pd.read_table(all_annot, header=0)
        # gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index('gene_id')
        # trans_annot_pd=all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript')
        # # for gene
        # gene_pd = pd.read_table(gene_exp_matrix, header=0, index_col=0)
        # # gene_annot_pd = pd.read_table(gene_annot, header=0, index_col=0)
        # #for count
        # gene_count_matrix=os.path.join(exp_output,"gene.count.matrix")
        # gene_count_pd=pd.read_table(gene_count_matrix, header=0, index_col=0)
        # gene_count_result = pd.concat([gene_count_pd, gene_annot_pd], axis=1)
        # gene_count_out = os.path.join(exp_output, 'gene.count.matrix.annot.xls')
        # gene_count_result.to_csv(gene_count_out, header=True, index=True, sep='\t')
        # # for tpm
        # gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1)
        # gene_out = os.path.join(exp_output, 'gene.tpm.matrix.annot.xls')
        # gene_result.to_csv(gene_out, header=True, index=True, sep='\t')
        # # for fpkm
        # if self.option("express_method") == "RSEM":
        #     gene_exp_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
        #     gene_pd = pd.read_table(gene_exp_matrix, header=0, index_col=0)
        #     gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1)
        #     gene_out = os.path.join(exp_output, 'gene.fpkm.matrix.annot.xls')
        #     gene_result.to_csv(gene_out, header=True, index=True, sep='\t')
        # # for transcript
        # gene_pd = pd.read_table(trans_exp_matrix, header=0, index_col=0)
        # # for count
        # gene_count_matrix = os.path.join(exp_output, "transcript.count.matrix")
        # gene_count_pd = pd.read_table(gene_count_matrix, header=0, index_col=0)
        # gene_count_result = pd.concat([gene_count_pd, gene_annot_pd], axis=1)
        # gene_count_out = os.path.join(exp_output, 'transcript.count.matrix.annot.xls')
        # gene_count_result.to_csv(gene_count_out, header=True, index=True, sep='\t')
        # # for tpm
        # gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1)
        # gene_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
        # gene_result.to_csv(gene_out, header=True, index=True, sep='\t')
        # # for fpkm
        # if self.option("express_method") == "RSEM":
        #     trans_exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
        #     gene_pd = pd.read_table(trans_exp_matrix, header=0, index_col=0)
        #     gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1)
        #     gene_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
        #     gene_result.to_csv(gene_out, header=True, index=True, sep='\t')

    def merge_annotation_diffexp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        annot_output = self.annot_stat.output_dir
        all_annot = os.path.join(annot_output, 'all_annot.xls')
        all_annot = pd.read_table(all_annot, header=0)
        pattern = '.*?[(](.*?)[)]'
        all_annot["nr_description"] = all_annot["nr"].str.extract(pattern)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index(
            'gene_id')
        trans_annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript')
        order = ["go", "KO_id", "KO_name", "paths", "cog", "cog_description", "swissprot", "pfam"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][["gene_id", 'nr_description']].set_index(
            'gene_id')
        trans_annot_pd = pd.DataFrame(trans_annot_pd, columns=order)
        trans_info_pd = all_annot[['transcript', 'gene_id', 'nr_description']].set_index('transcript')

        diff_output = self.diffexpress_gene.output_dir
        ## 防止重运行反复添加注释信息
        # remove_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls')
        # for file in remove_files:
        #     if os.path.exists(file):
        #         os.remove(file)
        duplicate_files = glob.glob(diff_output + '/' + '*.annot.xls') + glob.glob(
            diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(diff_output + '/' + '*_vs_*.sizeFactor.xls')
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + "/*")
        for each in target_files:
            if each.endswith("diff_summary.xls"):
                gene_out = os.path.join(diff_output,
                                        "diff_summary_{}.xls".format(self.option("diff_method"))) + '.annot.xls'
                gene_outraw = os.path.join(diff_output, "diff_summary_{}.xls".format(self.option("diff_method")))
                with open(each, "r")as raw, open(gene_out, "w")as new:
                    head0 = raw.readline()
                    compare_num = (len(head0.strip().split("\t")) - 2) / 2
                    newheadraw = "\t".join(head0.strip().split("\t")[0:compare_num + 1])
                    newhead01 = "gene_id\tnr_description\t" + "\t".join(head0.strip().split("\t")[1:compare_num + 1])
                    newhead0 = "gene_id\tnr_description\t" + "\t".join(
                        head0.strip().split("\t")[1:compare_num + 1]) + "\t" + "\t".join(order)
                    head1 = raw.readline()
                    newhead1 = head1.strip().split("\t")[0] + "\t\t" + "\t".join(
                        head1.strip().split("\t")[1:compare_num + 1])
                    new.write(newheadraw + "\n")
                    for line in raw.readlines():
                        line = line.strip().split("\t")
                        new.write(line[0])
                        for n, com in enumerate(line[1:compare_num + 1]):
                            new.write("\t" + com + "|" + line[n + compare_num + 2])
                        new.write("\n")
                summary_diff = pd.read_table(gene_out, header=0, index_col=0)
                gene_result = pd.concat([gene_info_pd, summary_diff, gene_annot_pd], join='inner', axis=1)
                gene_result1 = pd.concat([gene_info_pd, summary_diff], join='inner', axis=1)
                with open(gene_out, "w")as final:
                    final.write(newhead0 + "\n" + newhead1 + "\n")
                with open(gene_outraw, "w") as final1:
                    final1.write(newhead01 + "\n" + newhead1 + "\n")
                gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')
                gene_result1.to_csv(gene_outraw, header=False, index=True, sep='\t', mode='a')
                #
                #
                # gene_pd = pd.read_table(each)
                # gene_pd=gene_pd.loc[2:,].set_index("seq_id")
                # gene_result = pd.concat([gene_pd, gene_annot_pd], axis=1,join_axes=[gene_pd.index])
                # gene_out = each.split('.xls')[0] + '.annot.xls'
                # gene_result.to_csv(gene_out, header=False, index=True, sep='\t')
                # header=open(each).readlines()[0]+open(each).readlines()[1]
                # with open(gene_out,'r+') as sm:
                #     content=sm.read()
                #     sm.seek(0,0)
                #     sm.write(header+content)
            else:
                gene_pd = pd.read_table(each, header=0, index_col=0)
                gene_result = pd.concat([gene_info_pd, gene_pd, gene_annot_pd], axis=1)
                gene_out = each.split('.xls')[0] + '.annot.xls'
                header = ['gene_id']
                header.extend(gene_result.columns.tolist())
                with open(gene_out, "w") as w:
                    w.write("\t".join(header) + "\n")
                gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')

        if self.option("level").lower() == "transcript":
            # for transcript
            diff_output = self.diffexpress_trans.output_dir
            ## 防止重运行反复添加注释信息
            duplicate_files = glob.glob(diff_output + '/' + '*.annot.xls') + glob.glob(
                diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(diff_output + '/' + '*_vs_*.sizeFactor.xls')
            for file in duplicate_files:
                os.remove(os.path.join(diff_output, file))
            target_files = glob.glob(diff_output + "/*")
            for each in target_files:
                if each.endswith("diff_summary.xls"):
                    gene_out = os.path.join(diff_output,
                                            "diff_summary_{}.xls".format(self.option("diff_method"))) + '.annot.xls'
                    gene_outraw = os.path.join(diff_output, "diff_summary_{}.xls".format(self.option("diff_method")))
                    with open(each, "r")as raw, open(gene_out, "w")as new:
                        head0 = raw.readline()
                        compare_num = (len(head0.strip().split("\t")) - 2) / 2
                        newheadraw = "\t".join(head0.strip().split("\t")[0:compare_num + 1])
                        newhead01 = "transcript_id\tgene_id\tnr_description\t" + "\t".join(
                            head0.strip().split("\t")[1:compare_num + 1])
                        newhead0 = "transcript_id\tgene_id\tnr_description\t" + "\t".join(
                            head0.strip().split("\t")[1:compare_num + 1]) + "\t" + "\t".join(order)
                        head1 = raw.readline()
                        newhead1 = head1.strip().split("\t")[0] + "\t\t" + "\t".join(
                            head1.strip().split("\t")[1:compare_num + 1])
                        new.write(newheadraw + "\n")
                        for line in raw.readlines():
                            line = line.strip().split("\t")
                            new.write(line[0])
                            for n, com in enumerate(line[1:compare_num + 1]):
                                new.write("\t" + com + "|" + line[n + compare_num + 2])
                            new.write("\n")
                    summary_diff = pd.read_table(gene_out, header=0, index_col=0)
                    gene_result = pd.concat([trans_info_pd, summary_diff, trans_annot_pd], join='inner', axis=1)
                    gene_result1 = pd.concat([trans_info_pd, summary_diff], join='inner', axis=1)
                    with open(gene_out, "w")as final:
                        final.write(newhead0 + "\n" + newhead1 + "\n")
                    with open(gene_outraw, "w") as final1:
                        final1.write(newhead01 + "\n" + newhead1 + "\n")
                    gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')
                    gene_result1.to_csv(gene_outraw, header=False, index=True, sep='\t', mode='a')
                    # gene_pd = pd.read_table(each)
                    # gene_pd = gene_pd.loc[2:, ].set_index("seq_id")
                    # gene_result = pd.concat([gene_pd, trans_annot_pd], axis=1,join_axes=[gene_pd.index])
                    # gene_out = each.split('.xls')[0] + '.annot.xls'
                    # gene_result.to_csv(gene_out, header=False, index=True, sep='\t')
                    # header = open(each).readlines()[0] + open(each).readlines()[1]
                    # with open(gene_out, 'r+') as sm:
                    #     content = sm.read()
                    #     sm.seek(0, 0)
                    #     sm.write(header + content)
                else:
                    gene_pd = pd.read_table(each, header=0, index_col=0)
                    gene_result = pd.concat([trans_info_pd, gene_pd, trans_annot_pd], axis=1)
                    gene_out = each.split('.xls')[0] + '.annot.xls'
                    header = ['transcript_id']
                    header.extend(gene_result.columns.tolist())
                    with open(gene_out, "w") as w:
                        w.write("\t".join(header) + "\n")
                    gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')

    @time_count
    def build_seq_database(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_seq = self.api.api("denovo_rna_v2.seq_detail")
        cds = os.path.join(self.cds_predict.output_dir, 'all_predicted.cds.fa')
        pep = os.path.join(self.cds_predict.output_dir, 'all_predicted.pep.fa')
        fasta = os.path.join(self.assemble_filter.output_dir, 'Trinity.filter.fasta')
        trans2unigene = os.path.join(self.cds_predict.output_dir, 'all_tran2gen.txt')
        seq_db = os.path.join(self.work_dir, 'seq_db.sqlite3')
        self.export_seq.build_seq_database(seq_db, cds, pep, fasta, trans2unigene, task_id=self.task_id)
        self.seq_stat_seq = self.api.api("denovo_rna_v2.seq_stat_detail")
        self.logger.info("开始进行sg_seq_stat的导表")
        gene_stat = os.path.join(self.sequence_detail.work_dir, 'gene_stat')
        self.logger.info("gene_stat的路径是{}".format(gene_stat))
        trans_stat = os.path.join(self.sequence_detail.work_dir, 'tran_stat')
        self.logger.info("trans_stat的路径是{}".format(trans_stat))
        self.seq_stat_seq.add_seq_stat(gene_stat, trans_stat)
        self.logger.info("开始进行sg_seq_stat的结束")

    @time_count
    def export_snp(self):
        gevent.sleep()
        api_snpfinal = self.api.api("denovo_rna_v3.snp_api")
        self.logger.info("开始进行Snpfinal的导表")
        task_id = self.task_id
        project_sn = self.project_sn
        call_vcf_path = self.snp.work_dir + '/Snp/' + 'call.vcf'
        params = dict(
            task_id=task_id,
            submit_location="snp_detail",
            task_type=2,
            method=self.option("snp_method")
        )

        if self.snpfinal.work_dir + '/new_snp_rewrite' is None and self.snpfinal.work_dir + '/new_indel_rewrite' is None:
            self.logger.info("此次分析没有call出snp和indel")

        if self.snpfinal.work_dir + '/new_snp_rewrite' is not None and self.snpfinal.work_dir + '/new_indel_rewrite' is not None:
            new_snp_rewrite = self.snpfinal.work_dir + '/new_snp_rewrite'
            new_indel_rewrite = self.snpfinal.work_dir + '/new_indel_rewrite'
            depth_path = self.snpfinal.work_dir + '/depth_new_per'
            hh_path = self.snpfinal.work_dir + '/statis_hh'
            tt_new_per_path = self.snpfinal.work_dir + '/tt_new_per'
            cds_path = self.snpfinal.work_dir + "/statis_cds"
            anno_path = self.snpfinal.work_dir + "/snp_anno_stat"
            api_snpfinal.add_snp(new_snp_rewrite, new_indel_rewrite, depth_path, hh_path, tt_new_per_path, cds_path,
                                 anno_path, project_sn=project_sn, task_id=task_id, params=params, group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出snp和indel")

        if self.snpfinal.work_dir + '/new_snp_rewrite' is not None and self.snpfinal.work_dir + '/new_indel_rewrite' is None:
            new_snp_rewrite = self.snpfinal.work_dir + '/new_snp_rewrite'
            depth_path = self.snpfinal.work_dir + '/depth_new_per'
            hh_path = self.snpfinal.work_dir + '/statis_hh'
            tt_new_per_path = self.snpfinal.work_dir + '/tt_new_per'
            cds_path = self.snpfinal.work_dir + "/statis_cds"
            anno_path = self.snpfinal.work_dir + "/snp_anno_stat"
            api_snpfinal.add_snp(new_snp_rewrite=new_snp_rewrite, new_indel_rewrite=None, depth_path=depth_path,
                                 hh_path=hh_path, tt_new_per_path=tt_new_per_path, cds_path=cds_path,
                                 anno_path=anno_path, project_sn=project_sn, task_id=task_id, params=params, group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出snp,但是没有call出indel")

        if self.snpfinal.work_dir + '/new_snp_rewrite' is None and self.snpfinal.work_dir + '/new_indel_rewrite' is not None:
            new_indel_rewrite = self.snpfinal.work_dir + '/new_indel_rewrite'
            api_snpfinal.add_snp(new_snp_rewrite=None, new_indel_rewrite=new_indel_rewrite, project_sn=project_sn,
                                 task_id=task_id, params=params, group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出indel,但是没有call出snp")

    @time_count
    def export_ssr(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            submit_location="ssr",
            task_type=2,
            rept_1=10,
            rept_2=6,
            rept_3=5,
            rept_4=5,
            rept_5=5,
            rept_6=5,
            ssr_distance=100,
        )
        api_ssr = self.api.api("denovo_rna_v2.ssr")
        self.logger.info("开始进行ssr的导表")
        ssr_statistic_path = self.ssr.work_dir + '/' + 'ssr_type.txt'
        ssr_detail_path = self.ssr.work_dir + '/' + 'tmp.txt'
        ssr_class_path = self.ssr.work_dir + '/' + 'ssr_repeats_class.txt'
        api_ssr.add_ssr(ssr_statistic_path, ssr_detail_path, ssr_class_path, name=None, params=params,
                        project_sn=project_sn, task_id=task_id)

    @time_count
    def export_tf(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            task_type=2,
            search_pfam='True',
            p_length=50,
            Markov_length=3000,
            cpu=20,
            hmmcan1='noali',
            hmmcan2='acc',
            hmmcan3='notextw',
            E=self.option("tf_evalue")
        )
        # self.fasta_name = self.assemble.option("filter_fa").prop["path"].split("/")[-1]
        # bed = '{}.transdecoder.bed'.format(self.fasta_name)
        bed = "all_predicted.bed"
        bedpath = self.cds_predict.work_dir + '/' + bed
        if self.option("tf_database").lower() == "animal":
            tf_unigene_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_unigene_animal')
            tf_transcript_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_transcript_animal')
            api_tf_api = self.api.api("denovo_rna_v2.tf_api")
            self.logger.info("开始进行动物tf的导表")
            api_tf_api.add_tf_unigene(tf_unigene_path, bedpath=bedpath, name=None, params=params,
                                      project_sn=project_sn, task_id=task_id)
            self.logger.info("动物的tf基因水平的导表完成")
            api_tf_api.add_tf_transcript(tf_transcript_path, bedpath=bedpath, name=None, params=params,
                                         project_sn=project_sn, task_id=task_id)
            self.logger.info("动物的tf基因，转录水平的导表完成")

        if self.option("tf_database").lower() == "plant":
            tf_unigene_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_unigene_plant')
            tf_transcript_path = os.path.join(self.cds_predict.work_dir, 'Predict', 'merge_only_transcript_plant')
            api_tf_api = self.api.api("denovo_rna_v2.tf_api")
            self.logger.info("开始进行植物tf的导表")
            api_tf_api.add_tf_unigene(tf_unigene_path, bedpath=bedpath, name=None, params=params,
                                      project_sn=project_sn, task_id=task_id)
            self.logger.info("植物的tf基因水平的导表完成")
            api_tf_api.add_tf_transcript(tf_transcript_path, bedpath=bedpath, name=None, params=params,
                                         project_sn=project_sn, task_id=task_id)
            self.logger.info("植物的tf基因，转录水平的导表完成")

    @time_count
    def export_cds(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            task_type=2,
            search_pfam='True',
            p_length=50,
            Markov_length=3000,
            cpu=20,
            hmmcan1='noali',
            hmmcan2='acc',
            hmmcan3='notextw',
            E=self.option("tf_evalue")
        )
        api_cds = self.api.api("denovo_rna_v2.cdslen")
        self.logger.info("开始进行cds的导表")
        cds_unigene_length = self.cds_predict.output_dir + '/' + 'cds_len_unigene.txt'
        cds_transcript_length = self.cds_predict.output_dir + '/' + 'cds_len_transcript.txt'
        all_predicted = self.cds_predict.output_dir + '/' + 'all_predicted.xls'
        api_cds.cds(all_predicted, project_sn=self.project_sn, task_id=self.task_id,
                    cds_unigene_length=cds_unigene_length, cds_transcript_length=cds_transcript_length, params=params)
        # api_cds.add_cds_unigene_length(cds_unigene_length, name=None, params=params,
        # project_sn=project_sn, task_id=task_id)
        # self.logger.info("开始进行cds基因水平的导表完成")
        # api_cds.add_cds_transcript_length(cds_transcript_length, name=None, params=params,
        # project_sn=project_sn, task_id=task_id)
        # self.logger.info("开始进行cds基因，转录水平的导表完成")


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.denovo_rna_v2.denovorna import DenovornaWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "testdenovo" + str(random.randint(1, 10000)) + "yyyy",
            "type": "workflow",
            "name": "denovo_rna_v2.denovorna",
            "instant": False,
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/denovo_rna_v2/test_files_dirs/fastq_dir",
                group="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/denovo_rna_v2/test_files_dirs/group",
                control="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/denovo_rna_v2/test_files_dirs/control"
            )
        }
        wsheet = Sheet(data=data)
        wf = DenovornaWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    unittest.main()
