# -*- coding:utf-8 -*-
# __author__ = '封一统'
"""原核转录组一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
# import subprocess
import datetime
import shutil
import re
import time
import gevent
import functools
from biocluster.config import Config
import pandas as pd
import unittest
# # from biocluster.wpm.client import *
import json
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.lnc_rna.copy_file import CopyFile
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
import time
from mbio.packages.dna_evolution.send_email import SendEmail
import tarfile
from collections import OrderedDict


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
        workflow option参数设置
        """
        self._sheet = wsheet_object
        super(ProkrnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "strand_specific", "type": "bool", "default": True},  # 当为PE测序时，是否为链特异性
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "align_species", "type": "string", "default": "is_ncbi"},  # 选择参考物种是来着NCBI或用户自己上传
            {"name": "species_name", "type": "string", "default": ""},  # 客户研究对象在NCBI里的命名
            {"name": "genome_id", "type": "string", "default": ""},  # 客户研究对象的参考基因组等在NCBI里的accession号
            {"name": "genome_db", "type": "infile", "format": "prok_rna.fasta"},  # 客户自行上传的参考基因组
            {"name": "gff_or_gtf", "type": "string", "default": "gff"},  # 客户自行上传参考基因组时需要选择是上传gtf或者gff文件, gbk added on 202109
            {"name": "gff_or_gtf_file", "type": "infile", "format": "prok_rna.common"},  # 客户自行上传参考基因组时上传的gtf或者gff文件
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'}, # 上机名称
            {"name": "control_file", "type": "infile", "format": "sample.control_table"},  # 对照表
            {"name": "mapping_stop", "type": "bool", "default": False},  # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_sample_percent", "type": "string", "default": ""},  # 设定比对率不达标样本占总样本数多少时结束流程
            {"name": "mapping_ratio", "type": "string", "default": ""},  # 设定比对率小于多少时判定样本不达标
            {"name": "rrna_stop", "type": "bool", "default": False},  # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_sample_percent", "type": "string", "default": ""},  # 设定rRNA比率超标样本占总样本数多少时结束流程
            {"name": "rrna_ratio", "type": "string", "default": ""},  # 设定rRNA比率大于多少时判定样本不达标
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
            {"name": "kegg_version", "type": "string", "default": "202003"},
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-5},  # Swissprot比对使用的e值
            {"name": "cog_blast_evalue", "type": "float", "default": 1e-5},  # COG比对使用的e值
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},  # PFAM比对使用的e值
            # 表达量分析
            {"name": "express_method", "type": "string", "default": "rsem"},  # 表达量分析手段: Salmon, Kallisto, RSEM
            {"name": "exp_way", "type": "string", "default": "fpkm"},  # fpkm or tpm
            # 表达差异分析
            {"name": "diff_method", "type": "string", "default": "DESeq2"},  # 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  # 选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  # Bonferroni,Holm,BH,BY
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_2019"},
            {'name': 'report_img', 'type': 'bool', 'default': True},
        ]
        # 获取输出目录
        self.workflow_output_tmp = self._sheet.output
        self.logger.info(self.workflow_output_tmp)
        try:
            if re.match(r'tsanger:', self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
                self.workflow_output2 = self.workflow_output_tmp.replace('tsanger:', '')
            elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp
                self.workflow_output2 = self.workflow_output_tmp
            else:
                self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
                self.workflow_output2 = self.workflow_output_tmp.replace('sanger:', '')
        except:
            self.workflow_output = 'rere...'
            self.workflow_output2 = '/mnt/...'
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))
        if self.option('align_species') == 'is_ncbi':
            self.ncbi = True
        else:
            self.ncbi = False
        # 获取数据库中该物种相关文件路径
        # if self.option('align_species') == 'is_ncbi':
        #     self.db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        #     col = self.db["sg_genome_db"]
        #     genome_info = col.find_one(
        #         {"genome_id": self.option("genome_id")})
        #     if not genome_info:
        #         self.logger.info("数据库中不存在该物种注释信息，从NCBI官网中下载该物种注释信息")
        #         self.ncbi = True
        # if self.option('align_species') == 'is_ncbi':
        #     self.db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        #     col = self.db["sg_genome_db"]
        #     genome_info = col.find_one(
        #         {"organism_name": self.option("species_name"), "genome_id": self.option("genome_id")})
        #     if not genome_info:
        #         self.set_error("数据库中不存在该物种注释信息，程序退出", code="15000101")
        self.db_path = self.config.SOFTWARE_DIR + "/database/Transeq_DB_v1/"
        if self.option('align_species') == 'is_ncbi':
            # if self.ncbi is False:
            #     if not os.path.exists(os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna')):
            #         self.set_error("服务器内不存在您输入的NCBI accession号的参考基因组文件，请检查", code="15000102")
            #     if not os.path.exists(os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff')):
            #         self.set_error("服务器内不存在您输入的NCBI accession号的参考基因组的gff文件，请检查", code="15000103")
            # if self.ncbi is True:
            #     pass
            pass
        else:
            if not self.option('genome_db').is_set or not os.path.exists(self.option('genome_db').prop['path']):
                self.set_error("没有找到您上传的参考基因组", code="15000104")
            if not self.option('gff_or_gtf_file').is_set or not os.path.exists(
                    self.option('gff_or_gtf_file').prop['path']):
                self.set_error("没有找到您上传的gff或gtf文件", code="15000105")


        # if self.ncbi is False:
        #     self.ref_genome = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option(
        #         'align_species') == 'is_ncbi' else self.option('genome_db').prop['path']

        # 添加tool/module
        self.genomecheck = self.add_tool("prok_rna.check_genome")
        self.rock_index = self.add_tool("prok_rna.rockhopper_index")
        self.rockhopper = self.add_tool("prok_rna.rockhopper")
        self.filecheck = self.add_tool("prok_rna.file_check")
        # self.qc = self.add_module("prok_rna.hiseq_qc")
        self.qc = self.add_module('datasplit.fastp_rna')
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
        self.express_total = self.add_module("prok_rna.quant")
        self.exp_pca = self.add_tool("prok_rna.exp_pca")
        # self.ellipse = self.add_tool('graph.ellipse')
        self.ellipse = self.add_tool('ref_rna_v3.ellipse')
        self.exp_corr = self.add_tool("prok_rna.exp_corr")
        self.exp_venn = self.add_tool('prok_rna.exp_venn')
        # self.diffexpress = self.add_tool("prok_rna.diffexp")
        self.diffexpress = self.add_module("prok_rna.diffexp_batch_new")  # 新增，由原先的调用一个tool，改为调用一个module
        self.diff_geneset_analysis = self.add_module("prok_rna.workflow_diffgt.diff_geneset_all_pipline")#新增 工作流差异一键化新增模块
        self.srna = self.add_module("prok_rna.srna")
        self.promote = self.add_tool('prok_rna.promote')
        self.terminator = self.add_tool('prok_rna.terminator')  # 新增分析
        self.extract_biotype = self.add_tool('prok_rna.extract_biotype')
        self.extract_mrna = self.add_tool('prok_rna.extract_mrna')
        self.extract_pep = self.add_tool('prok_rna.extract_mrna')
        self.extract_mrna_srna = self.add_tool('prok_rna.extract_mrna')
        self.extract_ptt = self.add_tool('prok_rna.extract_ptt')
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        # 判断流程结束这次采用计数的方法，对于特定的tool，set_output之后这个数会加一，等这个数为7时，流程结束
        # self.end_number = 0
        # if self.option('strand_specific'):
        #     self.to_end = 7
        # else:
        #     self.to_end = 6

        self.final_tools = [self.qc_stat_before, self.map_qc, self.snp, self.exp_pca, self.ellipse, self.exp_corr, self.exp_venn,
                            self.promote, self.terminator, self.srna, self.annot_class, self.express_total, self.express,self.diff_geneset_analysis]
        # self.final_tools.remove(self.snp)

        # 添加step，显示在页面进度条
        all_steps = ['download4ncbi', 'gbk2gff', 'genomecheck',"rock_index", "filecheck", "rna_qc", "mapping", "express", "diffexpress", "snp_rna",
                     "map_qc", "annot_mapdb", "annot_orfpfam", "annot_filter", "annot_class", "exp_pca", "ellipse",
                     "exp_corr", "exp_venn", "qc_stat_before", "qc_stat_after", "rockhopper", "srna", "promote", 'terminator','geneset_analysis'
                     ]
        if self.option("group_table").prop["sample_number"] < 3:
            all_steps.remove('exp_pca')
            all_steps.remove('ellipse')
            self.final_tools.remove(self.exp_pca)
            self.final_tools.remove(self.ellipse)
        else:
            pass
        group_spname = self.option("group_table").get_group_spname()
        group_snum = [len(group_spname[g]) for g in group_spname]
        self.min_group_num = min(group_snum)
        if self.min_group_num >= 3:
            pass
        else:
            try:
                self.final_tools.remove(self.ellipse)
            except:
                pass
        if self.ncbi is True:
            pass
        else:
            all_steps.remove('download4ncbi')
        for step in all_steps:
            self.step.add_steps(step)
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self.task_id, 'prok_rna')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        if self.option('group_table').is_set:
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
            if self.option('control_file').is_set:
                vs_list = self.option("control_file").prop["vs_list"]
                for vs in vs_list:
                    dup = False
                    for i in vs:
                        if len(group_dict[i]) > 1:
                            dup = True
                    if not dup and self.option("diff_method").lower() in ['deseq2']:
                        raise OptionError("差异分析软件{}必须要每一个差异分组均存在重复，差异分组{}不存在重复，"
                                          "请重新选择差异分析软件进行分析。".format(self.option("diff_method"), vs))
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE", code="15000106")
        try:
            nr_evalue = float(self.option("nr_evalue"))
            cog_evalue = float(self.option("cog_evalue"))
            kegg_evalue = float(self.option("kegg_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范", code="15000107")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("cog_blast_evalue", cog_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围", code="15000108")
        if not self.option("cog_blast_evalue") > 0 and not self.option("cog_blast_evalue") < 1:
            raise OptionError("Cog比对的E值超出范围", code="15000109")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围", code="15000110")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围", code="15000111")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围", code="15000112")
        if self.option('align_species') != 'is_ncbi':
            if self.option('gff_or_gtf_file').prop['path'].endswith('gtf') and self.option('gff_or_gtf') == 'gff':
                raise OptionError("说好的你给我传gff，结果你给我gtf", code="15000113")
            if self.option('gff_or_gtf_file').prop['path'].endswith('gff') and self.option('gff_or_gtf') == 'gtf':
                raise OptionError("说好的你给我传gtf，结果你给我gff", code="15000114")
        if self.option('mapping_sample_percent') or self.option('rrna_sample_percent'):
            try:
                mapping_sample_percent = float(self.option("mapping_sample_percent"))
                mapping_ratio = float(self.option("mapping_ratio"))
                rrna_sample_percent = float(self.option("rrna_sample_percent"))
                rrna_ratio = float(self.option("rrna_ratio"))
            except:
                raise OptionError("判断终止那边不要传入非数字形式", code="15000115")
            else:
                if mapping_ratio and not mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整", code="15000116")
                if not mapping_ratio and mapping_sample_percent:
                    raise OptionError("判断终止mapping的参数不完整", code="15000117")
                if rrna_sample_percent and not rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整", code="15000118")
                if not rrna_sample_percent and rrna_ratio:
                    raise OptionError("判断终止rRNA的参数不完整", code="15000119")

        group_size = list()
        group_dict = self.option("group_table").prop['group_dict']
        for key in group_dict:
            group_size.append(len(group_dict[key]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if self.option("diff_method").lower == "deseq2":
                self.option("diff_method", "DEGseq")
                self.option("diff_fdr_ci", 0.001)
                self.logger.info("改项目没有生物学重复,不可以使用DESeq2,修改方法为DEGseq,阈值设置为0.001")
        elif group_size[0] == 1 and group_size[-1] >= 2:
            if self.option("diff_method").lower == "deseq2":
                self.option("diff_method", "edgeR")
                self.option("diff_fdr_ci", 0.05)
                self.logger.info("改项目部分组别有生物学重复,部分组别无生物学重复,不可以使用DESeq2,修改方法为edgeR,阈值设置为0.001")
        elif group_size[0] >= 2 and group_size[-1] >= 2:
            if self.option("diff_method").lower == "degseq":
                self.option("diff_method", "DESeq2")
                self.option("diff_fdr_ci", 0.05)
                self.logger.info("改项目有生物学重复,不可以使用DEGseq,DESeq2,阈值设置为0.001")
        else:
            pass
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.run_stage_1()
        super(ProkrnaWorkflow, self).run()

    def run_stage_1(self):
        self.logger.info("=====run stage one=====")
        self.rock_index.on('end', self.run_filecheck)
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_extract_biotype)  # 提取基因的biotype
        self.extract_biotype.on('end', self.run_extract_ptt)  # 过滤ptt文件
        self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.mapping.on('end', self.check_rrna_and_mapping)
        self.on_rely([self.qc_stat_after, self.qc_stat_before], self.run_mapping)
        if self.ncbi is False:
            if self.option('gff_or_gtf') == 'gbk':
                self.run_gbk2gff()
            else:
                self.run_check_genome()
        if self.ncbi is True:
            self.run_download4ncbi()

    def check_rrna_and_mapping(self):
        self.logger.info("=====check rrna and mapping=====")
        self.logger.info("=====check rrna=====")
        is_pass = True
        reason = ''
        if self.option('rrna_stop'):
            rrnastat_files = glob.glob(os.path.join(self.qc_stat_after.output_dir, 'RfamStat', '*.stat.*'))
            num = 0
            for file in rrnastat_files:
                with open(file, 'r') as f:
                    lines = f.readlines()
                    r_rna_ratio = float(lines[1].split("\t")[1])
                    if r_rna_ratio > float(self.option('rrna_ratio')):
                        num += 1
            rrna_sample_percent = float(num)/float(len(rrnastat_files)) * 100
            if rrna_sample_percent > float(self.option('rrna_sample_percent')):
                is_pass = False
                reason = 'rrna'
        self.logger.info("=====check mapping=====")
        mappingstat_files = glob.glob(os.path.join(self.mapping.output_dir, 'stat/*.stat'))
        num = 0
        flag = False
        for file in mappingstat_files:
            with open(file, 'r') as f:
                lines = f.readlines()
                if 'bam_sort_core' not in lines[-2]:
                    mapping_ratio = lines[-2].split(" ")[0].rstrip("%")
                else:
                    mapping_ratio = lines[-3].split(" ")[0].rstrip("%")
                if self.option('mapping_ratio') != '' and float(mapping_ratio) * 100 < float(self.option('mapping_ratio')):
                    num += 1
                if mapping_ratio < 10:
                    flag = True
        mapping_percent = float(num) / float(len(mappingstat_files)) * 100
        if (self.option('mapping_stop') and mapping_percent > float(self.option('mapping_sample_percent'))) or flag:
            is_pass = False
            if reason is not None:
                reason = 'both'
            else:
                reason = 'mapping'
        if not is_pass:
            self.stop(reason)
        else:
            self.run_stage_2()

    def stop(self, reason):
        assert reason in ['rrna', 'mapping', 'both']
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_genome_info()
        self.export_rock_index()
        self.export_qc()
        self.export_productive_name()
        self.export_mapping_only()
        if reason == 'rrna':
            msg = 'Workflow will stop because the rRNA ratio is not up to par'
            self.logger.warn(msg)
        elif reason == 'mapping':
            msg = 'Workflow will stop because the alignment rate is not up to par'
            self.logger.warn(msg)
        else:
            msg = 'Workflow will stop because the alignment rate and rRNA ratio both are not up to par'
            self.logger.warn(msg)
        receiver = ['caiping.shi@majorbio.com', 'meta@majorbio.com']
        a = SendEmail("897236887@qq.com", "smtp.qq.com", "fhwuvcclstjqbfga", "897236887@qq.com", ','.join(receiver),
                      "WARNING - project_sn ({}), task_id ({})".format(self._sheet.project_sn, self._sheet.id), 465)
        a.send_msg(msg)
        a.send_email()
        # db = Config().get_mongo_client(mtype='project', dydb_forbid=True)[Config().get_mongo_dbname(mtype='project', dydb_forbid=True)]
        # email = db['sg_task_email']
        # email.update({'task_id': self.task_id}, {'$set': {'status': '5'}}, upsert=True)
        super(ProkrnaWorkflow, self).end()

    def run_stage_2(self):
        self.logger.info("=====run stage two=====")
        self.rockhopper.on('end', self.run_extract_mrna)
        self.extract_mrna.on('end', self.run_extract_pep)
        self.extract_mrna.on('end', self.run_annot_mapdb)
        self.extract_pep.on('end', self.run_annot_orfpfam)
        self.extract_mrna.on('end', self.run_promote)
        self.extract_mrna.on('end', self.run_terminator)
        self.on_rely([self.annot_mapdb, self.annot_orfpfam], self.run_annot_filter)
        self.on_rely([self.annot_filter, self.extract_ptt], self.run_annot_class)
        if self.option("strand_specific"):
            self.on_rely([self.rockhopper, self.extract_biotype], self.run_extract_mrna_srna)
            # self.on_rely([self.rockhopper, self.extract_biotype], self.run_extract_mrna)
            self.extract_mrna_srna.on('end', self.run_express)
            self.extract_mrna_srna.on('end', self.run_express_total)
            self.run_rockhopper()
            self.run_srna()
        else:
            self.extract_mrna_srna.on('end', self.run_express)
            self.extract_mrna_srna.on('end', self.run_express_total)
            self.final_tools.remove(self.srna)
            self.run_extract_mrna_srna()
        self.express.on('end', self.run_diffexpress)
        self.diffexpress.on('end', self.run_exp_corr)
        if not self.option("group_table").prop["sample_number"] < 3:
            self.diffexpress.on('end', self.run_exp_pca)
            if self.min_group_num >= 3:
                self.exp_pca.on('end', self.run_ellipse)
        self.diffexpress.on('end', self.run_exp_venn)
        self.on_rely([self.annot_class,self.diffexpress],self.run_geneset_analysis)
        self.on_rely(self.final_tools, self.run_chart)
        self.run_map_assess()
        self.run_snp()

    def run_chart(self):
        '''
        绘图步骤插入在导表前
        '''
        self.chart = self.add_tool('prok_rna.chart')
        chart_dict = {
            "type": "workflow",
            "qc_file_before": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=self.qc_stat_before.output_dir, sample_name='{sample_name}'),
            "qc_file_after": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=self.qc_stat_after.output_dir, sample_name='{sample_name}'),
        }

        if self.option("group_table").is_set:
            # group_dict = self.option("group_table").prop["group_dict"]
            samples = self.option("group_table").prop["sample"]
            group_dict = OrderedDict()
            group_table = pd.read_table(self.option('group_table').prop['path'], header=0)
            for each in group_table.values.tolist():
                if each[1] not in group_dict.keys():
                    group_dict[each[1]] = list()
                group_dict[each[1]].append(each[0])

            chart_dict.update({
                "group_dict": group_dict,
                "samples": samples,
            })

        if self.option('fq_type') == "SE":
            chart_dict.update({
                "qc_file_before": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.qc_stat_before.output_dir, sample_name='{sample_name}'),
                "qc_file_after": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.qc_stat_after.output_dir, sample_name='{sample_name}')
            })

        analysis = self.option("map_assess_method").split(",")
        if "saturation" in analysis:
            chart_dict.update({
                'map_saturation': "{table_dir}/saturation/satur_{sample_name}.eRPKM.xls".format(table_dir=self.map_qc.output_dir, sample_name='{sample_name}')
            })
        if 'coverage' in analysis:
            chart_dict.update({
                'map_coverage': '{table_dir}/coverage/coverage_{sample_name}.geneBodyCoverage.txt'.format(table_dir=self.map_qc.output_dir, sample_name='{sample_name}'),
            })

        # annotation venn
        chart_dict.update({
            'annot_venn': '{table_dir}/anno_stat/venn'.format(table_dir=self.annot_class.output_dir),
        })

        # annotation stats
        chart_dict.update({
            'annot_stats': '{table_dir}/anno_stat/all_annotation_statistics.xls'.format(table_dir=self.annot_class.output_dir),
        })

        # cog annotation
        chart_dict.update({
            'cog_annot': '{table_dir}/cog/cog_summary.xls'.format(table_dir=self.annot_class.output_dir),
        })

        # go annotation
        chart_dict.update({
            'go_annot_level2': '{table_dir}/go/go_detail_stat.tsv'.format(table_dir=self.annot_class.output_dir),
        })

        # kegg annotation
        chart_dict.update({
            'kegg_annot': '{table_dir}/kegg/kegg_layer.xls'.format(table_dir=self.annot_class.output_dir),
        })

        # expression
        chart_dict.update({
            'exp_matrix': '{table_dir}/transcript.{exp_type}.matrix'.format(table_dir=self.express.output_dir, exp_type=self.option('exp_way').lower()),
            # 'group_table': self.option('group_table').prop['path'],
        })

        # expression venn
        if len(self.option("group_table").prop['group_dict']) > 1:
            chart_dict.update({
                "exp_venn": self.exp_venn.output_dir + "/venn_graph.xls"
            })

        # expression correlation heatmap & pca
        if self.option('group_table').prop['sample_number'] > 2:
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
        if self.ellipse in self.final_tools:
            chart_dict.update({
                'exp_pca_ellipse': "{table_dir}/ellipse_out.xls".format(table_dir=self.ellipse.output_dir),
            })

        # differential analysis
        cmp_list = list()
        if self.option("control_file").is_set:
            with open(self.option("control_file").prop['path'], 'r') as c:
                c.readline()
                for line in c:
                    if len(line.strip().split('\t')) > 1:
                        cmp_list.append(line.strip().split('\t'))
            # cmp_list = self.option("control_file").prop["cmp_list"]
            chart_dict.update({
                "cmp_list": cmp_list
            })

        chart_dict.update({
            'diff_method': self.option('diff_method'),
            'gene_diff_scatter': '{table_dir}/{cmp}.{diff_method}.xls'.format(table_dir=self.diffexpress.output_dir,
                                                                              cmp='{cmp1}_vs_{cmp2}',
                                                                              diff_method=self.option('diff_method').lower()),
            'gene_diff_summary': '{table_dir}/diff_summary_{diff_method}.xls'.format(table_dir=self.diffexpress.output_dir,
                                                                                     diff_method=self.option('diff_method').lower()),
        })

        # sRNA
        if self.option("strand_specific"):
            chart_dict.update({
                'srna_length': '{table_dir}/rockhopper/genome.predicted_RNA.bed.xls'.format(table_dir=self.srna.output_dir),
                'srna_venn': '{table_dir}/srna_annot/'.format(table_dir=self.srna.output_dir),
                'srna_rfam': '{table_dir}/srna_annot/rfam_stat.xls'.format(table_dir=self.srna.output_dir),
            })

        # operon
        chart_dict.update({
            'operon_xls': '{table_dir}/operon.xls'.format(table_dir=self.rockhopper.output_dir),
        })

        # UTR
        chart_dict.update({
            'utr_xls': '{table_dir}/UTR.xls'.format(table_dir=self.rockhopper.output_dir),
        })

        # SNP/InDel
        chart_dict.update({
            'snp_anno': '{table_dir}/snp_anno.xls'.format(table_dir=self.snp.output_dir),
        })

        # geneset pipeline
        geneset_outdir = self.diff_geneset_analysis.output_dir
        if os.path.exists(geneset_outdir) and not os.path.exists(os.path.join(geneset_outdir, 'results_info')):
            # venn
            chart_dict.update({
                'geneset_detail': '{table_dir}/Uniform/output/all_detail.txt'.format(table_dir=self.diffexpress.work_dir),
            })
            # cluster
            chart_dict.update({
                'cluster_exp': '{table_dir}/cluster/All_Diff_mRNA/expression_matrix.xls'.format(table_dir=geneset_outdir),
                'cluster_tree': '{table_dir}/cluster/All_Diff_mRNA/seq.cluster_tree.txt'.format(table_dir=geneset_outdir),
                'sample_tree': '{table_dir}/cluster/All_Diff_mRNA/sample.cluster_tree.txt'.format(table_dir=geneset_outdir),
                'subclusters': '{table_dir}/cluster/All_Diff_mRNA/seq.subcluster*xls'.format(table_dir=geneset_outdir),
            })
            # geneset annotation
            chart_dict.update({
                'analysis_json': '{table_dir}/All_Diff_mRNA/analysis_json'.format(table_dir=geneset_outdir),
                'kegg_level': '{table_dir}/DiffGenesetAll/AnnotPrepare/output/gene_kegg_level_table.xls'.format(table_dir=self.diff_geneset_analysis.work_dir)
            })

        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
            json.dump(chart_dict, json_f, indent=4)
        opts = {"file_json": self.work_dir + "/chart_workflow.json"}
        if self.option('report_img'):
            opts.update({"chart_report": True})
        self.chart.set_options(opts)
        self.chart.on('end', self.end)
        self.chart.run()

    def run_gbk2gff(self):
        self.gbk2gff = self.add_tool("prok_rna.gbk_to_gff")
        opts = {
            'gbk_file': self.option('gff_or_gtf_file').prop['path'],
        }
        self.gbk2gff.set_options(opts)
        self.gbk2gff.on('start', self.set_step, {'start': self.step.gbk2gff})
        self.gbk2gff.on('end', self.set_step, {'end': self.step.gbk2gff})
        self.gbk2gff.on('end', self.run_check_genome)
        self.gbk2gff.run()

    def run_check_genome(self):
        if self.option("align_species") == 'is_ncbi':
            if self.ncbi is False:
                opts = {
                    "genome": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                    "in_file": os.path.join(self.db_path, 'genome_gff', self.option('genome_id') + '.gff'),
                    "in_type": "gff"
                }
            else:
                ncbi_file = os.path.join(self.download4ncbi.work_dir, 'file.txt')
                with open(ncbi_file, 'r') as n:
                    genome_id, self.fna, self.gff, self.report = n.readline().strip().split('\t')
                if 'gff' in self.gff:
                    in_type = 'gff'
                else:
                    in_type = 'gtf'
                opts = {
                    "genome": self.fna,
                    "in_file": self.gff,
                    "in_type": in_type
                }
        else:
            if self.option('gff_or_gtf') == 'gbk':
                in_type = 'gff'
                in_file = self.gbk2gff.option('gff_file').prop['path']
            else:
                in_type = self.option('gff_or_gtf')
                in_file = self.option('gff_or_gtf_file').prop['path']

            opts = {
                "genome": self.option('genome_db').prop['path'],
                "in_file": in_file,
                "in_type": in_type
            }
        self.genomecheck.set_options(opts)
        self.genomecheck.on('start', self.set_step, {'start': self.step.genomecheck})
        self.genomecheck.on('end', self.set_step, {'end': self.step.genomecheck})
        self.genomecheck.on('end', self.run_rock_index)
        self.genomecheck.run()

    def run_download4ncbi(self):
        self.download4ncbi = self.add_tool("prok_rna.download4ncbi.download")
        opts = {
            "id": self.option('genome_id')
        }
        self.download4ncbi.set_options(opts)
        self.download4ncbi.on('end', self.set_output, 'download4ncbi')
        self.download4ncbi.on('start', self.set_step, {'start': self.step.download4ncbi})
        self.download4ncbi.on('end', self.set_step, {'end': self.step.download4ncbi})
        self.download4ncbi.on('end', self.run_check_genome)
        self.download4ncbi.run()

    def run_rock_index(self):
        gtf_gff = self.genomecheck.option("out_file").prop['path']
        if gtf_gff.endswith("gtf"):
            type = 'gtf'
        else:
            type = 'gff'
        if self.option("align_species") == 'is_ncbi':
            if self.ncbi is False:
                opts = {
                    "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                    "input_file": gtf_gff,
                    "type": type,
                }
                self.ref_genome = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna')
                self.db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
                col = self.db["sg_genome_db"]
                genome_info = col.find_one(
                    {"genome_id": self.option("genome_id")})
                self.species_name = genome_info['organism_name']
                self.fna = self.ref_genome
            else:
                opts = {
                    "fna": self.fna,
                    "input_file": gtf_gff,
                    "type": type
                }
                self.ref_genome = self.fna
                with open(self.report, 'r') as r:
                    lines = r.read()
                    self.species_name = re.findall(r'Organism name\:\s\s(.*?)\r', lines)[0]
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": gtf_gff,
                "type": type
            }
            self.ref_genome = self.option('genome_db').prop['path']
            self.species_name = self.option('species_name')
        self.rock_index.set_options(opts)
        self.rock_index.on('end', self.set_output, 'rock_index')
        self.rock_index.on('start', self.set_step, {'start': self.step.rock_index})
        self.rock_index.on('end', self.set_step, {'end': self.step.rock_index})
        self.rock_index.run()
        pass

    def run_extract_biotype(self):
        gtf_gff = self.genomecheck.option("out_file").prop['path']
        if gtf_gff.endswith("gtf"):
            type = 'gtf'
        else:
            type = 'gff'
        if type == "gff":
            opts = {
                "gff": gtf_gff,
            }
        else:
            opts = {
                "gtf": gtf_gff,
            }
        self.extract_biotype.set_options(opts)
        self.extract_biotype.run()

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            "in_gtf": self.rock_index.option('gtf_exon'),
            "ref_genome": self.ref_genome,
            'is_duplicate': self.option('is_duplicate')
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
        # self.fastp_rna = self.add_module('datasplit.fastp_rna')
        fq_dir = self.option('fastq_dir').path
        sample_path = os.path.join(fq_dir, 'abs.list.txt')
        open(sample_path, 'w').writelines(
            '{}/{}'.format(fq_dir, line) for line in open(os.path.join(fq_dir, 'list.txt'))
        )
        if self.option('fq_type') == "PE":
            self.qc.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
            })
        else:
            self.qc.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
            })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.rna_qc})
        self.qc.on('end', self.set_step, {'end': self.step.rna_qc})
        self.qc.run()

    def run_qc_stat(self, event):
        if self.option('quality_score_system').lower() == "phred+33" or self.option(
                'quality_score_system').lower() == "phred 33":
            quality = 33
        else:
            quality = 64
        if event['data']:
            options = {
                'fastq_dir': self.qc.option('sickle_dir'),
                'fq_type': self.option('fq_type'),
                'quality': quality,
                'rfam': True,
                'dup': True
            }
            if os.path.exists(self.rock_index.output_dir + "/rrna.fa"):
                options.update({
                    'rrna': self.rock_index.output_dir + "/rrna.fa"
                })
            self.qc_stat_after.set_options(options)
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
            "seq_method": self.option("fq_type"),  # PE or SE
            "fastq_dir": self.qc.option("sickle_dir"),
            "ref_gtf": self.rock_index.option("gtf"),
        }
        if self.option("strand_specific"):
            opts.update(
                {
                    "strand_specific": True,
                }
            )
        self.mapping.set_options(opts)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_rockhopper(self):
        self.logger.info("开始运行rockhopper")
        gtf_gff = self.genomecheck.option("out_file").prop['path']
        if gtf_gff.endswith("gtf"):
            type = 'gtf'
        else:
            type = 'gff'
        if self.option("align_species") == 'is_ncbi':
            if self.ncbi is False:
                opts = {
                    "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                    "input_file": gtf_gff,
                    "type": type
                }
            if self.ncbi is True:
                opts = {
                    "fna": self.fna,
                    "input_file": gtf_gff,
                    "type": type
                }
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": gtf_gff,
                "type": type
            }
        opts.update({"group_list": self.option('group_table').prop['path'],
                     "trimPairFq": self.qc.option('fq_list').prop['path']})
        self.rockhopper.set_options(opts)
        self.rockhopper.on("end", self.set_output, "rockhopper")
        self.rockhopper.on('start', self.set_step, {'start': self.step.rockhopper})
        self.rockhopper.on('end', self.set_step, {'end': self.step.rockhopper})
        self.rockhopper.run()

    def run_annot_mapdb(self, event):
        self.logger.info("开始运行diamond注释")
        opts = {
            "query": self.extract_mrna.option("out_fasta"),
            "method": "diamond",
            "nr_db": "bacteria",
            'kegg_version': self.annot_config_dict['kegg']['version'],
            'nr_version': self.annot_config_dict['nr']['version'],
            'version': self.annot_config_dict['diamond_all']['version'],
            'eggnog_version': self.annot_config_dict['eggnog']['version'],
            'string_version': self.annot_config_dict['string']['version'],
            "swissprot_version": self.annot_config_dict['swissprot']['version'],
            "pir_version": self.annot_config_dict['pir']['version'],
            "lines": 5000,
        }

        if self.option("annot_group") == "REFRNA_GROUP_202110":
            opts.update({
                "database": "nr;kegg;swissprot;cog",
                'cog_version': self.annot_config_dict['cog']['version'],
                'diamond_version' : "v2.0.13"
            })
        self.annot_mapdb.set_options(opts)
        self.annot_mapdb.on("end", self.set_output, "annot_mapdb")
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.logger.info("开始运行pfam注释")
        opts = {
            "lines": 5000,
            "pep": self.extract_pep.option("out_fasta"),
            "pfam_version": self.annot_config_dict['pfam']['version'],
        }
        self.annot_orfpfam.set_options(opts)
        self.annot_orfpfam.on("end", self.set_output, "annot_orfpfam")
        self.annot_orfpfam.on('start', self.set_step, {'start': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_step, {'end': self.step.annot_orfpfam})
        self.annot_orfpfam.run()

    def run_annot_filter(self):
        options = {
            "blast_nr_xml": self.annot_mapdb.output_dir + "/nr/blast.xml",
            "blast_eggnog_xml": self.annot_mapdb.output_dir + "/eggnog/blast.xml",
            "blast_kegg_xml": self.annot_mapdb.output_dir + "/kegg/blast.xml",
            "blast_swissprot_xml": self.annot_mapdb.output_dir + "/swissprot/blast.xml",
            "pfam_domain": self.annot_orfpfam.output_dir + "/pfam_domain",
            "blast2go_annot": self.annot_mapdb.output_dir + "/GO/go_annot.xls",
            'nr_evalue': self.option('nr_blast_evalue'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        } 
        if self.option("annot_group") == "REFRNA_GROUP_202110":
            options.update({"blast_eggnog_xml": self.annot_mapdb.output_dir + "/cog/blast.xml"})

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
            
            'blast_swissprot_xml': filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': filter_dir + "/pfam/pfam_domain.filter.xls",
            "kegg_version": self.annot_config_dict['kegg']['version'],
            "nr_version": self.annot_config_dict['nr']['version'],
            "eggnog_version": self.annot_config_dict['eggnog']['version'],
            "string_version": self.annot_config_dict['string']['version'],
            "go_version": self.annot_config_dict['go']['version'],
            "pir_version": self.annot_config_dict['pir']['version'],
            "swissprot_version": self.annot_config_dict['swissprot']['version'],
            "blast2go_annot": filter_dir + "/go/go_annot.xls.filter.xls",
            "des": self.extract_ptt.option('out_ptt').prop['path'],
            "novel": self.rockhopper.output_dir + "/genome.predicted_cds.bed.xls"
        }
        if self.option("annot_group") == "REFRNA_GROUP_202110":
            options.update({
                'blast_ncbicog_xml': filter_dir + "/eggnog/blast.xml.filter.xml",
                'cog_verion': self.annot_config_dict['cog']['version']
            })
        else:
            options.update({
                'blast_eggnog_xml': filter_dir + "/eggnog/blast.xml.filter.xml"
            })
            
        self.annot_class.set_options(options)
        self.annot_class.on('start', self.set_step, {'start': self.step.annot_class})
        self.annot_class.on('end', self.set_step, {'end': self.step.annot_class})
        self.annot_class.on('end', self.set_output, "annot_class")
        self.annot_class.run()

    def run_extract_mrna_srna(self):
        self.logger.info("开始提取mRNA+sRNA序列")
        if self.option("strand_specific") == True:
            opts = {
                "fasta": self.rockhopper.option("feature_fa").prop['path'],
                "predicted_cds": self.rockhopper.output_dir + "/genome.predicted_cds.fa",
                "biotype": self.extract_biotype.option("gene_biotype"),
                "rna_type": "mRNA+sRNA",
            }
        else:
            opts = {
                "predicted_cds": self.rockhopper.output_dir + "/genome.predicted_cds.fa",
                "fasta": self.rock_index.option("query").prop["path"],
                "biotype": self.extract_biotype.option("gene_biotype"),
                "rna_type": "mRNA+sRNA",
            }
        self.extract_mrna_srna.set_options(opts)
        self.extract_mrna_srna.run()

    def run_extract_ptt(self):
        self.logger.info("开始过滤ptt文件")
        opts = {
            "ptt": self.rock_index.option('ptt').prop['path'],
            "biotype": self.extract_biotype.option("gene_biotype"),
            "rna_type": "mRNA+sRNA",
        }
        self.extract_ptt.set_options(opts)
        self.extract_ptt.run()

    def run_extract_mrna(self):
        self.logger.info("开始提取mRNA序列")
        opts = {
            "fasta": self.rock_index.option("query"),
            "predicted_cds": self.rockhopper.output_dir + "/genome.predicted_cds.fa",
            "biotype": self.extract_biotype.option("gene_biotype"),
            "rna_type": "mRNA",
        }
        self.extract_mrna.set_options(opts)
        self.extract_mrna.run()

    def run_extract_pep(self):
        self.logger.info("开始提取蛋白序列")
        opts = {
            "fasta": self.rock_index.option("querypep"),
            "biotype": self.extract_biotype.option("gene_biotype"),
            "rna_type": "mRNA",
        }
        self.extract_pep.set_options(opts)
        self.extract_pep.run()

    def run_express(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option("fq_type") == "PE":
                self.libtype = "rf"
            else:
                self.libtype = "r"
        else:
            self.libtype = None
        opts = {
            "fastq": self.qc.option("fq_list").path,
            "method": self.option("express_method").lower(),
            "libtype": self.libtype,
            "transcriptome": self.extract_mrna_srna.option("out_fasta"),
            "id2name": self.rock_index.option("ptt").prop["path"],
            "biotype": self.extract_biotype.option("gene_biotype"),
        }
        self.express.set_options(opts)
        self.express.on("end", self.set_output, "express")
        self.express.on('start', self.set_step, {'start': self.step.express})
        self.express.on('end', self.set_step, {'end': self.step.express})
        self.express.run()

    def run_express_total(self):
        self.logger.info("开始运行表达定量分析")
        if self.option("strand_specific") == True:
            if self.option("fq_type") == "PE":
                self.libtype = "rf"
            else:
                self.libtype = "r"
            transcriptome = self.rockhopper.option("feature_fa").prop['path']
        else:
            self.libtype = None
            transcriptome = self.rock_index.option("query")
        opts = {
            "fastq": self.qc.option("fq_list").path,
            "method": self.option("express_method").lower(),
            "libtype": self.libtype,
            "transcriptome": transcriptome,
            "id2name": self.rock_index.option("ptt").prop["path"],
            "biotype": self.extract_biotype.option("gene_biotype"),
        }
        self.express_total.set_options(opts)
        self.express_total.run()

    def deal_express_file(self):
        self.express_deal = self.express.output_dir + "/transcript.{}.matrix.dealed".format(
            self.option("exp_way").lower())
        exp_df = pd.read_table(self.express.output_dir + "/ref_gene.{}.matrix".format(self.option("exp_way").lower()),
                               header=0)
        exp_df.drop(['gene_name', 'type', 'biotype'], axis=1, inplace=True)
        exp_df.to_csv(self.express_deal, sep="\t", index=False)

    def run_exp_pca(self):
        # time.sleep(60)
        self.logger.info("开始运行pca")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp": self.express_deal
            }
        else:
            opts = {
                "exp": self.express_deal
            }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.run()

    def run_ellipse(self):
        self.ellipse.set_options({
            'analysis': 'pca',
            'group_table': self.option('group_table').prop['path'],
            'pc_table': os.path.join(self.exp_pca.output_dir, 'PCA.xls'),
        })
        self.ellipse.on('start', self.set_step, {'start': self.step.ellipse})
        self.ellipse.on('end', self.set_step, {'end': self.step.ellipse})
        self.ellipse.on('end', self.set_output, 'exp_pca')
        self.ellipse.run()

    def run_exp_corr(self):
        # time.sleep(60)
        self.logger.info("开始运行聚类分析")
        if self.option("express_method").lower() == "rsem" and self.option("exp_way").lower() == "fpkm":
            opts = {
                "exp": self.express_deal
            }
        else:
            opts = {
                "exp": self.express_deal
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
                "express_matrix": self.express_deal,
                "group_table": self.option('group_table')
            }
        else:
            opts = {
                "express_matrix": self.express_deal,
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
            count_file = self.express.output_dir + "/ref_gene.count.matrix"
            fpkm_file = self.express_deal
            # opts = {
            #     "count": count_file,
            #     "exp": fpkm_file,
            #     "group": self.option("group_table"),
            #     "cmp": self.option("control_file"),
            #     "pvalue_padjust": self.option("pvalue_padjust"),
            #     "pvalue": float(self.option("diff_fdr_ci")),
            #     "fc": float(self.option("fc")),
            #     "padjust_way": self.option("padjust_way"),
            #     "method": self.option("diff_method")
            # }
            opts = {
                "is_workflow": "yes",
                "count": count_file,
                "exp_matrix": fpkm_file,
                "group": self.option("group_table").path,
                "cmp": self.option("control_file").path,
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
                    self.option('prob'): float(self.option('diff_fdr_ci')),
                })
        else:
            count_file = self.express.output_dir + "/ref_gene.count.matrix"
            tpm_file = self.express_deal
            # opts = {
            #     "count": count_file,
            #     "exp": tpm_file,
            #     "group": self.option("group_table"),
            #     "cmp": self.option("control_file"),
            #     "pvalue_padjust": self.option("pvalue_padjust"),
            #     "pvalue": float(self.option("diff_fdr_ci")),
            #     "fc": float(self.option("fc")),
            #     "padjust_way": self.option("padjust_way"),
            #     "method": self.option("diff_method")
            # }
            opts = {
                "is_workflow": "yes",
                "count": count_file,
                "exp_matrix": tpm_file,
                "group": self.option("group_table").path,
                "cmp": self.option("control_file").path,
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
                    self.option('prob'): float(self.option('diff_fdr_ci')),
                })
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

    def run_snp(self):
        self.logger.info("开始运行snp步骤")
        # opts = {
        #     "ref_genome_custom": os.path.join(self.db_path, 'genome_fna',
        #                                       self.option('genome_id') + '.fna') if self.option(
        #         'align_species') == 'is_ncbi' else self.option('genome_db').prop['path'],
        #     "ref_genome": "customer_mode",
        #     "ref_gtf": self.rock_index.option('gtf'),
        #     "fq_list": self.qc.option("fq_list").prop["path"],
        #     "bamlist": self.mapping.option("bamlist"),
        #     "id2name": self.rock_index.option("ptt").prop['path'],
        # }
        opts = {
            "ref_genome_custom": self.ref_genome,
            "ref_genome": "customer_mode",
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
        gtf_gff = self.genomecheck.option("out_file").prop['path']
        if gtf_gff.endswith("gtf"):
            type = 'gtf'
        else:
            type = 'gff'
        if self.option("align_species") == 'is_ncbi':
            if self.ncbi is False:
                opts = {
                    "fna": os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna'),
                    "input_file": gtf_gff,
                    "type": type
                }
            if self.ncbi is True:
                opts = {
                    "fna": self.fna,
                    "input_file": gtf_gff,
                    "type": type
                }
        else:
            opts = {
                "fna": self.option('genome_db').prop['path'],
                "input_file": gtf_gff,
                "type": type
            }
        opts.update({"group_list": self.option('group_table').prop['path'],
                     "trimPairFq": self.qc.option('fq_list').prop['path']})
        self.srna.set_options(opts)
        self.srna.on("end", self.set_output, "srna")
        self.srna.on('start', self.set_step, {'start': self.step.srna})
        self.srna.on('end', self.set_step, {'end': self.step.srna})
        self.srna.run()

    def run_promote(self):
        self.logger.info("开始运行promote")
        opts = {
            "sequence": self.extract_mrna.option('out_fasta'),
            "assemble": self.fna if self.option(
                'align_species') == 'is_ncbi' else self.option('genome_db').prop['path'],
            "rock_index": self.rock_index.option('ptt'),
            "sample": "all"
        }
        self.promote.set_options(opts)
        self.promote.on("end", self.set_output, "promote")
        self.promote.on('start', self.set_step, {'start': self.step.promote})
        self.promote.on('end', self.set_step, {'end': self.step.promote})
        self.promote.run()

    def run_terminator(self):
        self.logger.info("开始运行terminator")
        opts = {
            "rock_index": os.path.join(self.rock_index.work_dir, 'rock_index'),
        }
        self.terminator.set_options(opts)
        self.terminator.on("end", self.set_output, "terminator")
        self.terminator.on('start', self.set_step, {'start': self.step.terminator})
        self.terminator.on('end', self.set_step, {'end': self.step.terminator})
        self.terminator.run()

    def run_geneset_analysis(self):
        diff_option_file = os.path.join(self.diffexpress.work_dir, "option_file")
        with open(diff_option_file, "r") as r:
            diff_options = json.load(r)
        self.option("diff_method", diff_options["method"])
        if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            self.option('pvalue_padjust', diff_options["pvalue_padjust"])
            self.option("fc", diff_options["fc"])
        else:
            self.option("diff_fdr_ci",diff_options["prob"])

        diff_dir = self.diffexpress.output_dir
        diff_dir = self.diffexpress.uniform.output_dir
        exp = os.path.join(self.express.output_dir, 'transcript.tpm.matrix')
        annot_result = self.annot_class.output_dir
        group_file = self.option("group_table").prop["path"]
        level = "T"
        kegg_version = self.annot_config_dict["kegg"]["version"]
        species = self.species_name
        opts = {
            'diff_path': diff_dir,
            'annot_result': annot_result,
            'diff_method': self.option('diff_method'),
            'transcript_exp_file': exp,
            'group': group_file,
            'level': level,
            "kegg_version": kegg_version,
            "go_version": self.annot_config_dict['go']['version'],
            'species': species,
        }
        self.diff_geneset_analysis.set_options(opts)
        self.diff_geneset_analysis.on('start', self.set_step, {'start': self.step.geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_step, {'end': self.step.geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_output, 'genesets_analysis')
        self.diff_geneset_analysis.run()


    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error("需要移动到output目录的文件夹不存在", code="15000120")
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
        if event['data'] == 'terminator':
            self.move2outputdir(obj.output_dir, 'terminator')
        if event['data'] == 'annot_class':
            self.move2outputdir(obj.output_dir, 'annot_class')
            # self.end_number += 1
        if event['data'] == 'rockhopper':
            self.move2outputdir(obj.output_dir, 'rockhopper')
        if event['data'] == 'rock_index':
            self.move2outputdir(obj.output_dir, 'rock_index')

    def end(self):
        intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        self.run_api()
        ## 更新一系列主表的字段，用于页面交互分析
        self.fa = os.path.join(intermediate_dir, "genome.feature.fa") if self.option(
            'strand_specific') else os.path.join(self.workflow_output, "Sequence_database/cds.fa")
        db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        col = db["sg_task"]
        col.update({"task_id": self.task_id}, {"$set": {"fastq": os.path.join(intermediate_dir, "fq_list.txt")}}, upsert=True)
        col.update({"task_id": self.task_id},
                   {"$set": {"rock_index": os.path.join(intermediate_dir, "Sequence_database/")}}, upsert=True)
        col.update({"task_id": self.task_id},
                   {"$set": {"ref_gtf": os.path.join(self.workflow_output, "Sequence_database/reshape.gtf")}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"ref_genome": self.ref_genome}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"strand_specific": self.option('strand_specific')}},
                   upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"genome_id": self.option('genome_id')}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"assemble_fa": self.fa}}, upsert=True)
        # col.update({"task_id": self.task_id},
        #            {"$set": {"fastq": self.workflow_output + "/QC/cleandata/fq_list.txt"}}, upsert=True)
        col.update({"task_id": self.task_id},
                   {"$set": {"id2name": os.path.join(intermediate_dir, "ptt.bed")}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"srna_files": os.path.join(self.workflow_output2, "sRNA/")}},
                   upsert=True)
        col.update({'task_id': self.task_id},
                   {'$set': {'database_version': {"kegg": self.option("kegg_version")}}}, upsert=True)
        
        annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        if self.annot_config_dict['kegg']['version'] > "2020":
            pass
            # if self.option("kegg_org") not in [None, ""]:
            #     annot_version_dict['kegg'] += "_spe"
        else:
            del annot_version_dict['kegg']
        col.update({'task_id': self.task_id},
                    {'$set': {'database_version': annot_version_dict,
                              'annot_group': self.option("annot_group")}}, upsert=True)

        count_matrix = os.path.join(self.express_total.output_dir, 'transcript.count.matrix')
        biotype = self.biotype_show(count_matrix)
        col.update(
            {'task_id': self.task_id},
            {'$set': {'biotype': biotype}}
        )

        col1 = db["sg_annotation_stat"]
        col1.update({"task_id": self.task_id}, {"$set": {"result_dir": os.path.join(intermediate_dir, "Annotation")}},
                    upsert=True)

        col2 = db["sg_exp"]
        if self.option("express_method").lower() == "rsem":
            col2.update({"task_id": self.task_id}, {
                "$set": {"count_file": os.path.join(self.workflow_output, "Express/ExpAnnalysis/transcript.count.matrix.xls"),
                         "exp_level": self.exp_level, "result_dir": os.path.join(intermediate_dir, "Express")}}, upsert=True)
        else:
            col2.update({"task_id": self.task_id}, {
                "$set": {"count_file": os.path.join(self.workflow_output, "Express/ExpAnnalysis/transcript.count.matrix.xls"),
                         "exp_level": self.exp_level}}, upsert=True)

        col3 = db["sg_specimen"]
        # col3.update({"task_id": self.task_id}, {"$set": {"fq_type": self.option('fq_type').lower()}}, upsert=True)
        col4 = db["sg_srna_fold"]
        col4.update({"task_id": self.task_id}, {"$set": {"result_dir": os.path.join(self.workflow_output2, "sRNA/srna_fold")}},
                    upsert=True)
        col5 = db["sg_snp"]
        col5.update({"task_id": self.task_id},
                    {"$set": {"result_dir": os.path.join(self.workflow_output, "SNP/snp_annotation_detail.xls")}},
                    upsert=True)
        self.merge_annotation_exp_matrix()  # 表达量表增加注释信息
        self.merge_annotation_exp_matrix_total()  # 表达量表增加注释信息
        self.merge_annotation_diffexp_matrix()  # 差异表达量表增加注释信息
        if self.option("report_img"):
            self.export_report_img()
        self.modify_output()
        super(ProkrnaWorkflow, self).end()

    def biotype_show(self, count_matrix):
        count_df = pd.read_table(count_matrix, header=0)
        biotype_sets = set(count_df['biotype'])
        if(len(biotype_sets) >1):
            return True
        else:
            return False

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def modify_output(self):
        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"

        ## 上传用于交互分析的中间结果文件
        self.intermediate_result = os.path.join(self.work_dir, 'intermediate_results')
        if os.path.isdir(self.intermediate_result):
            shutil.rmtree(self.intermediate_result)
        os.mkdir(self.intermediate_result)
        # AlignBam
        os.makedirs(os.path.join(self.intermediate_result, "Align/AlignBam"))
        for file in os.listdir(os.path.join(origin_dir, "mapping/bam")):
            file_path = os.path.join(origin_dir, "mapping/bam", file)
            os.link(file_path, os.path.join(self.intermediate_result, "Align/AlignBam", file))
        # modify fq_list.txt
        with open(origin_dir + "/QC_stat/fastq/fq_list.txt", "r") as f, open(
                self.intermediate_result + "/fq_list.txt", "w") as w:
            for line in f:
                items = line.strip().split("\t")
                if len(items) == 3:
                    w.write(items[0] + "\t" + os.path.join(self.workflow_output, os.path.basename(
                        items[1])) + "\t" + os.path.join(self.workflow_output, os.path.basename(items[2])) + "\n")
                else:
                    w.write(items[0] + "\t" + os.path.join(self.workflow_output, os.path.basename(items[1])) + "\n")
        # genome.feature.fa
        if self.option("strand_specific"):
            genome_feature_fa = os.path.join(self.rockhopper.output_dir, "genome.feature.fa")
            os.link(genome_feature_fa, os.path.join(self.intermediate_result, 'genome.feature.fa'))
        # seq_db
        os.mkdir(self.intermediate_result + "/Sequence_database")
        for file in os.listdir(os.path.join(origin_dir, 'rock_index')):
            os.link(os.path.join(origin_dir, 'rock_index', file),
                    os.path.join(self.intermediate_result, 'Sequence_database', file))
        # ptt.bed
        ptt_bed = os.path.join(origin_dir, 'rock_index', 'ptt.bed')
        if os.path.exists(ptt_bed):
            os.link(ptt_bed, os.path.join(self.intermediate_result, 'ptt.bed'))
        # Annotation
        os.makedirs(self.intermediate_result + "/Annotation")
        os.makedirs(self.intermediate_result + "/Annotation/go")
        os.makedirs(self.intermediate_result + "/Annotation/cog")
        os.makedirs(self.intermediate_result + "/Annotation/nr")
        os.makedirs(self.intermediate_result + "/Annotation/swissprot")
        os.makedirs(self.intermediate_result + "/Annotation/pfam")
        os.makedirs(self.intermediate_result + "/Annotation/kegg")
        os.makedirs(self.intermediate_result + "/Annotation/summary")
        genome_stat = self.genome_stat
        if os.path.exists(self.intermediate_result + "/Annotation/" + os.path.basename(genome_stat)):
            os.remove(self.intermediate_result + "/Annotation/" + os.path.basename(genome_stat))
        os.link(genome_stat, self.intermediate_result + "/Annotation/" + os.path.basename(genome_stat))
        for file in glob.glob(origin_dir + "/annot_class/cog/*"):
            os.link(file, self.intermediate_result + "/Annotation/cog/" + os.path.basename(file))
        for file in glob.glob(origin_dir + "/annot_class/go/*"):
            os.link(file, self.intermediate_result + "/Annotation/go/" + os.path.basename(file))
        for file in glob.glob(origin_dir + "/annot_class/kegg/*"):
            CopyFile().linkdir(file, self.intermediate_result + "/Annotation/kegg/" + os.path.basename(file))
        for file in os.listdir(origin_dir + "/annot_class/anno_stat/venn"):
            if file:
                dir = file.split('_')[0]
                os.link(origin_dir + "/annot_class/anno_stat/venn/" + file,
                        self.intermediate_result + "/Annotation/" + dir + '/' + file)
        for file in os.listdir(origin_dir + "/annot_class/anno_stat/blast"):
            if file:
                dir = file.split('.')[0]
                os.link(origin_dir + "/annot_class/anno_stat/blast/" + file,
                        self.intermediate_result + "/Annotation/" + dir + '/' + file)
        if os.path.exists(origin_dir + "/annot_orfpfam/pfam_domain"):
            os.link(origin_dir + "/annot_orfpfam/pfam_domain",
                    self.intermediate_result + "/Annotation/pfam/pfam_domain")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/all_anno_detail.xls"):
            os.link(origin_dir + "/annot_class/anno_stat/all_anno_detail.xls",
                    self.intermediate_result + "/Annotation/summary/all_anno_detail.xls")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/all_annotation_statistics.xls"):
            os.link(origin_dir + "/annot_class/anno_stat/all_annotation_statistics.xls",
                    self.intermediate_result + "/Annotation/summary/all_annotation_statistics.xls")
        if os.path.exists(origin_dir + "/annot_class/ref.txt"):
            os.link(origin_dir + "/annot_class/ref.txt", self.intermediate_result + "/Annotation/ref.txt")
        CopyFile().linkdir(origin_dir + "/annot_mapdb", self.intermediate_result + "/Annotation/annot_mapdb")
        CopyFile().linkdir(origin_dir + "/annot_orfpfam", self.intermediate_result + "/Annotation/annot_mapdb/pfam")
        # Express
        if self.option("express_method").lower() == "rsem":
            os.makedirs(self.intermediate_result + "/Express")
            for file in glob.glob(self.express.output_dir + "/*"):
                os.link(file, os.path.join(self.intermediate_result + "/Express/" + os.path.basename(file)))

        # Background
        os.mkdir(target_dir + "/Background")
        ## genome_stat
        genome_stat = self.annot_class.output_dir + "/genome_stat.xls"
        if os.path.exists(genome_stat):
            os.link(genome_stat, target_dir + "/Background/" + os.path.basename(genome_stat))
        ## software_info
        software_info = os.path.join(target_dir, 'Background', "software_info.xls")
        db = Config().get_mongo_client(mtype='prok_rna')[Config().get_mongo_dbname('prok_rna')]
        my_collection = db['sg_software_database']
        my_results = my_collection.find({})
        with open(software_info, "w") as w:
            w.write("\t".join(["Soft/Database", "Version", "Analysis", "Source"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["software_database"]), str(collection["version"]), str(collection["usage"]),
                     str(collection["source"])]) + "\n")
        ## sample_info
        sample_info = os.path.join(target_dir, 'Background', "sample_info.xls")
        productive_names = dict()
        if self.option("productive_table").is_set:
            with open(self.option("productive_table").path, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    items = line.strip().split("\t")
                    if len(items) >= 2:
                        productive_names[items[0]] = items[1]
        my_collection = db['sg_specimen']
        my_results = my_collection.find({"task_id": self.task_id})
        with open(sample_info, "w") as w:
            if self.option("productive_table").is_set:
                w.write("\t".join(
                    ["Sample Productive Name", "Sample Initial Name", "Sample Analysis Name", "Library", "Read Length (nt)", "Raw Data File"]) + "\n")
            else:
                w.write("\t".join(
                    ["Sample Initial Name", "Sample Analysis Name", "Library", "Read Length (nt)", "Raw Data File"]) + "\n")
            for collection in my_results:
                if str(collection["old_name"]) in productive_names:
                    w.write("\t".join([str(productive_names[collection["old_name"]]), str(collection["old_name"]),
                            str(collection["new_name"]), str(collection["fq_type"]), str(collection["read_length"]),
                            str(collection["raw_file"])]) + "\n")
                else:
                    w.write("\t".join([str(collection["old_name"]), str(collection["new_name"]), str(collection["fq_type"]),
                            str(collection["read_length"]), str(collection["raw_file"])]) + "\n")

        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        os.mkdir(target_dir + "/QC")
        if os.path.exists(fq_stat_before):
            os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        if os.path.exists(fq_stat_after):
            os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        rfam_stat = target_dir + "/QC/rrna_assessment.xls"
        rrna_stat = glob.glob(self.qc_stat_after.output_dir + '/RfamStat/' + "*_vs_rfam.stat.xls")
        with open(rfam_stat, "w") as w:
            w.write("Sample" + "\t" + "rRNA(%)" + "\n")
            for fs in sorted(rrna_stat):
                f = open(fs, "r")
                f.readline()
                items = f.readline().strip().split("\t")
                specimen_name = items[0]
                r_rna_ratio = items[1]
                w.write(specimen_name + "\t" + r_rna_ratio + "\n")

        # QC chart
        qc_pdfs = glob.glob(os.path.join(self.chart.work_dir, '*_qc_base.line.pdf'))
        qc_pdfs += glob.glob(os.path.join(self.chart.work_dir, '*_qc_error.line.pdf'))
        qc_pdfs += glob.glob(os.path.join(self.chart.work_dir, '*_qc_qual.box.pdf'))
        for each in qc_pdfs:
            s_name, type_info = os.path.basename(each).split('.', 1)
            new_name = s_name + '_' + type_info.split('_', 1)[0]
            if each.endswith('_qc_base.line.pdf'):
                new_name += '_bases_distribution.pdf'
            elif each.endswith('_qc_error.line.pdf'):
                new_name += '_bases_error_rate.pdf'
            elif each.endswith('_qc_qual.box.pdf'):
                new_name += '_bases_quality.pdf'
            new_path = os.path.join(target_dir, 'QC', new_name)
            self.move_pdf(each, new_path)

        # Align
        os.mkdir(target_dir + '/Align')
        mapping_stat_file = target_dir + '/Align/align_stat.xls'
        my_collection = db['sg_specimen_mapping']
        my_results = my_collection.find({"task_id": self.task_id})
        with open(mapping_stat_file, "w") as w:
            w.write("\t".join(
                ["Sample Name", "Total Reads", "Genome Mapped Reads", "Genome Mapped Ratio(%)", "Unmapped Reads",
                 "Unmapped Reads Ratio(%)", "Uniq Mapped Reads", "Uniq Mapped Reads Ratio(%)", "CDS Mapped Reads",
                 "CDS Mapped Ratio(%)"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["specimen_name"]), str(collection["total_reads"]), str(collection["mapping_reads"]),
                     str(collection["mapping_reads_ratio"]), str(collection["unmapping_reads"]),
                     str(collection["unmapping_reads_ratio"]), str(collection["uniq_mapped"]),
                     str(collection["uniq_mapped_ratio"]), str(collection.get("cds", "")),
                     str(collection["cds_ratio"])]) + "\n")
        ## QualityAssessment
        os.makedirs(target_dir + "/Align/QualityAssessment")
        chr_stat = target_dir + "/Align/QualityAssessment/chr_distribution.xls"
        chr_stat_data = dict()
        samples = list()
        chrs = list()
        for file in sorted(glob.glob(self.map_qc.output_dir + "/chr_stat/" + "*bam_chr_stat.xls")):
            sample_name = os.path.basename(file).split(".bam_chr_stat")[0]
            samples.append(sample_name)
            with open(file, "r") as f:
                f.readline()
                chr_stat_data[sample_name] = dict()
                for line in f:
                    item = line.strip().split("\t")
                    if item[0] not in chrs:
                        chrs.append(item[0])
                    if item[0] not in chr_stat_data[sample_name]:
                        chr_stat_data[sample_name][item[0]] = item[1]
        with open(chr_stat, "w") as w:
            w.write("Chromosome" + "\t" + "\t".join(samples) + "\n")
            for chr in chrs:
                w.write(chr)
                for sample in samples:
                    if chr not in chr_stat_data[sample]:
                        chr_stat_data[sample][chr] = 0
                    w.write("\t" + str(chr_stat_data[sample][chr]))
                w.write("\n")

        # Align Chart
        align_pdfs = glob.glob(os.path.join(self.chart.work_dir, '*map_saturation.line.pdf'))
        for each in align_pdfs:
            s_name = os.path.basename(each).split('.')[0]
            new_path = os.path.join(target_dir, 'Align', s_name + '_saturation_curve.pdf')
            self.move_pdf(each, new_path)

        coverage_pdf = os.path.join(self.chart.work_dir, 'map_coverage.line.pdf')
        new_path = os.path.join(target_dir, 'Align', 'coverage_distribution.pdf')
        self.move_pdf(coverage_pdf, new_path)


        # Sequence_database
        os.mkdir(target_dir + "/Sequence_database")
        cds_fa = origin_dir + "/rock_index/cds.fa"
        cds_faa = origin_dir + "/rock_index/cds.faa"
        reshape_gtf = origin_dir + "/rock_index/reshape.gtf"
        if os.path.exists(cds_fa):
            os.link(cds_fa, target_dir + "/Sequence_database/cds.fa")
        if os.path.exists(cds_faa):
            os.link(cds_faa, target_dir + "/Sequence_database/cds.faa")
        if os.path.exists(reshape_gtf):
            os.link(reshape_gtf, target_dir + "/Sequence_database/reshape.gtf")

        # Annotation
        os.makedirs(target_dir + "/Annotation")
        gene_name = {}
        gene_description = {}
        go_lev2 = {}
        kegg_info = {}
        with open(origin_dir + "/annot_class/anno_stat/all_anno_detail.xls", "r") as f:
            f.readline()
            for line in f:
                item = line.strip("\n").split("\t")
                if item[0] not in gene_name:
                    gene_name[item[0]] = item[1]
                    gene_description[item[0]] = item[2]
                    go_lev2[item[0]] = item[18]
                    if "".join([item[11], item[10], item[12], item[12]]) != "":
                        kegg_info[item[0]] = "\t".join([item[11], item[10], item[12], item[12]])
        ## nr
        os.makedirs(target_dir + "/Annotation/nr")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/blast/nr.xls"):
            os.link(origin_dir + "/annot_class/anno_stat/blast/nr.xls", target_dir + "/Annotation/nr/nr.xls")
        ## go
        os.makedirs(target_dir + "/Annotation/go")
        go_stat = target_dir + "/Annotation/go/go_stat.xls"
        go_detail = target_dir + "/Annotation/go/go_detail.xls"
        go_level_statistics = target_dir + "/Annotation/go/go_level_statistics.xls"
        if os.path.exists(origin_dir + "/annot_class/go/GO.stat.xls"):
            with open(origin_dir + "/annot_class/go/GO.stat.xls", "r") as f, open(go_stat, "w") as w:
                w.write("\t".join(["GO Category No.", "Gene No. of CC", "Gene No. of MF", "Gene No. of BP", "Gene No.",
                                   "Percent of All Genes (%)"]) + "\n")
                lines = f.readlines()
                go_num = lines[0].strip().split("\t")[1]
                go_cc = lines[1].strip().split("\t")[1]
                go_mf = lines[2].strip().split("\t")[1]
                go_bp = lines[3].strip().split("\t")[1]
                go_total = lines[4].strip().split("\t")[1]
                go_percent = float(lines[4].strip().split("\t")[1]) / (len(gene_name)) * 100
                w.write(
                    "\t".join([str(go_num), str(go_cc), str(go_mf), str(go_bp), str(go_total), str(go_percent)]) + "\n")
        with open(go_detail, "w") as w:
            w.write("\t".join(["Gene ID", "Gene Name", "Description", "GO ID（Level2）"]) + "\n")
            for gene_id in sorted(gene_name):
                if gene_id in go_lev2:
                    if go_lev2[gene_id] != "":
                        w.write("\t".join(
                            [gene_id, gene_name[gene_id], gene_description[gene_id], go_lev2[gene_id]]) + "\n")
                # else:
                #     w.write("\t".join([gene_id, gene_name[gene_id], gene_description[gene_id], ""]) + "\n")
        if os.path.exists(origin_dir + "/annot_class/go/go12level_statistics.xls"):
            os.link(origin_dir + "/annot_class/go/go12level_statistics.xls", go_level_statistics)
        ## cog
        os.makedirs(target_dir + "/Annotation/cog")
        cog_stat = target_dir + "/Annotation/cog/cog_stat.xls"
        cog_summary = target_dir + "/Annotation/cog/cog_summary.xls"
        cog_detail = target_dir + "/Annotation/cog/cog_detail.xls"
        if os.path.exists(origin_dir + "/annot_class/cog/cog_stat.xls"):
            os.link(origin_dir + "/annot_class/cog/cog_stat.xls", cog_stat)
        if os.path.exists(origin_dir + "/annot_class/cog/cog_summary.xls"):
            os.link(origin_dir + "/annot_class/cog/cog_summary.xls", cog_summary)
        if os.path.exists(origin_dir + "/annot_class/cog/cog.xls"):
            os.link(origin_dir + "/annot_class/cog/cog.xls", cog_detail)
        ## swissprot
        os.makedirs(target_dir + "/Annotation/swissprot")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/blast/swissprot.xls"):
            os.link(origin_dir + "/annot_class/anno_stat/blast/swissprot.xls",
                    target_dir + "/Annotation/swissprot/swissprot.xls")
        ## pfam
        os.makedirs(target_dir + "/Annotation/pfam")
        if os.path.exists(origin_dir + "/annot_class/pfam_domain"):
            os.link(origin_dir + "/annot_class/pfam_domain", target_dir + "/Annotation/pfam/pfam_details.xls")
        pfam_stat = target_dir + "/Annotation/pfam/pfam_statistics.xls"
        gene_pfam = {}
        with open(origin_dir + "/annot_class/pfam_domain", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if gene_pfam.has_key(line[0]):
                    gene_pfam[line[0]][0].append(line[2])
                    gene_pfam[line[0]][1].append(line[3])
                else:
                    gene_pfam.update({line[0]: ([line[2]], [line[3]])})
        with open(pfam_stat, "w") as w:
            w.write("\t".join(["Gene ID", "Gene Name", "Description", "Pfam No.", "Pfam ID", "Domain"]) + "\n")
            for gene_id in gene_pfam:
                w.write("\t".join(
                    [gene_id, gene_name[gene_id], gene_description[gene_id], str(len(gene_pfam[gene_id][0])),
                     ";".join(gene_pfam[gene_id][0]), ";".join(gene_pfam[gene_id][1])]) + "\n")
        ## kegg
        os.makedirs(target_dir + "/Annotation/kegg")
        kegg_layer = origin_dir + "/annot_class/kegg/kegg_layer.xls"

        # added by zhangyitong on 20211011
        pathway_file = os.path.join(origin_dir, 'annot_class/kegg/pathways')
        pngs = os.listdir(pathway_file)
        tar_file = pathway_file + ".tar.gz"
        with tarfile.open(tar_file, mode='w:gz') as f:
            for png in pngs:
                f.add(pathway_file + "/" + png, arcname=png)
        # shutil.rmtree(pathway_file)

        pathway = origin_dir + "/annot_class/kegg/pathways.tar.gz"
        if os.path.exists(kegg_layer):
            with open(target_dir + "/Annotation/kegg/" + "kegg_layer.xls", "w") as w, open(kegg_layer, "r") as f:
                w.write("\t".join(["First Category", "Second Category", "Gene No.", "Gene list"]) + "\n")
                for line in f:
                    w.write(line)
        with open(target_dir + "/Annotation/kegg/" + "kegg_table.xls", "w") as w:
            w.write("\t".join(
                ["Gene ID", "Gene Name", "Description", "KEGG Name", "KO ID", "KO Description", "Pathway ID"]) + "\n")
            for gene_id in gene_name:
                if gene_id in kegg_info:
                    w.write(
                        "\t".join([gene_id, gene_name[gene_id], gene_description[gene_id], kegg_info[gene_id]]) + "\n")
        if os.path.exists(pathway):
            os.link(pathway, target_dir + "/Annotation/kegg/" + "pathways.tar.gz")
        ## summary
        os.makedirs(target_dir + "/Annotation/summary")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/all_anno_detail.xls"):
            os.link(origin_dir + "/annot_class/anno_stat/all_anno_detail.xls",
                    target_dir + "/Annotation/summary/all_anno_detail.xls")
        if os.path.exists(origin_dir + "/annot_class/anno_stat/all_annotation_statistics.xls"):
            with open(origin_dir + "/annot_class/anno_stat/all_annotation_statistics.xls", "r") as f, open(
                    target_dir + "/Annotation/summary/all_annotation_statistics.xls", "w") as w:
                head = f.readline()
                w.write("\t".join(["Type", "genes", "genes_percent"]) + "\n")
                for line in f:
                    w.write(line)

        # Annotation Chart
        cog_pdf = os.path.join(self.chart.work_dir, 'annot_cog.cog_bar.pdf')
        new_path = os.path.join(target_dir, 'Annotation/cog', 'cog_level.pdf')
        self.move_pdf(cog_pdf, new_path)

        go_bar_pdf = os.path.join(self.chart.work_dir, 'annot_go_bar.multi_bar.pdf')
        new_path = os.path.join(target_dir, 'Annotation/go', 'go_level_bar.pdf')
        self.move_pdf(go_bar_pdf, new_path)

        go_pie_pdf = os.path.join(self.chart.work_dir, 'annot_go_pie.multi_pie.pdf')
        new_path = os.path.join(target_dir, 'Annotation/go', 'go_level_pie.pdf')
        self.move_pdf(go_pie_pdf, new_path)

        kegg_pdf = os.path.join(self.chart.work_dir, 'annot_kegg_bar.kegg_bar.pdf')
        new_path = os.path.join(target_dir, 'Annotation/kegg', 'kegg_level.pdf')
        self.move_pdf(kegg_pdf, new_path)

        for i in [['stats_annot_venn.venn.pdf', 'all_anno_venn.pdf'],
                  ['stats_annot_num_bar.bar.pdf', 'all_anno_number.pdf'],
                  ['stats_annot_percent_bar.bar.pdf', 'all_anno_percent.pdf']]:
            pdf = os.path.join(self.chart.work_dir, i[0])
            new_path = os.path.join(target_dir, 'Annotation/summary', i[1])
            self.move_pdf(pdf, new_path)


        # Express
        os.mkdir(target_dir + "/Express")
        ## ExpAnnalysis
        os.mkdir(target_dir + "/Express/ExpAnnalysis")
        if os.path.exists(self.express.output_dir + "/transcript.count.matrix"):
            trans_count = os.path.join(self.express.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
        if os.path.exists(self.express.output_dir + "/transcript.tpm.matrix"):
            trans_tpm = os.path.join(self.express.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
        if os.path.exists(self.express.output_dir + "/transcript.tpm.matrix.annot.xls"):
            trans_tpm_anno = os.path.join(self.express.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")
        if os.path.exists(self.express.output_dir + "/transcript.fpkm.matrix.annot.xls"):
            trans_fpkm_anno = os.path.join(self.express.output_dir + "/transcript.fpkm.matrix.annot.xls")
            os.link(trans_fpkm_anno, target_dir + "/Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls")
        if os.path.exists(self.express_total.output_dir + "/transcript.count.matrix"):
            trans_count = os.path.join(self.express_total.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript_total.count.matrix.xls")
        if os.path.exists(self.express_total.output_dir + "/transcript.tpm.matrix"):
            trans_tpm = os.path.join(self.express_total.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript_total.tpm.matrix.xls")
        if os.path.exists(self.express_total.output_dir + "/transcript.tpm.matrix.annot.xls"):
            trans_tpm_anno = os.path.join(self.express_total.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript_total.tpm.matrix.annot.xls")
        if os.path.exists(self.express_total.output_dir + "/transcript.fpkm.matrix.annot.xls"):
            trans_fpkm_anno = os.path.join(self.express_total.output_dir + "/transcript.fpkm.matrix.annot.xls")
            os.link(trans_fpkm_anno, target_dir + "/Express/ExpAnnalysis/transcript_total.fpkm.matrix.annot.xls")
        ## ExpDistribution
        os.mkdir(target_dir + "/Express/ExpDistribution")
        ## ExpVenn
        os.mkdir(target_dir + "/Express/ExpVenn")
        ## ExpCorr
        os.mkdir(target_dir + "/Express/ExpCorr")
        gene_corr = os.path.join(self.exp_corr.output_dir + "/sample_correlation.xls")
        os.link(gene_corr, target_dir + "/Express/ExpCorr/sample_correlation.xls ")
        ## ExpPCA
        if self.option("group_table").prop["sample_number"] > 2:
            os.mkdir(target_dir + "/Express/ExpPCA")
            gene_pca = os.path.join(self.exp_pca.output_dir + "/PCA.xls")
            os.link(gene_pca, target_dir + "/Express/ExpPCA/PCA.xls")
            gene_ratio = os.path.join(self.exp_pca.output_dir + "/Explained_variance_ratio.xls")
            with open(gene_ratio, "r") as f, open(target_dir + "/Express/ExpPCA/Explained_variance_ratio.xls",
                                                  "w") as w:
                w.write("\tProportion of Variance\n")
                for line in f:
                    w.write(line)

        # Express Chart
        distribution_pdfs = glob.glob(os.path.join(self.chart.work_dir, "*exp_distribution*.pdf"))
        for each in distribution_pdfs:
            dis_type = os.path.basename(each).split('.')[0]
            if each.endswith('.density.pdf'):
                new_file = dis_type + '_exp_density.pdf'
            elif each.endswith('.box.pdf'):
                new_file = dis_type + '_exp_box.pdf'
            elif each.endswith('.violin.pdf'):
                new_file = dis_type + '_exp_violin.pdf'
            new_path = os.path.join(target_dir, 'Express/ExpDistribution', new_file)
            self.move_pdf(each, new_path)

        venn_pdf = os.path.join(self.chart.work_dir, 'all.exp.venn.pdf')
        new_path = os.path.join(target_dir, 'Express/ExpVenn', 'sample_venn.pdf')
        self.move_pdf(venn_pdf, new_path)

        corr_pdf = os.path.join(self.chart.work_dir, 'exp.heatmap.heat_corr.pdf')
        new_path = os.path.join(target_dir, 'Express/ExpCorr', 'sample_correlation.pdf')
        self.move_pdf(corr_pdf, new_path)

        pca_pdf = os.path.join(self.chart.work_dir, 'all.exp_relation_pca.scatter.pdf')
        new_path = os.path.join(target_dir, 'Express/ExpPCA', 'sample_pca.pdf')
        self.move_pdf(pca_pdf, new_path)
        if os.path.exists(os.path.join(self.chart.work_dir, 'all.exp_relation_pca_ell.scatter.pdf')):
            new_path = os.path.join(target_dir, 'Express/ExpPCA', 'sample_pca_ell.pdf')
            self.move_pdf(os.path.join(self.chart.work_dir, 'all.exp_relation_pca_ell.scatter.pdf'), new_path)


        # DiffExpress
        os.mkdir(target_dir + "/DiffExpress")
        if self.option("diff_method").lower() == "deseq2":
            deseq2_xls = glob.glob(self.diffexpress.output_dir + "/*.deseq2.xls")
            deseq2_annot_xls = glob.glob(self.diffexpress.output_dir + "/*.deseq2.annot.xls")
            # deseq2_sizeFactor_xls = glob.glob(self.diffexpress.output_dir + "/*.deseq2.sizeFactor.xls")
            deseq2_sizeFactor_xls = glob.glob(self.diffexpress.tool.work_dir + "/*.deseq2.sizeFactor.xls")
            # deseq2_normalize_xls = glob.glob(self.diffexpress.output_dir + "/*.deseq2.normalize.xls")
            deseq2_normalize_xls = glob.glob(self.diffexpress.tool.work_dir + "/*.deseq2.normalize.xls")
            for file1 in deseq2_xls:
                os.link(file1, target_dir + "/DiffExpress/" + os.path.basename(file1))
            for file2 in deseq2_annot_xls:
                os.link(file2, target_dir + "/DiffExpress/" + os.path.basename(file2))
            for file3 in deseq2_sizeFactor_xls:
                os.link(file3, target_dir + "/DiffExpress/" + os.path.basename(file3))
            for file4 in deseq2_normalize_xls:
                os.link(file4, target_dir + "/DiffExpress/" + os.path.basename(file4))
            # DESeq2_diff_summary = self.diffexpress.output_dir + "/DESeq2_diff_summary.xls"
            DESeq2_diff_summary = self.diffexpress.output_dir + "/diff_summary_deseq2.xls"
            os.link(DESeq2_diff_summary, target_dir + "/DiffExpress/diff_summary_deseq2.xls")
        else:
            results = glob.glob(self.diffexpress.output_dir + "/*.xls")
            for file in results:
                os.link(file, target_dir + "/DiffExpress/" + os.path.basename(file))

        # DiffExpress Chart
        volcano_pdfs = glob.glob(os.path.join(self.chart.work_dir, "*.diffexp.*.pdf"))
        for each in volcano_pdfs:
            cmp, _, suffix = os.path.basename(each).split('.', 2)
            new_path = os.path.join(target_dir, 'DiffExpress', cmp + '.' + suffix)
            self.move_pdf(each, new_path)

        bar_pdf = glob.glob(os.path.join(self.chart.work_dir, '*.diffexp_summary.bar.pdf'))
        if bar_pdf:
            diff_method = os.path.basename(bar_pdf[0]).split('.')[0]
            new_path = os.path.join(target_dir, 'DiffExpress', diff_method + '_diff_bar.pdf')
            self.move_pdf(bar_pdf[0], new_path)

        stacked_pdf = glob.glob(os.path.join(self.chart.work_dir, '*.diffexp_summary.stacked_bar.pdf'))
        if stacked_pdf:
            new_path = os.path.join(target_dir, 'DiffExpress', diff_method + '_diff_summary.pdf')
            self.move_pdf(stacked_pdf[0], new_path)

        # srna
        if self.option("strand_specific"):
            os.makedirs(target_dir + "/sRNA")
            os.makedirs(target_dir + "/sRNA/srna_predict")
            os.makedirs(target_dir + "/sRNA/srna_annot")
            os.makedirs(target_dir + "/sRNA/srna_fold")
            os.makedirs(target_dir + "/sRNA/srna_target")
            # srna/srna_predicted
            genome_predicted_RNA_fa = origin_dir + "/rockhopper/genome.predicted_RNA.fa"
            if os.path.exists(genome_predicted_RNA_fa):
                os.link(genome_predicted_RNA_fa, target_dir + "/sRNA/srna_predict/genome.predicted_RNA.fa")
            srna_detail = origin_dir + "/rockhopper/genome.predicted_RNA.bed.xls"
            if os.path.exists(srna_detail):
                os.link(srna_detail, target_dir + "/sRNA/srna_predict/sran_detail.xls")
            # srna/srna_annot
            annotation_stat_xls = origin_dir + "/srna/srna_annot/annotation_stat.xls"
            if os.path.exists(annotation_stat_xls):
                os.link(annotation_stat_xls, target_dir + "/sRNA/srna_annot/annotation_stat.xls")
            rfam_stat_xls = origin_dir + "/srna/srna_annot/rfam_stat.xls"
            if os.path.exists(rfam_stat_xls):
                with open(rfam_stat_xls, "r") as f, open(target_dir + "/sRNA/srna_annot/rfam_stat.xls", "w") as w:
                    w.write("\t".join(["Type", "Total Number", "Total Percent (%)"]) + "\n")
                    for line in f:
                        w.write(line)
            rfam_detail_xls = origin_dir + "/srna/srna_annot/genome.predicted_RNA.fa_vs_rfam.xls"
            if os.path.exists(rfam_detail_xls):
                os.link(rfam_detail_xls, target_dir + "/sRNA/srna_annot/rfam_detail.xls")
            sRNAMap_detail_xls = origin_dir + "/srna/srna_annot/genome.predicted_RNA.fa_vs_sRNAMap.xls"
            if os.path.exists(sRNAMap_detail_xls):
                os.link(sRNAMap_detail_xls, target_dir + "/sRNA/srna_annot/sRNAMap_detail.xls")
            sRNATarBase_detail_xls = origin_dir + "/srna/srna_annot/genome.predicted_RNA.fa_vs_sRNATarBase.xls"
            if os.path.exists(sRNATarBase_detail_xls):
                os.link(sRNATarBase_detail_xls, target_dir + "/sRNA/srna_annot/sRNATarBase_detail.xls")
            SIPHT_detail_xls = origin_dir + "/srna/srna_annot/genome.predicted_RNA.fa_vs_SIPHI.xls"
            if os.path.exists(SIPHT_detail_xls):
                os.link(SIPHT_detail_xls, target_dir + "/sRNA/srna_annot/SIPHT_detail.xls")
            bsrd_detail_xls = origin_dir + "/srna/srna_annot/genome.predicted_RNA.fa_vs_BSRD.xls"
            if os.path.exists(bsrd_detail_xls):
                os.link(bsrd_detail_xls, target_dir + "/sRNA/srna_annot/BSRD_detail.xls")

            # srna/srna_fold
            pdf = glob.glob(origin_dir + "/srna/srna_fold/*pdf")
            for file in pdf:
                os.link(file, target_dir + "/sRNA/srna_fold/" + os.path.basename(file))
            sRNA_stat = origin_dir + "/srna/srna_fold/RNAfold_stat.txt"
            if os.path.exists(sRNA_stat):
                os.link(sRNA_stat, target_dir + "/sRNA/srna_fold/sRNA_stat.xls")
            sRNA_str = origin_dir + "/srna/srna_fold/RNAfold.str"
            if os.path.exists(sRNA_str):
                os.link(sRNA_str, target_dir + "/sRNA/srna_fold/RNAfold.str")

            # srna/srna_target
            for file in glob.glob(origin_dir + "/srna/srna_target/*"):
                os.link(file, target_dir + "/sRNA/srna_target/" + os.path.basename(file) + ".xls")

            # srna chart
            for i in [['srna_length_bar.bar.pdf', 'srna_predict/srna_length.pdf'],
                      ['srna_annot_venn.venn.pdf', 'srna_annot/annotation_venn.pdf'],
                      ['srna_rfam_pie.pie.pdf', 'srna_annot/rfam_stat.pdf']]:
                new_file = os.path.join(target_dir, 'sRNA', i[1])
                self.move_pdf(os.path.join(self.chart.work_dir, i[0]), new_file)

        # GeneStructure
        os.makedirs(target_dir + "/GeneStructure")
        if self.option("strand_specific"):
            os.makedirs(target_dir + "/GeneStructure/operon/")
            os.makedirs(target_dir + "/GeneStructure/TSS_and_TTS")
            os.makedirs(target_dir + "/GeneStructure/UTR/")
            os.makedirs(target_dir + "/GeneStructure/promoter/")
            os.makedirs(target_dir + "/GeneStructure/terminator/")

            operon_xls = origin_dir + "/rockhopper/operon.xls"
            os.link(operon_xls, target_dir + "/GeneStructure/operon/operon.xls")
            tss_and_tts = origin_dir + "/rockhopper/TSS_and_TTS.xls"
            os.link(tss_and_tts, target_dir + "/GeneStructure/TSS_and_TTS/TSS_and_TTS.xls")
            UTR_xls = origin_dir + "/rockhopper/UTR.xls"
            os.link(UTR_xls, target_dir + "/GeneStructure/UTR/UTR.xls")
            UTR3_fa = origin_dir + "/rockhopper/UTR3.fa"
            os.link(UTR3_fa, target_dir + "/GeneStructure/UTR/UTR3.fa")
            UTR5_fa = origin_dir + "/rockhopper/UTR5.fa"
            os.link(UTR5_fa, target_dir + "/GeneStructure/UTR/UTR5.fa")
            terminator = origin_dir + "/terminator/terminator_prediction.xls"
            os.link(terminator, target_dir + "/GeneStructure/terminator/terminator_detail.xls")
            promote = origin_dir + "/promote/all_promoter_result.xls"
            promote_result = target_dir + "/GeneStructure/promoter/promoter_detail.xls"
            with open(promote, "r") as f, open(promote_result, "w") as w:
                head = f.readline()
                items = head.split("\t")
                w.write("Gene ID\tGene name\tGene description\t" + "\t".join(items[1:]))
                for line in f:
                    items = line.split("\t")
                    gene_id = items[0]
                    if gene_id in gene_name:
                        name = gene_name[gene_id]
                    else:
                        name = "-"
                    if gene_id in gene_description:
                        description = gene_description[gene_id]
                    else:
                        description = "-"
                    w.write("\t".join([gene_id, name, description]) + "\t".join(items[1:]))

        # SNP
        os.makedirs(target_dir + "/GeneStructure/SNP/")
        CopyFile().linkdir(self.work_dir + "/snp_tmp", target_dir + "/GeneStructure/SNP/")
        self.merge1(target_dir + "/GeneStructure/SNP/data_anno_pre.xls",
                    target_dir + "/GeneStructure/SNP/snp_annotation.xls")
        if os.path.exists(target_dir + "/GeneStructure/SNP/snp_annotation.xls"):
            os.remove(target_dir + "/GeneStructure/SNP/snp_annotation.xls")
        if os.path.exists(target_dir + "/GeneStructure/SNP/data_anno_pre.xls"):
            os.remove(target_dir + "/GeneStructure/SNP/data_anno_pre.xls")

        # GeneStructure Chart
        for i in [['operon_gene_num_bar.bar.pdf', 'operon/operon_gene.num.pdf'],
                  ['operon_len_bar.bar.pdf', 'operon/operon_length.pdf'],
                  ['UTR3_len_bar.bar.pdf', 'UTR/UTR3_length.pdf'],
                  ['UTR5_len_bar.bar.pdf', 'UTR/UTR5_length.pdf']]:
            new_file = os.path.join(target_dir, 'GeneStructure', i[1])
            self.move_pdf(os.path.join(self.chart.work_dir, i[0]), new_file)

        for i in ['pie', 'bar']:
            pie_pdfs = glob.glob(os.path.join(self.chart.work_dir, '*_stats_{}.{}.pdf'.format(i, i)))
            for each in pie_pdfs:
                sample, suffix = os.path.basename(each).split('.', 1)
                if 'type_' in suffix:
                    new_file = sample + '_snp_type_{}.pdf'.format(i)
                elif 'frequency_' in suffix:
                    new_file = sample + '_snp_frequency_{}.pdf'.format(i)
                elif 'depth_' in suffix:
                    new_file = sample + '_snp_depth_{}.pdf'.format(i)
                new_path = os.path.join(target_dir, 'GeneStructure/SNP', new_file)
                self.move_pdf(each, new_path)

        pdfs = glob.glob(os.path.join(self.chart.work_dir, '*_regions_pie.pie.pdf'))
        for each in pdfs:
            sample, suffix = os.path.basename(each).split('.', 1)
            if 'indel_' in suffix:
                new_file = 'indel_distribution_pie_{}.pdf'.format(sample)
            if 'snp_' in suffix:
                new_file = 'snp_distribution_pie_{}.pdf'.format(sample)
            new_path = os.path.join(target_dir, 'GeneStructure/SNP', new_file)
            self.move_pdf(each, new_path)

        # Geneset Pipeline
        if not os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir,"results_info")):
            os.makedirs(target_dir + "/Geneset")
            os.makedirs(target_dir + "/Geneset/GenesetVenn")
            os.makedirs(target_dir + "/Geneset/Geneset_Cluster")
            os.makedirs(target_dir + "/Geneset/GenesetCogClass")
            os.makedirs(target_dir + "/Geneset/GenesetGoClass")
            os.makedirs(target_dir + "/Geneset/GenesetKeggClass")
            os.makedirs(target_dir + "/Geneset/GenesetGoEnrich")
            os.makedirs(target_dir + "/Geneset/GenesetKeggEnrich")
            os.makedirs(target_dir + "/Geneset/GenesetIpath")

            if os.path.exists(os.path.join(self.chart.work_dir, 'geneset.venn.pdf')):
                new_path = os.path.join(target_dir, 'Geneset/GenesetVenn/venn.pdf')
                self.move_pdf(os.path.join(self.chart.work_dir, 'geneset.venn.pdf'), new_path)

            cluster_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'cluster/All_Diff_mRNA/*'))
            for each in cluster_files:
                if each.endswith('tre'):
                    continue
                else:
                    self.move_pdf(each, os.path.join(target_dir, 'Geneset/Geneset_Cluster', os.path.basename(each)))
            cluster_pdfs = glob.glob(os.path.join(self.chart.work_dir, 'subcluster*pdf'))
            for each in cluster_pdfs:
                num = os.path.basename(each).split('.')[1]
                new_path = os.path.join(target_dir, 'Geneset/Geneset_Cluster/seq.subcluster_{}.pdf'.format(num))
                self.move_pdf(each, new_path)
            pdf = os.path.join(self.chart.work_dir, 'geneset.cluster.heat_corr.pdf')
            self.move_pdf(pdf, os.path.join(target_dir, 'Geneset/Geneset_Cluster/expression_matrix.pdf'))

            cog_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_cog_class/*'))
            for each in cog_files:
                self.move_pdf(each, os.path.join(target_dir, 'Geneset/GenesetCogClass', os.path.basename(each)))
            pdf = os.path.join(self.chart.work_dir, 'cog_annot.gene_set.column.pdf')
            self.move_pdf(pdf, os.path.join(target_dir, 'Geneset/GenesetCogClass/cog_class.pdf'))

            goclass_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_go_class/*'))
            for each in goclass_files:
                self.move_pdf(each, os.path.join(target_dir, 'Geneset/GenesetGoClass', os.path.basename(each)))
            pdf = glob.glob(os.path.join(self.chart.work_dir, '*.go_bar_geneset.go_bar.pdf'))
            if pdf:
                self.move_pdf(pdf[0], os.path.join(target_dir, 'Geneset/GenesetGoClass/go_class.pdf'))

            pathway_file = os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_kegg_class/pathways')
            pngs = os.listdir(pathway_file)
            tar_file = pathway_file + ".tar.gz"
            with tarfile.open(tar_file, mode='w:gz') as f:
                for png in pngs:
                    f.add(pathway_file + "/" + png, arcname=png)
            # shutil.rmtree(pathway_file)
            self.move_pdf(tar_file, os.path.join(target_dir, 'Geneset/GenesetKeggClass/pathways.tar.gz'))
            keggclass_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_kegg_class/kegg_stat.xls'))
            self.move_pdf(keggclass_files[0], os.path.join(target_dir, 'Geneset/GenesetKeggClass', 'kegg_stat.xls'))
            pdf = glob.glob(os.path.join(self.chart.work_dir, '*annot_kegg_bar.kegg_bar.pdf'))[0]
            self.move_pdf(pdf, os.path.join(target_dir, 'Geneset/GenesetKeggClass/kegg_class.pdf'))

            goenrich_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_go_enrich/*'))
            for each in goenrich_files:
                self.move_pdf(each, os.path.join(target_dir, 'Geneset/GenesetGoEnrich', os.path.basename(each)))
            for i in [['*_all_go_enrich_dense_bubble.densebubble.pdf', 'go_enrich_bubble.pdf'],
                      ['*_all_go_enrich_scatter_bubble.scatterbubble.pdf', 'go_enrich_disbubble.pdf'],
                      ['*_go_enrich_bar.shadowbar.pdf', 'go_enrich_bar.pdf']]:
                pdf = glob.glob(os.path.join(self.chart.work_dir, i[0]))
                if pdf:
                    self.move_pdf(pdf[0], os.path.join(target_dir, 'Geneset/GenesetGoEnrich', i[1]))

            pathway_file = os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_kegg_enrich/class/pathways')
            pngs = os.listdir(pathway_file)
            tar_file = pathway_file + ".tar.gz"
            with tarfile.open(tar_file, mode='w:gz') as f:
                for png in pngs:
                    f.add(pathway_file + "/" + png, arcname=png)
            self.move_pdf(tar_file, os.path.join(target_dir, 'Geneset/GenesetKeggEnrich/pathways.tar.gz'))
            keggenrich_files = glob.glob(
                os.path.join(self.diff_geneset_analysis.output_dir, 'All_Diff_mRNA/diff_kegg_enrich/enrich/*.xls'))
            for each in keggenrich_files:
                self.move_pdf(each, os.path.join(target_dir, 'Geneset/GenesetKeggEnrich', 'kegg_enrich_{}.xls'.format(os.path.basename(each).split('.')[0])))
            for i in [['*_enrichkegg.densebubble.pdf', 'pathway_enrich_bubble.pdf'],
                      ['*_enrichkegg.shadowbar.pdf', 'pathway_enrich_bar.pdf']]:
                pdf = glob.glob(os.path.join(self.chart.work_dir, i[0]))
                if pdf:
                    self.move_pdf(pdf[0], os.path.join(target_dir, 'Geneset/GenesetKeggEnrich', i[1]))

            ipath_files = glob.glob(os.path.join(self.diff_geneset_analysis.output_dir, 'ipath/All_Diff_mRNA/*'))
            for each in ipath_files:
                self.move_pdf(each, os.path.join(target_dir, 'Geneset/GenesetIpath', os.path.basename(each)))

        # upload png for report_model
        if self.option("report_img"):
            s3 = self._sheet.output.split(":")[0]
            report_img_dir = self.chart.work_dir + '/png/'
            report_img_s3 = s3 + "://commonbucket/files/report_img/prokrna/" + self.task_id + "/"
            self.upload_to_s3(report_img_dir, report_img_s3)

        #差异一键化
        # shutil.copytree(self.diff_geneset_analysis.output_dir, os.path.join(target_dir, "diff_geneset_analysis"))
        sdir = self.add_upload_dir(target_dir)
        regexp_path = [
            [r"QC/.*_raw_bases_distribution\.pdf", "pdf", "Raw Reads碱基组成分布图"],
            [r"QC/.*_raw_bases_error_rate\.pdf", "pdf", "Raw Reads碱基错误率分布图"],
            [r"QC/.*_raw_bases_quality\.pdf", "pdf", "Raw Reads碱基质量分布图"],
            [r"QC/.*_clean_bases_distribution\.pdf", "pdf", "Clean Reads碱基组成分布图"],
            [r"QC/.*_clean_bases_error_rate\.pdf", "pdf", "Clean Reads碱基错误率分布图"],
            [r"QC/.*_clean_bases_quality\.pdf", "pdf", "Clean Reads碱基质量分布图"],
            [r"Align/.*_saturation_curve\.pdf", 'pdf', '测序饱和度曲线图'],
            [r"Express/ExpDistribution/.*_exp_density\.pdf", "pdf", "表达量分布密度图"],
            [r"Express/ExpDistribution/.*_exp_box\.pdf", "pdf", "表达量分布盒形图"],
            [r"Express/ExpDistribution/.*_exp_violin\.pdf", "pdf", "表达量分布小提琴图"],
            [r"DiffExpress/.*_vs_.*\.xls", "xls", "差异分析结果表"],
            [r"DiffExpress/.*summary\.xls", "xls", "差异统计结果表"],
            [r"DiffExpress/.*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表"],
            [r"DiffExpress/.*_vs_.*\..*annot\.xls", "xls", "差异统计注释结果表"],
            [r"DiffExpress/.*_vs_.*\..*normalize\.xls", "xls", "差异矩阵标准化结果表"],
            [r"DiffExpress/.*_vs_.*\..*sizeFactor\.xls", "xls", "差异矩阵标准化因子表"],
            [r"DiffExpress/.*_vs_.*\..*normFactor\.xls", "xls", "差异矩阵标准化因子表"],
            [r"DiffExpress/.*_vs_.*\.scatter\.pdf", "pdf", "表达量差异散点图"],
            [r"DiffExpress/.*_vs_.*\.volcano\.pdf", "pdf", "表达量差异火山图"],
            [r"DiffExpress/.*_diff_bar\.pdf", "pdf", "差异统计统计柱状图"],
            [r"DiffExpress/.*_diff_summary\.pdf", "pdf", "差异统计统计堆积图"],
            [r"Align/AlignBam/.*\.bam", "", "样本比对bam文件"],
            [r"GeneStructure/SNP/snp_distribution_pie_.*\.pdf", "pdf", "SNP不同区域分布饼图"],
            [r"GeneStructure/SNP/indel_distribution_pie_.*\.pdf", "pdf", "InDel不同区域分布饼图"],
            [r"GeneStructure/SNP/.*_snp_type_bar\.pdf", "pdf", "SNP类型统计柱状图"],
            [r"GeneStructure/SNP/.*_snp_type_pie\.pdf", "pdf", "SNP类型统计饼图"],
            [r"GeneStructure/SNP/.*_snp_frequency_bar\.pdf", "pdf", "SNP频率统计柱状图"],
            [r"GeneStructure/SNP/.*_snp_frequency_pie\.pdf", "pdf", "SNP频率统计饼图"],
            [r"GeneStructure/SNP/.*_snp_depth_bar\.pdf", "pdf", "SNP深度统计柱状图"],
            [r"GeneStructure/SNP/.*_snp_depth_pie\.pdf", "pdf", "SNP深度统计饼图"],
            [r"Geneset/Geneset_Cluster/seq\.subcluster_.*\.xls", "", "子聚类分析表"],
            [r"Geneset/Geneset_Cluster/seq\.subcluster_.*\.pdf", "", "子聚类趋势图"],
            [r"Geneset/GenesetGoEnrich/adjust_lineage\.(svg|png|pdf)", "", "显著富集GO的层级结构图"],
            [r"Geneset/GenesetGoEnrich/go_lineage\.(svg|png|pdf)", "", "GO 的层级结构图"],
            [r"Geneset/GenesetGoEnrich/go_enrich_.*_gene\.xls", "", "GO富集分析统计表"],
            [r"Geneset/GenesetKeggEnrich/kegg_enrich_.*_gene\.xls", "", "KEGG富集分析统计表"],
        ]
        rel_path = [
            [".", "", "流程分析结果目录"],
            ["Background", "", "项目背景目录"],
            ["Background/genome_stat.xls", "", "参考基因组注释信息表"],
            ["Background/sample_info.xls", "", "样本信息表"],
            ["Background/software_info.xls", "", "软件信息表"],
            ["Sequence_database", "", "序列文件数据库目录", 1],
            ["Sequence_database/cds.fa", "", "参考基因组核酸序列文件"],
            ["Sequence_database/reshape.gtf", "", "参考基因组gtf文件"],
            ["Sequence_database/cds.faa", "", "参考基因组蛋白序列文件"],
            ["QC", "", "测序数据统计与质控结果目录"],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表"],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表"],
            ["QC/rrna_assessment.xls", "", "核糖体RNA污染率评估统计表"],
            ["Align", "", "比对结果目录"],
            ["Align/align_stat.xls", "", "比对结果统计表"],
            ["Align/QualityAssessment", "", "比对结果整体评估结果文件"],
            ["Align/QualityAssessment/chr_distribution.xls", "", "Reads在不同染色体的分布统计表"],
            ["Align/coverage_distribution.pdf", 'pdf', '测序覆盖度分布图'],
            ["Annotation", "", "参考基因组注释结果目录", 1],
            ["Annotation/nr", "", "nr注释结果文件"],
            ["Annotation/nr/nr.xls", "", "nr注释详情表"],
            ["Annotation/swissprot", "", "swissprot注释结果文件"],
            ["Annotation/swissprot/swissprot.xls", "", "swissprot注释详情表"],
            ["Annotation/pfam", "", "pfam注释结果文件"],
            ["Annotation/pfam/pfam_details.xls", "", "pfam注释详情表"],
            ["Annotation/pfam/pfam_statistics.xls", "", "pfam注释统计表"],
            ["Annotation/cog", "", "cog注释结果文件"],
            ["Annotation/cog/cog_summary.xls", "", "cog注释详情表"],
            ["Annotation/cog/cog_stat.xls", "", "cog注释统计表"],
            ["Annotation/cog/cog_detail.xls", "", "基因与cog对应关系详情表"],
            ["Annotation/cog/cog_level.pdf", "pdf", "COG分类统计柱状图"],
            ["Annotation/go", "", "go注释结果文件"],
            ["Annotation/go/go_stat.xls", "", "go注释统计表"],
            ["Annotation/go/go_detail.xls", "", "基因与go对应关系详情表"],
            ["Annotation/go/go_level_statistics.xls", "", "go分类统计表"],
            ["Annotation/go/go_level_bar.pdf", "pdf", "GO注释分类统计柱形图"],
            ["Annotation/go/go_level_pie.pdf", "pdf", "GO注释分类统计饼图"],
            ["Annotation/kegg", "", "kegg注释结果目录"],
            ["Annotation/kegg/kegg_table.xls", "", "基因与Pathway对应关系详情表"],
            ["Annotation/kegg/pathways.tar.gz", "", "Pathway通路图压缩文件"],
            ["Annotation/kegg/kegg_layer.xls", "", "Pathway分类统计表"],
            ["Annotation/kegg/kegg_level.pdf", "pdf", "Pathway分类统计柱状图"],
            ["Annotation/summary", "", "注释结果汇总文件"],
            ["Annotation/summary/all_annotation_statistics.xls", "", "注释结果统计表"],
            ["Annotation/summary/all_anno_detail.xls", "", "注释结果汇总详情表"],
            ["Annotation/summary/all_anno_venn.pdf", "pdf", "基础注释统计Venn图"],
            ["Annotation/summary/all_anno_percent.pdf", "pdf", "基础注释统计柱状占比图"],
            ["Annotation/summary/all_anno_number.pdf", "pdf", "基础注释统计柱状图"],
            ["Express", "", "表达量分析结果目录"],
            ["Express/ExpAnnalysis", "", "表达定量结果文件"],
            ["Express/ExpAnnalysis/transcript.count.matrix.xls", "", "transcript count表达定量结果文件"],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls", "", "transcript tpm表达定量结果注释文件"],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量结果注释文件"],
            ["Express/ExpAnnalysis/transcript.tpm.matrix.xls", "", "transcript tpm表达定量结果文件"],
            ["Express/ExpAnnalysis/transcript.fpkm.matrix.xls", "", "transcript fpkm表达定量结果文件"],
            ["Express/ExpAnnalysis/transcript_total.count.matrix.xls", "", "transcript tpm表达定量结果总表"],
            ["Express/ExpAnnalysis/transcript_total.tpm.matrix.xls", "", "transcript tpm表达定量结果总表"],
            ["Express/ExpAnnalysis/transcript_total.fpkm.matrix.xls", "", "transcript fpkm表达定量结果总表"],
            ["Express/ExpAnnalysis/transcript_total.tpm.matrix.annot.xls", "", "transcript tpm表达定量结果总表注释文件"],
            ["Express/ExpAnnalysis/transcript_total.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量结果总表注释文件"],
            ["Express/ExpDistribution", "", "表达量分布结果目录"],
            ["Express/ExpVenn", "", "样本间Venn分析结果目录"],
            ["Express/ExpVenn/sample_venn.pdf", "pdf", "样本间Venn图"],
            ["Express/ExpCorr", "", "样本间相关性分析结果文件"],
            ["Express/ExpCorr/sample_correlation.xls ", "", "样本间相关性分析矩阵表"],
            ["Express/ExpCorr/sample_correlation.pdf", "pdf", "样本间相关性热图"],
            ["Express/ExpPCA", "", "样本间PCA分析结果文件"],
            ["Express/ExpPCA/PCA.xls", "", "样本间PCA分析结果表"],
            ["Express/ExpPCA/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表"],
            ["Express/ExpPCA/sample_pca.pdf", "pdf", "样本间PCA图"],
            ["Express/ExpPCA/sample_pca_ell.pdf", "pdf", "样本间PCA图（置信圈）"],
            ["DiffExpress", "", "表达量差异分析目录"],
            ["GeneStructure", "", "转录本结构分析结果目录"],
            ["GeneStructure/operon", "", "操纵子分析结果文件"],
            ["GeneStructure/operon/operon.xls", "", "操纵子详情表"],
            ["GeneStructure/operon/operon_length.pdf", "pdf", "操纵子长度统计图"],
            ["GeneStructure/operon/operon_gene.num.pdf", "pdf", "操纵子基因数量统计图"],
            ["GeneStructure/UTR", "", "UTR分析结果文件"],
            ["GeneStructure/UTR/UTR.xls", "", "UTR结果文件表"],
            ["GeneStructure/UTR/UTR3.fa", "", "3UTR序列文件"],
            ["GeneStructure/UTR/UTR5.fa", "", "5UTR序列文件"],
            ["GeneStructure/UTR/UTR3_length.pdf", "pdf", "3UTR长度分布图"],
            ["GeneStructure/UTR/UTR5_length.pdf", "pdf", "5UTR长度分布图"],
            ["GeneStructure/promoter", "", "启动子预测结果文件"],
            ["GeneStructure/promoter/promoter_detail.xls", "", "启动子预测结果详情表"],
            ["GeneStructure/terminator", "", "终止子预测结果文件"],
            ["GeneStructure/terminator/terminator_detail.xls", "", "终止子预测结果详情表"],
            ["GeneStructure/TSS_and_TTS", "", "转录起始与终止位点结果文件"],
            ["GeneStructure/TSS_and_TTS/TSS_and_TTS.xls", "", "转录起始与终止位点预测结果表"],
            ["GeneStructure/SNP", "", "snp分析结果文件"],
            ["GeneStructure/SNP/snp_annotation_statistics.xls", "", "snp分析结果注释统计表格"],
            ["GeneStructure/SNP/snp_annotation_detail.xls", "", "snp分析结果注释详情表格"],
            ["GeneStructure/SNP/snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格"],
            ["GeneStructure/SNP/snp_freq_statistics.xls", "", "SNP频率统计结果表格"],
            ["GeneStructure/SNP/snp_depth_statistics.xls", "", "SNP深度统计结果表格"],
            ["GeneStructure/SNP/snp_position_distribution.xls", "", "SNP不同区域分析结果表格"],
            ["GeneStructure/SNP/indel_position_distribution.xls", "", "InDel不同区域分析结果表格"],
            ["Geneset", "", "基因集分析结果目录"],
            ["Geneset/GenesetVenn", "", "基因集venn分析结果"],
            ["Geneset/GenesetVenn/venn.pdf", "", "venn分析结果图"],
            ["Geneset/Geneset_Cluster", "", "基因集聚类分析结果文件"],
            ["Geneset/Geneset_Cluster/expression_matrix.pdf", "", "基因聚类热图"],
            ["Geneset/Geneset_Cluster/expression_matrix.xls", "", "聚类分析结果表"],
            ["Geneset/Geneset_Cluster/seq.cluster_tree.txt", "", "基因聚类树文件"],
            ["Geneset/Geneset_Cluster/sample.cluster_tree.txt", "", "样本间聚类树文件"],
            ["Geneset/GenesetCogClass", "", "基因集COG分类结果文件"],
            ["Geneset/GenesetCogClass/cog_class.pdf", "", "COG分类统计柱状图"],
            ["Geneset/GenesetCogClass/cog_class_table.xls", "", "COG分类统计表"],
            ["Geneset/GenesetCogClass/cog_analysis_detail.xls", "", "基因与COG对应关系详情表"],
            ["Geneset/GenesetGoClass", "", "基因集GO分类结果文件"],
            ["Geneset/GenesetGoClass/go_class.pdf", "", "GO分类统计柱形图"],
            ["Geneset/GenesetGoClass/go_class_table.xls", "", "GO分类统计表"],
            ["Geneset/GenesetGoClass/go_detail_table.xls", "", "基因与GO对应关系详情表"],
            ["Geneset/GenesetKeggClass", "", "基因集KEGG功能分类结果文件"],
            ["Geneset/GenesetKeggClass/kegg_class.pdf", "", "Pathway分类统计柱状图"],
            ["Geneset/GenesetKeggClass/kegg_stat.xls", "", "Pathway分类统计表"],
            ["Geneset/GenesetKeggClass/pathways_analysis_detail.xls", "", "基因与Pathway对应关系详情表"],
            ["Geneset/GenesetKeggClass/pathways.tar.gz", "", "Pathway代谢通路压缩文件"],
            ["Geneset/GenesetGoEnrich", "", "基因集GO富集分析结果文件"],
            ["Geneset/GenesetGoEnrich/go_enrich_bar.pdf", "", "GO富集分析分类统计柱形图"],
            ["Geneset/GenesetGoEnrich/go_enrich_bubble.pdf", "", "GO富集分析分类统计气泡图"],
            ["Geneset/GenesetGoEnrich/go_enrich_disbubble.pdf", "", "GO富集分析分类统计分散型气泡图"],
            ["Geneset/GenesetGoEnrich/go_gene_stat.xls", "", "基因与GO对应关系详情表"],
            ["Geneset/GenesetGoEnrich/adjust_lineage.svg", "", "显著富集GO的层级结构图"],
            ["Geneset/GenesetKeggEnrich", "", "基因集KEGG富集分析结果文件"],
            ["Geneset/GenesetKeggEnrich/pathway_enrich_bar.pdf", "", "KEGG富集分析分类统计柱形图"],
            ["Geneset/GenesetKeggEnrich/pathway_enrich_bubble.pdf", "", "KEGG富集分析分类统计气泡图"],
            ["Geneset/GenesetKeggEnrich/pathway_analysis_detail.xls", "", "基因与Pathway对应关系详情表"],
            ["Geneset/GenesetKeggEnrich/pathways.tar.gz", "", "Pathway通路图压缩文件"],
            ["Geneset/GenesetIpath", "", "基因集Ipath结果文件"],
            ["Geneset/GenesetIpath/Biosynthesis_of_secondary_metabolities.svg", "", "次级代谢产物合成通路图"],
            ["Geneset/GenesetIpath/Regulatory_pathways.svg", "", "调控通路图"],
            ["Geneset/GenesetIpath/Metabolic_pathways.svg", "", "iPath代谢通路图"],
            ["Geneset/GenesetIpath/gene_ipath_input.xls", "", "代谢通路数据表"],

        ]
        if self.option("strand_specific"):
            extra_rel_path = [
                [r"sRNA", "", "sRNA分析结果目录"],
                [r"sRNA/srna_predict", "", "sRNA预测结果文件"],
                [r"sRNA/srna_predict/genome.predicted_RNA.fa", "", "sRNA序列文件"],
                [r"sRNA/srna_predict/srna_detail.xls", "", "sRNA预测详情表"],
                [r"sRNA/srna_predict/srna_length.pdf", "pdf", "sRNA长度统计图"],
                [r"sRNA/srna_annot", "", "sRNA注释分析结果文件"],
                [r"sRNA/srna_annot/annotation_stat.xls", "", "sRNA注释汇总表"],
                [r"sRNA/srna_annot/rfam_stat.xls", "", "Rfam注释统计表"],
                [r"sRNA/srna_annot/rfam_detail.xls", "", "Rfam注释结果详情表"],
                [r"sRNA/srna_annot/sRNAMap_detail.xls", "", "sRNAMap注释结果详情表"],
                [r"sRNA/srna_annot/sRNATarBase_detail.xls", "", "sRNATarBase注释结果详情表"],
                [r"sRNA/srna_annot/SIPHT_detail.xls", "", "SIPHT注释结果详情表"],
                [r"sRNA/srna_annot/BSRD_detail.xls", "", "BSRD注释结果详情表"],
                [r"sRNA/srna_annot/annotation_venn.pdf", "pdf", "sRNA注释Venn图"],
                [r"sRNA/srna_annot/rfam_stat.pdf", "pdf", "Rfam注释统计图"],
                [r"sRNA/srna_fold", "", "sRNA二级结构预测结果文件"],
                [r"sRNA/srna_fold/sRNA_stat.xls", "", "sRNA二级结构预测结果表"],
                [r"sRNA/srna_fold/RNAfold.str", "", "sRNA二级结构文件"],
                [r"sRNA/srna_target", "", "sRNA靶基因预测结果文件"],
                [r"sRNA/srna_target/combine_RNAplex_RNAhybrid.xls", "", "sRNA靶基因预测统计表"],
                [r"sRNA/srna_target/RNAhybrid_merge.xls", "", "RNAhybrid预测详情表"],
                [r"sRNA/srna_target/RNAplex_merge.xls", "", "RNAplex预测详情表"],
                [r"sRNA/srna_fold/*\.pdf", "pdf", "sRNA二级结构图"],
                [r"sRNA/srna_fold/*\.png", "png", "sRNA二级结构图"],
            ]
            regexp_path += extra_rel_path
        sdir.add_relpath_rules(rel_path)
        sdir.add_regexp_rules(regexp_path)

        intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        if not intermediate_dir.endswith("/"):
            intermediate_dir = intermediate_dir + "/"
        if not self.intermediate_result.endswith("/"):
            self.intermediate_result = self.intermediate_result + "/"
        self.upload_to_s3(self.intermediate_result, intermediate_dir)

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
        self.export_productive_name()
        self.export_map_assess()
        if self.option("strand_specific") == True:
            self.export_srna()
            self.export_structure(srna=True)
        self.export_structure(promote=True)
        self.export_structure(terminator=True)
        self.export_structure(antisense=True)
        self.export_structure(novelgene=True)
        self.export_snp()
        # self.export_promote()
        self.export_annotation()
        self.export_expression()
        self.export_rna_type()
        self.export_genesets_analysis()
        self.logger.info("导表完成")

    def run_api_qc_mapping(self):
        self.logger.info('经判断，序列比对或rRNA比例不达标，流程提前终止')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.prok_rna')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.export_qc()
        self.export_productive_name()
        self.export_mapping_only()

    @time_count
    def export_rna_type(self):
        self.api_rna_type = self.api.api("prok_rna.rna_type")
        rna_type = self.extract_biotype.output_dir + "/gene_biotype.txt"
        self.api_rna_type.add_rna_type(rna_type)

    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("prok_rna.genome_info")
        # genome_path = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option(
        #     'align_species') == 'is_ncbi' else self.option('genome_db').prop['path']
        if self.option('align_species') == 'is_ncbi':
            genome_path = self.fna
        else:
            genome_path = self.option('genome_db').prop['path']
        all_seq = ''
        with open(genome_path, 'r') as fna_r:
            for seq in fna_r.read().split('\n>'):
                all_seq += ''.join(seq.strip().split('\n')[1:])
        self.genome_size = len(all_seq)
        gc = 0
        for i in all_seq:
            if i.lower() == 'g' or i.lower() == 'c':
                gc += 1
        self.gc_content = float(gc) / self.genome_size * 100
        self.gc_content = round(self.gc_content, 2)
        # self.gc_content = (all_seq.lower().count('g') + all_seq.lower().count('c'))/self.genome_size*100
        # gc_content = int(gc_content)
        file_path = genome_path
        # if self.option('align_species') == 'is_ncbi':
        #     if self.ncbi is False:
        #         self.db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        #         col = self.db["sg_genome_db"]
        #         genome_info = col.find_one(
        #             {"genome_id": self.option("genome_id")})
        #         self.species_name = genome_info['organism_name']
        #     if self.ncbi is True:
        #
        # self.species_name = self.option('species_name')
        self.accession = self.option('genome_id')
        self.api_geno.add_genome_info(file_path=file_path, species_name=self.species_name, accession=self.accession,
                                      genome_size=self.genome_size, gc_content=self.gc_content, major=False)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("prok_rna.prok_rna_qc")
        qc_stat_before = self.qc_stat_before.output_dir
        qc_stat_after = self.qc_stat_after.output_dir
        list_txt = os.path.join(self.option('fastq_dir').prop['path'], 'list.txt')
        fq_type = self.option("fq_type").upper()
        self.api_qc.add_samples_info(qc_stat_before, qc_stat_after, list_txt, fq_type=fq_type, group=self.option('group_table').path)
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat"
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat"  # 将qc前导表加于该处
        self.api_qc.add_gragh_info(quality_stat_before, "before")
        # self.api_qc.add_bam_path(self.mapping.output_dir + "/bam")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(
            self.option("group_table").prop["path"])
        self.logger.info("group_detail为：" + str(self.group_detail))
        self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control_file").prop["path"],
                                                                        self.group_id)
        self.compare_detail = compare_detail
        self.api_qc.add_bam_path(self._sheet.output.replace('workflow_results', 'intermediate_results'))
        # map_stat = os.path.join(self.mapping.output_dir, 'stat')
        # self.api_qc.add_bowtie2_mapping_stat(map_stat)
        rrna_stat = os.path.join(self.qc_stat_after.output_dir, 'RfamStat')
        rrna_stat2 = os.path.join(self.qc_stat_after.output_dir, 'RrnaStat')
        if os.path.exists(self.rock_index.output_dir + "/rrna.fa"):
            self.api_qc.add_rfam_stat(rrna_stat, rrna_stat2, group=self.option('group_table').path)
        else:
            self.api_qc.add_rfam_stat(rrna_stat, group=self.option('group_table').path)

    def export_productive_name(self):
        api = self.api.api('prok_rna.prok_rna_qc')
        if self.option('productive_table').is_set:
            api.add_productive_name(samples=self.option('group_table').prop['sample'],
                                    productive_table=self.option('productive_table').path)
        else:
            pass

    @time_count
    def export_mapping_only(self):
        self.api_map = self.api.api("prok_rna.prok_rna_qc")
        stat_dir = self.mapping.output_dir + "/stat"
        self.api_map.add_bowtie2_mapping_stat(stat_dir)

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
        bsrd = glob.glob(os.path.join(annot_dir, '*BSRD*'))
        if len(bsrd) > 0:
            self.api_srna.add_srna_annot(annot_id, annot_dir, rfam=True, bsrd=True)
        else:
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
        db_outpath = os.path.join(self.output_dir, 'rock_index')
        db_path = self.rock_index.output_dir
        # genome_path = os.path.join(self.db_path, 'genome_fna', self.option('genome_id') + '.fna') if self.option(
        #     'align_species') == 'is_ncbi' else self.option('genome_db').prop['path']
        genome_path = self.fna if self.option(
            'align_species') == 'is_ncbi' else self.option('genome_db').prop['path']
        self.api_srna.build_seq_database(db_outpath + '/cds.fa.db.sqlite3', db_path + '/cds.fa', table_name='seq_gene',
                                         type_='gene')
        self.api_srna.build_seq_database(db_outpath + '/cds.faa.db.sqlite3', db_path + '/cds.faa',
                                         table_name='seq_protein', type_='protein')
        self.api_srna.build_seq_database(db_outpath + '/genome_fna.db.sqlite3', genome_path, table_name='seq_genome',
                                         type_='genome')
        # self.api_srna.add_

    @time_count
    def export_map_assess(self):
        gevent.sleep()
        self.api_map = self.api.api("prok_rna.prok_rna_qc")
        stat_dir = self.mapping.output_dir + "/stat"
        self.api_map.add_bowtie2_mapping_stat(stat_dir, group=self.option('group_table').path)
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
        with open(self.extract_mrna.option("out_fasta").prop["path"], 'r') as cds_r:
            cds_num = len(cds_r.read().strip().split('\n>'))
        with open(self.genome_stat, 'w') as genome_s:
            genome_s.write('Organism' + '\t' + self.species_name + '\n')
            genome_s.write('Accession' + '\t' + self.accession + '\n')
            genome_s.write('Size' + '\t' + str(self.genome_size) + '\n')
            genome_s.write('GC' + '\t' + str(self.gc_content) + '\n')
            genome_s.write('CDS' + '\t' + str(cds_num) + '\n')
        params = {
            "nr_evalue": str(self.option("nr_blast_evalue")),
            "swissprot_evalue": str(self.option("swissprot_blast_evalue")),
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
                                  task_id=task_id, method_type="samtools", project_sn=project_sn,
                                  new_output=self.work_dir + "/snp_tmp", group=self.option('group_table').path)

    @time_count
    def export_expression(self):
        gevent.sleep()
        self.stop_timeout_check()
        all_exp = self.api.api("prok_rna.all_exp")
        group_dict = self.option('group_table').prop['group_dict']
        ### 此处和邱于涵商量好存储字段大小写格式
        # 工作流跑的默认存到mongo库params的method是RSEM，salmon，kallisto；exp_type是TPM，FPKM
        # 交互分析运行后存到mongo库params的method是RSEM，salmon，kallisto；exp_type是TPM，FPKM
        if self.option('express_method').lower() == "rsem":
            quant_method = self.option('express_method').upper()
        else:
            quant_method = self.option('express_method').lower()
        group_id = self.group_id
        control_id = self.control_id
        task_id = self.task_id
        project_sn = self.project_sn

        # add exp matrix
        ## 还需确认路径信息，因为后续交互需要用到express.output_dir
        exp_output = self.express.output_dir
        # self.exp_level = 'mRNA+sRNA' if self.option('strand_specific') else 'mRNA'
        self.exp_level = 'ref_gene'
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
                exp_matrix, quant_method=quant_method, group_dict=group_dict, main_id_m=None, main_id_s=None,
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
                                     quant_method=quant_method, project_sn=project_sn, task_id=task_id,
                                     exp_level=self.exp_level)
            self.exp_id = trans_exp_id
        else:
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, group_dict=group_dict, main_id_m=None,
                                           main_id_s=None,
                                           main_id_ms=None, group_id=group_id,
                                           add_distribution=True, exp_type='FPKM', project_sn=project_sn,
                                           task_id=task_id,
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
                                     quant_method=quant_method, project_sn=project_sn, task_id=task_id,
                                     exp_level=self.exp_level)
            self.exp_id = trans_exp_id

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
            Draw_in_groups='no',
            log_base='10'
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
            threshold="1",
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
            main_id = all_exp.add_exp_pca2(pca_output, exp_level=self.exp_level,
                                           params=params, project_sn=project_sn, task_id=task_id)
            if self.min_group_num >= 3:
                all_exp.insert_ellipse_table(os.path.join(self.ellipse.work_dir, 'ellipse_out.xls'), main_id)

        # add gene diff
        diff_output = self.diffexpress.tool.output_dir
        uniform_output = self.diffexpress.uniform.output_dir
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
            # correct_method=self.option('padjust_way'),
            # stat_type=stat_type,
            # stat_cutoff=self.option('diff_fdr_ci'),
            # quant_method=quant_method,
            diff_method=diff_method,
            is_batch='False'
        )
        if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            params.update({
                'stat_type': stat_type,
                'stat_cutoff': self.option('diff_fdr_ci')
            })
            if stat_type == "padjust":
                params.update({
                    'correct_method': self.option('padjust_way')
                })
        if diff_method.lower() in ['noiseq']:
            params.update({
                'prob': self.option('diff_fdr_ci')
            })
        all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id, exp_matrix=exp_matrix,
                                group_dict=group_dict, group_id=group_id,
                                exp_level=exp_level, quant_method=quant_method,
                                diff_method=diff_method,
                                project_sn=project_sn, task_id=task_id, params=params,
                                pvalue_padjust=stat_type
                                )
        # all_exp.add_diffexp(diff_output, exp_id=exp_id, group_dict=group_dict, group_id=group_id,
        #                     exp_level=exp_level, quant_method=quant_method, diff_method=diff_method,
        #                     project_sn=project_sn, task_id=task_id, params=params,
        #                     pvalue_padjust=stat_type)

    def export_structure(self, srna=False, promote=False, terminator=False, antisense=False, novelgene=False):
        structure = self.api.api("prok_rna.gene_structure")
        if srna:
            structure.add_rock_structure(self.rockhopper.output_dir)
        if promote:
            structure.add_promote_structure(self.promote.output_dir + '/all_promoter_result.xls')
        if terminator:
            structure.add_terminator_structure(self.terminator.output_dir + '/terminator_prediction.xls')
        if antisense:
            structure.add_antisense(self.rockhopper.output_dir + '/antisense.xls')
        if novelgene:
            structure.add_novelgene(self.rockhopper.output_dir + '/genome.predicted_cds.xls')

    def export_genesets_analysis(self):
        if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir,"results_info")):
            pass
        else:
            if os.path.exists(os.path.join(self.work_dir, "temporary")):
                shutil.rmtree(os.path.join(self.work_dir, "temporary"))
            os.makedirs(os.path.join(self.work_dir, "temporary"))
            self.export_temporary = os.path.join(self.work_dir, "temporary")
            api = self.api.api('prok_rna.diff_geneset_work_pipline')
            diff_geneset_pipline_result = self.diff_geneset_analysis.output_dir
            diff_id = None
            task_id = self.task_id
            analysis_names = ["kegg", "go", "cog"]
            file_json_path = os.path.join(self.diff_geneset_analysis.file_prepare.output_dir, "prepare_json")
            with open(file_json_path, "r") as j:
                file_dict = json.load(j)
            geneset_file = file_dict["genesets"]["All_Diff_mRNA"]['file_path']["kegg_class"]['multi_gene_list_path']
            kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
            task_infos = {"task_id":self.task_id,"project_sn":self.project_sn}
            api.add_diff_genest_pipline_table(diff_geneset_pipline_result, diff_id=diff_id, task_id=task_id,
                                              analysis_names=analysis_names,
                                              kegg_level_path=kegg_level_path, inter_path=self.export_temporary,
                                              group_id=self.group_id,task_infos = task_infos,geneset_file = geneset_file,exp_id=self.exp_id)

    def export_report_img(self):
        report_config = os.path.join(self.chart.work_dir, 'report_config.json')
        api =  self.api.api('prok_rna.report_model')
        s3 = self._sheet.output.split(":")[0]
        report_img_s3 = s3 + ":commonbucket/files/report_img/prokrna/" + self.task_id
        api.add_report_image(self.task_id, report_config, report_img_s3)

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

        if os.path.exists(os.path.join(exp_output, 'transcript.tpm.matrix')):
            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            trans_pd = pd.read_table(exp_matrix, header=0, index_col=0)
            trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
            trans_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            trans_result.to_csv(trans_out, header=True, index=True, sep='\t')
        if os.path.exists(os.path.join(exp_output, 'transcript.fpkm.matrix')):
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            trans_pd = pd.read_table(exp_matrix, header=0, index_col=0)
            trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
            trans_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
            trans_result.to_csv(trans_out, header=True, index=True, sep='\t')

    def merge_annotation_exp_matrix_total(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        exp_output = self.express_total.output_dir
        annot = os.path.join(self.output_dir, 'annot_class/anno_stat/all_anno_detail.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        trans_annot_pd = all_annot.reset_index().set_index('gene_id')

        if os.path.exists(os.path.join(exp_output, 'transcript.tpm.matrix')):
            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            trans_pd = pd.read_table(exp_matrix, header=0, index_col=0)
            trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
            trans_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            trans_result.to_csv(trans_out, header=True, index=True, sep='\t')
        if os.path.exists(os.path.join(exp_output, 'transcript.fpkm.matrix')):
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            trans_pd = pd.read_table(exp_matrix, header=0, index_col=0)
            trans_result = pd.concat([trans_pd, trans_annot_pd], axis=1)
            trans_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
            trans_result.to_csv(trans_out, header=True, index=True, sep='\t')

    def merge_annotation_diffexp_matrix(self):
        """
        给表达矩阵添加注释信息
        :return:
        """
        annot = os.path.join(self.output_dir, 'annot_class/anno_stat/all_anno_detail.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        diff_output = self.diffexpress.output_dir

        duplicate_files = glob.glob(diff_output + '/' + '*_vs_*.annot.xls')  ## 为了防止流程重跑的时候反复增加注释结果表
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
                                        df_anno_pre1['chrom'].apply(str) + df_anno_pre1["end"].apply(str) + \
                                        df_anno_pre1["start"].apply(str) + df_anno_pre1["ref"].apply(str)
        df_anno_pre2 = pd.read_table(y, header=0, sep="\t", low_memory=False)
        df_anno_pre2_select_pre = df_anno_pre2.loc[:, df_anno_pre2.columns[:13]]
        df_anno_pre2_select_pre.rename(
            columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno", "END": "End",
                     "START": "Start", "REF": "Ref",
                     "MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
        df_anno_pre2_select = df_anno_pre2_select_pre[
            ["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
             "Total depth", "QUAL",
             "Anno", "MUT type", "MUT info"]]
        df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2[
            'CHROM'].apply(str) + \
                                        df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + \
                                        df_anno_pre2["REF"].apply(str)
        df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
        df_join_pre.drop(columns=['index1'], inplace=True)
        df_join_pre.to_csv(self.work_dir + "/upload_results/GeneStructure/SNP/snp_annotation_detail.xls", sep="\t",
                           index=False)


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
        # worker = worker_client()
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
                express_method='rsem',
                exp_way='tpm'
            )
        }

        info = worker.add_task(data)
        print(info)


if __name__ == '__main__':
    unittest.main()
