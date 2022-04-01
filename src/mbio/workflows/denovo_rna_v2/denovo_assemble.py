# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""无参转录组工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import json
import shutil
import re
import time
#from gevent.monkey import patch_all
import gevent
import functools
from bson.son import SON
from biocluster.config import Config
import pandas as pd
from collections import OrderedDict
from mbio.packages.denovo_rna_v2.nr_stat import nr_stat




#patch_all()


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


class DenovoAssembleWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(DenovoAssembleWorkflow, self).__init__(wsheet_object)
        options = [
            ## 单样本 OR 多样本
            {"name": "sample_num", "type": "string", "default": "multiple"},  # 测序类型，single OR multiple
            ##基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "datatype", "type": "string", "default": "rawdata"},
            {"name": "qc_soft", "type" : "string", "default": "fastp"},
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "strand_specific", "type": "bool", "default": False},  # 是否有链特异性, 默认是False, 无特异性
            {"name": "strand_dir", "type": "string", "default": ""},  # 当链特异性时为True时，正义链为forward，反义链为reverse
            {"name": "is_duplicate", "type": "bool", "default": ""},  # 是否有生物学重复
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹
            {"name": "qc_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹
            {"name": "group_table", "type": "infile", "format": "denovo_rna_v2.group_table"},  # 分组文件
            {"name": "control", "type": "infile", "format": "denovo_rna_v2.compare_table"},  # 对照文件
            # 从头组装
            {"name": "assemble", "type": "bool", "default": True},  #是否进行从头组装
            {"name": "assemble_soft", "type": "string", "default": "trinity"},  # 选择拼接软件trinity, spades
            {"name": "assemble_method", "type": "string", "default": "total"}, # 选择拼接软件sample, group, total

            {'name': 'level', 'type': 'string', 'default': 'Transcript'},
            {"name": "k_mer", "type": "int", "default": 25},  # 设置Trinity软件的kmer值
            {"name": "min_kmer_cov", "type": "int", "default": 5},  # 设置Trinity软件的kmer值
            {"name": "normalize", "type": "bool", "default": True},  # 是否进行reads均一化
            {"name": "normalize_max_read_cov", "type": "int", "default": 50},  # reads均一化最大覆盖倍数
            {"name": "jaccard_clip", "type": "bool", "default": False},  #分割高密度基因区间基因，建议物种为Fungi时True，其他物种时False
            {"name": "min_contig_length", "type": "int", "default": 200},  #设置最短contig的长度
            {"name": "remove_samples", "type": "infile", 'format': "denovo_rna_v2.common"},  #不参与组装的样本
            {"name": "optimize", "type": "bool", "default": True},  #是否进行组装结果优化
            {"name": "transrate", "type": "bool", "default": True},  #是否选择TransRate进行组装结果评估
            {"name": "cd_hit", "type": "float", "default": 0.99},  # 利用cd-hit软件进行过滤时identity的取值，值为1时不进行过滤
            {"name": "tpm", "type": "float", "default": 0},  #过滤掉tpm小于给定值的转录本

            # spades 运行参数
            {'type': 'string', 'name': 'spades_k', 'default': "auto"},
            {'type': 'string', 'name': 'spades_c', 'default': "auto"},

            {'type': 'string', 'name': 'merge_soft', 'default': "tgicl"},

            {'type': 'float', 'name': 'tgicl_p', 'default': 88},
            {'type': 'int', 'name': 'tgicl_c', 'default': 10},

            # 注释
            {"name": "nr_evalue", "type": "string", "default": "1e-5"},
            {"name": "string_evalue", "type": "string", "default": "1e-5"},
            {"name": "kegg_evalue", "type": "string", "default": "1e-5"},
            {"name": "swissprot_evalue", "type": "string", "default": "1e-5"},
            {"name": "pfam_evalue", "type": "string", "default": "1e-5"},
            {"name": "tf_evalue", "type": "string", "default": "1e-5"},
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            {"name": "nr_sub_database", "type": "string", "default": None},  # nr库类型
            {"name": "kegg_database", "type": "string", "default": "All"},  # kegg注释库类型
            {"name": "tf_database", "type": "string", "default": ""},  # tf_database注释库类型,Plant,Animal,Other
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "tf_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg,swissprot,pfam'},
            # 表达及差异表达
            {"name": "express_method", "type": "string", "default": "Salmon"},  #选择表达量计算软件，RSEM,Kallisto,Salmon，定量指标为TPM
            # {"name": "fastp", "type": "int", "default": 0},  #是否选用fastp进行质控
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"},
            # 基因和转录本对应关系文件
            {"name": "gene_to_trans", "type": "infile", 'format': "denovo_rna_v2.common"},

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
            self.set_error("json output wrong", code="12002801")

        self.project_sn = self._sheet.project_sn #获取project_sn
        self.task_id = self._sheet.id #获取task_id
        self.add_option(options)
        self.set_options(self._sheet.options())

        #添加tool/module
        self.filecheck = self.add_tool("denovo_rna_v2.filecheck_denovo")
        if self.option('qc_soft') == "fastp":
            self.qc = self.add_module("denovo_rna_v2.fastp_rna")
        else:
            self.qc = self.add_module("denovo_rna_v2.hiseq_qc")
        self.qc_stat_before = self.add_module("lnc_rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("lnc_rna.hiseq_reads_stat")
        self.assemble = self.add_module("denovo_rna_v2.denovo_assemble3")
        self.assemble_filter = self.add_module("denovo_rna_v2.denovo_assemble2_filter")
        self.filter_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")
        self.unigene_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")
        self.filter_unigene_evaluation = self.add_module("denovo_rna_v2.assemble_evalution")

        #self.diamond = self.add_module("denovo_rna_v2.diamond")
        #self.blast = self.add_module("denovo_rna_v2.blast")
        self.annotation = self.add_module("denovo_rna_v2.denovo_annotation")
        if self.option("express_method") == "RSEM":
            self.align = self.add_module("denovo_rna_v2.quant")
        else:
            self.align = self.add_module("denovo_rna_v2.quant")
            # self.express = self.add_module("denovo_rna_v2.quant")

        self.cds_predict = self.add_module("denovo_rna_v2.cds_tf")
        self.nr_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.swissprot_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.cog_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.kegg_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.pfam_filter = self.add_tool("denovo_rna_v2.filter_annot")

        #判断流程结束tool/module list
        if self.option("sample_num") == "multiple":
            self.final_tools = [self.annotation,
                                self.cds_predict,
                                self.align]
        else:
            self.final_tools = [self.annotation, self.cds_predict, self.align]

        if self.option("assemble_soft").lower() == "trinity" and self.option("assemble_method") == "total":
            self.final_tools.append(self.unigene_evaluation)
            self.final_tools.append(self.filter_unigene_evaluation)
        '''
        if self.option("express_method").lower() != "rsem":
            self.final_tools.append(self.express)
        '''
        if self.option("optimize") == True:
            self.final_tools.append(self.filter_evaluation)

        # 添加step，显示在页面进度条
        self.step.add_steps( "qc_stat_after", "qc_stat_before")
        if self.option("sample_num") == "multiple":
            self.step.add_steps("filecheck", "qualitycontrol", "assemble",
                                "evaluation", "annotation", "express",
                                "cds_predict")
        else:
            if self.option("optimize") == True:
                self.step.add_steps("filecheck", "qualitycontrol", "assemble",
                                    "evaluation", "annotation", "express",
                                    "cds_predict", "unigene_evaluation", "filter_unigene_evaluation")
            else:
                self.step.add_steps("filecheck", "qualitycontrol", "assemble", "annotation", "express", "cds_predict")

        if self.option("assemble_soft").lower() == "trinity" and self.option("assemble_method") == "total":
            self.step.add_steps("unigene_evaluation", "filter_unigene_evaluation")

    def check_options(self):
        """
        检查选项
        """
        # 基础参数
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE", code = "12002801")

        # 组装相关参数
        if not self.option("k_mer") >= 20 and not self.option("string_evalue") <= 32:
            raise OptionError("kmer超出范围，值的范围为[25-32]", code = "12002802")

        # 注释相关参数
        try:
            nr_evalue = float(self.option("nr_evalue"))
            string_evalue = float(self.option("string_evalue"))
            kegg_evalue = float(self.option("string_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
            tf_evalue = float(self.option("tf_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范", code = "12002803")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("string_blast_evalue", string_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)

        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围", code = "12002804")
        if not self.option("string_blast_evalue") > 0 and not self.option("string_blast_evalue") < 1:
            raise OptionError("String比对的E值超出范围", code = "12002805")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围", code = "12002806")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围", code = "12002807")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围", code = "12002808")


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

    def run_filecheck(self):
        self.logger.info("开始运行文件检查")
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'group_table': self.option('group_table'),
            'sample_num': self.option('sample_num'),
        }
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        self.logger.info("开始运行质控")
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
                'fq_type': self.option('fq_type')
            })

        self.qc.on('start', self.set_step, {'start': self.step.qualitycontrol})
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('end', self.set_step, {'end': self.step.qualitycontrol})
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
                    'rfam': False,
                    'dup': True
                })
            else:
                self.qc_stat_after.set_options({
                    'fastq_dir': self.qc.option('sickle_dir'),
                    'fq_type': self.option('fq_type'),
                    'quality': quality,
                    'rfam': False,
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

    def run_assemble(self):
        self.logger.info("开始运行组装")
        opts = {
            "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "fq_type" : self.option("fq_type"),
            "strand_direct" : self.option("strand_dir"),
            "group": self.option("group_table"),
            "assemble_soft": self.option("assemble_soft").lower(),
            "assemble_method": self.option("assemble_method")
        }
        if self.option("assemble_soft").lower() == "trinity":
            opts.update({
                "filter_sample" : self.option("remove_samples"),
                "min_contig_length" : self.option("min_contig_length"),
                "kmer_size" : self.option("k_mer"),
                "min_kmer_cov" : self.option("min_kmer_cov"),
                "jaccard_clip" : self.option("jaccard_clip"),
                "no_normalize_reads" : not self.option("normalize"),
                "normalize_max_read_cov" : self.option("normalize_max_read_cov"),
            })
        else:
            opts.update({
                "spades_k" : self.option("spades_k"),
            })

        if self.option("assemble_method") != "total":
            opts.update({
                "tgicl_p" : self.option("tgicl_p"),
                "tgicl_c" : self.option("tgicl_c"),
            })
        self.assemble.set_options(opts)
        self.assemble.on("end", self.set_output, "assemble")
        self.assemble.on('start', self.set_step, {'start': self.step.assemble})
        self.assemble.on('end', self.set_step, {'end': self.step.assemble})
        self.assemble.run()


    def run_assemble_filter(self):
        if self.option("assemble") == True:
            opts = {
                "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
                "fq_type": self.option("fq_type"),
                "assemble_fa": self.assemble.output_dir + '/assemble_raw.fasta',
                "assemble_g2t": self.assemble.output_dir + '/assemble_raw.gene_trans_map',
                "min_contig_length": self.option("min_contig_length"),
                "species" : self.option("kegg_database"),
                "transrate_filter" : self.option("transrate"),
                "cdhit_filter" : True,
                "TPM_filter" : False,
                "cdhit_identity_threshold" : self.option("cd_hit"),
                "TPM_threshold" : self.option("tpm"),
            }
        else:
            opts = {
                "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
                "fq_type": self.option("fq_type"),
                "assemble_fa": self.option("assembly_file").prop['path'],
                "assemble_g2t": self.option("gene_to_trans").prop['path'],
                "species" : self.option("kegg_database"),
                "transrate_filter" : self.option("transrate"),
                "cdhit_filter" : True,
                "TPM_filter" : False,
                "cdhit_identity_threshold" : self.option("cd_hit"),
                "TPM_threshold" : self.option("tpm"),
            }

        self.assemble_filter.set_options(opts)
        self.assemble_filter.on("end", self.set_output, "assemble")
        self.assemble_filter.on('start', self.set_step, {'start': self.step.assemble})
        self.assemble_filter.on('end', self.set_step, {'end': self.step.assemble})
        self.assemble_filter.run()

    def run_filter_evaluation(self):
        self.logger.info("开始运行组装优化后评估")
        opts = {
            "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "fq_type" : self.option("fq_type"),
            "min_contig_length": self.option("min_contig_length"),
            "species" : self.option("kegg_database"),
            "assemble_fa" : self.assemble_filter.option("filter_fa"),
            "assemble_g2t" : os.path.join(self.assemble_filter.output_dir,"Trinity.filter.gene_trans_map"),
        }
        self.filter_evaluation.set_options(opts)
        self.filter_evaluation.on("end", self.set_output, "filter_evaluation")
        self.filter_evaluation.on('start', self.set_step, {'start': self.step.evaluation})
        self.filter_evaluation.on('end', self.set_step, {'end': self.step.evaluation})
        self.filter_evaluation.run()

    def run_unigene_evaluation(self):
        self.logger.info("开始运行组装优化后评估")
        seq = glob.glob(self.assemble_filter.output_dir + '/trinity_stat/unigene.fasta')[0]
        opts = {
            "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "fq_type" : self.option("fq_type"),
            "species" : self.option("kegg_database"),
            "min_contig_length": self.option("min_contig_length"),
            "assemble_fa" : seq,
        }
        self.unigene_evaluation.set_options(opts)
        self.unigene_evaluation.on("end", self.set_output, "unigene_evaluation")
        self.unigene_evaluation.on('start', self.set_step, {'start': self.step.unigene_evaluation})
        self.unigene_evaluation.on('end', self.set_step, {'end': self.step.unigene_evaluation})
        self.unigene_evaluation.run()

    def run_filter_unigene_evaluation(self):
        self.logger.info("开始运行组装优化后评估")
        seq = glob.glob(self.assemble_filter.output_dir + '/Trinity.filter.unigene.fasta')[0]
        opts = {
            "sample_fq_list": os.path.join(self.qc.option("sickle_dir").prop['path'], "fq_list.txt"),
            "fq_type" : self.option("fq_type"),
            "min_contig_length": self.option("min_contig_length"),
            "species" : self.option("kegg_database"),
            "assemble_fa" : seq,
        }
        self.filter_unigene_evaluation.set_options(opts)
        self.filter_unigene_evaluation.on("end", self.set_output, "filter_unigene_evaluation")
        self.filter_unigene_evaluation.on('start', self.set_step, {'start': self.step.filter_unigene_evaluation})
        self.filter_unigene_evaluation.on('end', self.set_step, {'end': self.step.filter_unigene_evaluation})
        self.filter_unigene_evaluation.run()

    def run_nr_stat(self):
        nr_stat().nr_stat_species(
            tran_taxonfile = self.annotation.output_dir + '/blast_nr_taxon/query_taxons_detail.xls',
            gene_list = "",
            exp = self.align.output_dir + '/transcript.count.matrix',
            outpath = self.work_dir + '/transcript_species_count.xls'
        )

        nr_stat().nr_stat_species(
            tran_taxonfile = self.annotation.output_dir + '/blast_nr_taxon/query_taxons_detail.xls',
            gene_list = self.assemble_filter.output_dir + '/Trinity.filter.gene_trans_map',
            exp = self.align.output_dir + '/gene.count.matrix',
            outpath = self.work_dir + '/gene_species_count.xls'
        )


    def run_align(self):
        self.logger.info("开始运行样本比对和表达定量分析")
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
        opts = {
            "fastq" : self.qc.option("fq_list"),
            "method" : self.option("express_method"),
            "libtype" : self.libtype,
            "transcriptome" : self.assemble_filter.option("filter_fa"),
            "t2g" : os.path.join(self.assemble_filter.output_dir,"Trinity.filter.gene_trans_map"),
        }
        self.align.set_options(opts)
        self.align.on("end", self.set_output, "align")
        self.align.on('start', self.set_step, {'start': self.step.express})
        self.align.on('end', self.set_step, {'end': self.step.express})
        self.align.run()


    def run_cds_predict(self):
        self.logger.info("开始运行cds预测")
        opts = {
            "fasta" : self.assemble_filter.option("filter_fa"),
            "e_value" : 1e-3,
            "species_type" : self.option("tf_database"),
            "isoform_unigene" : os.path.join(self.assemble_filter.output_dir,"Trinity.filter_t2g2u"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.set_output, "cds_predict")
        self.cds_predict.on("end", self.run_pfam_filter)
        self.cds_predict.on('start', self.set_step, {'start': self.step.cds_predict})
        self.cds_predict.on('end', self.set_step, {'end': self.step.cds_predict})
        self.cds_predict.run()

    def run_diamond(self):
        self.logger.info("开始运行blast注释")
        self.blast_modules = []
        #blast_lines = int(self.new_trans_abs.option('query').prop['seq_number']) / 10 + 1
        #self.logger.info('.......blast_lines:%s' % blast_lines)
        blast_opts = {
            'query': self.assemble_filter.option("filter_fa"),
            'query_type': 'nucl',
            'database': None,
            'blast': 'blastx',
            'evalue': None,
            'outfmt': 5,
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.diamond_nr = self.add_module("denovo_rna_v2.diamond")
            if self.option("nr_sub_database") in ["Spermatophyta", "Reptilia", "OtherPlant", "Mammalia", "Invertebrate", "Fishes", "Aves", "Amphibia", "Algae"]:
                blast_opts.update(
                    {
                        'database': self.option("nr_sub_database"),
                        'evalue': 1e-3
                    }
                )
            else:
                blast_opts.update(
                    {
                        'database': self.option("nr_database"),
                        'evalue': 1e-3
                    }
                )
            self.diamond_nr.set_options(blast_opts)
            self.blast_modules.append(self.diamond_nr)
            self.diamond_nr.on('end', self.set_output, 'diamond_nr')
            self.diamond_nr.on('end', self.run_nr_blast_filter)
            self.diamond_nr.run()
        if 'cog' in self.option('database'):
            self.diamond_string = self.add_module("denovo_rna_v2.diamond")
            blast_opts.update(
                {'database': 'eggnog', 'evalue': 1e-3}
            )
            self.diamond_string.set_options(blast_opts)
            self.blast_modules.append(self.diamond_string)
            self.diamond_string.on('end', self.set_output, 'diamond_string')
            self.diamond_string.on('end', self.run_cog_blast_filter)
            self.diamond_string.run()
        if 'kegg' in self.option('database'):
            self.diamond_kegg = self.add_module("denovo_rna_v2.diamond")
            blast_opts.update(
                {'database': 'kegg', 'evalue': 1e-3}
            )
            self.diamond_kegg.set_options(blast_opts)
            self.blast_modules.append(self.diamond_kegg)
            self.diamond_kegg.on('end', self.set_output, 'diamond_kegg')
            self.diamond_kegg.on('end', self.run_kegg_blast_filter)
            self.diamond_kegg.run()
        if 'swissprot' in self.option('database'):
            self.blast_swissprot = self.add_module('denovo_rna_v2.blast')
            blast_opts.update(
                {'database': 'swissprot', 'evalue': 1e-3}
            )
            self.blast_swissprot.set_options(blast_opts)
            self.blast_modules.append(self.blast_swissprot)
            self.blast_swissprot.on('end', self.set_output, 'blast_swissprot')
            self.blast_swissprot.on('end', self.run_swissprot_blast_filter)
            self.blast_swissprot.run()
        # self.on_rely([self.diamond_nr, self.diamond_kegg, self.diamond_string, self.blast_swissprot, self.cds_predict], self.run_annotation)

    def run_nr_blast_filter(self):
        options = {
            'xml': self.diamond_nr.option('outxml').prop['path'],
            'types': "xml",
            'evalue': self.option('nr_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.nr_filter.set_options(options)
        self.nr_filter.run()

    def run_swissprot_blast_filter(self):
        options = {
            'xml': self.blast_swissprot.option('outxml').prop['path'],
            'types': "xml",
            'evalue': self.option('swissprot_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.swissprot_filter.set_options(options)
        self.swissprot_filter.run()

    def run_cog_blast_filter(self):
        options = {
            'xml': self.diamond_string.option('outxml').prop['path'],
            'types': "xml",
            'evalue': self.option('string_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.cog_filter.set_options(options)
        self.cog_filter.run()

    def run_kegg_blast_filter(self):
        options = {
            'xml': self.diamond_kegg.option('outxml').prop['path'],
            'types': "xml",
            'evalue': self.option('kegg_evalue'),
            'identity': 0,
            'similarity': 0
        }
        self.kegg_filter.set_options(options)
        self.kegg_filter.run()

    def run_pfam_filter(self):
        options = {
            'hmm': self.cds_predict.output_dir + "/pfam_domain",
            'types': "hmm",
            'evalue': self.option('pfam_evalue'),
        }
        self.pfam_filter.set_options(options)
        self.pfam_filter.run()

    def run_annotation(self):
        self.logger.info("开始运行注释统计")
        anno_opts = {
            "gene2trans" : os.path.join(self.assemble_filter.output_dir,"Trinity.filter.gene_trans_map"),
            "go_annot" : True,
            "nr_annot" : True,
            "taxonomy" : self.option("kegg_database"),
            "blast_nr_xml" : self.nr_filter.option('outxml').prop['path'],
            "blast_kegg_xml" : self.kegg_filter.option('outxml').prop['path'],
            "blast_eggnog_xml" : self.cog_filter.option('outxml').prop['path'],
            "blast_swissprot_xml" : self.swissprot_filter.option('outxml').prop['path'],
            "pfam_domain" : self.pfam_filter.option('outtable').prop['path']
        }
        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        self.annotation.run()


    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "12002802")
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
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'annotation')
        if event['data'] == 'assemble':
            self.move2outputdir(obj.output_dir, 'assemble')
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

    def run(self):
        """
        denovo-rna workflow run方法
        """
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False

        if self.option("fastq_dir").is_set:
            self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计

        if self.option("qc_dir").is_set:
            self.filecheck.on('end', self.run_assemble )  # 质控前统计
            self.filecheck.on('end', self.run_qc_stat, True)  # 质控后统计
            if self.option("assemble") == True:
                self.filecheck.on('end', self.run_assemble)
            else:
                self.filecheck.on('end', self.run_assemble_filter)

        else:
            self.filecheck.on('end', self.run_qc)
            self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
            if self.option("assemble") == True:
                self.qc.on('end', self.run_assemble)
            else:
                self.qc.on('end', self.run_assemble_filter)

        # if self.option('fastp'):
        #     self.qc_stat_before.on('end', self.run_qc)
        # else:
        #     self.filecheck.on('end', self.run_qc)
        if self.option("assemble") == True:
            self.assemble.on('end', self.run_assemble_filter)
            
        if self.option("optimize") == True:
            self.assemble_filter.on('end', self.run_filter_evaluation)
        else:
            pass
        if self.option("express_method") == "RSEM":
            self.assemble_filter.on('end', self.run_align)
        else:
            self.assemble_filter.on('end', self.run_align)

        self.assemble_filter.on('end', self.run_diamond)
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

        self.assemble_filter.on('end', self.run_cds_predict)
        if self.option("assemble_soft").lower() == "trinity" and self.option("assemble_method") == "total":
            self.assemble_filter.on('end', self.run_filter_unigene_evaluation)
            self.assemble_filter.on('end', self.run_unigene_evaluation)
        self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(DenovoAssembleWorkflow, self).run()

    def end(self):
        self.run_nr_stat()
        self.run_api() # 运行导表函数

        ## 导表后，修改文件的绝对路径
        db = Config().get_mongo_client(mtype="denovo_rna_v2")[Config().get_mongo_dbname("denovo_rna_v2")]
        col1 = db["sg_annotation_stat"]
        col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_task"]
        col2.update({"task_id" : self.task_id}, {"$set": {"seq_db": self.workflow_output + "/Sequence_database/seq_db.sqlite3"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"fastq": self.workflow_output + "/QC/fq_list.txt"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.workflow_output + "/Assemble/AssembleResult/Trinity.filter.fasta"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"assemble_t2g": self.workflow_output + "/Assemble/AssembleResult/Trinity.filter.gene_trans_map"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"bedpath": self.workflow_output + "/CDS/Trinity.filter.fasta.transdecoder.bed"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"unigene_fa": self.workflow_output + "/Assemble/AssembleResult/Trinity.filter.unigene.fasta"}}, upsert=True)
        col3 = db["sg_exp"]
        col3.update({"task_id" : self.task_id, "exp_level" : "G"}, {"$set": {"count_file": self.workflow_output + "/Express/ExpAnnalysis/unigene.count.matrix.xls"}}, upsert=True)
        col3.update({"task_id" : self.task_id, "exp_level" : "T"}, {"$set": {"count_file": self.workflow_output + "/Express/ExpAnnalysis/transcript.count.matrix.xls"}}, upsert=True)
        if self.option("sample_num") == "multiple":
            col4 = db["sg_snp"]
            col4.update({"task_id" : self.task_id}, {"$set": {"call_vcf_path": self.workflow_output + "/SNP/call.vcf"}}, upsert=True)

        self.modify_output() # 修改文件目录结构
        super(DenovoAssembleWorkflow, self).end()

    def modify_output(self):
        # 文件上传到页面后，需要修改文件的内容，才能用于交互分析
        # old_name = self.work_dir + "/HiseqQc/output/sickle_dir/fq_list.txt"
        old_name = self.qc.option('fq_list').prop['path']
        a = open(old_name,"r").read()
        start_dir = self.work_dir + "/HiseqQc/output/sickle_dir/"
        end_dir = self.workflow_output + "/QC/cleandata/"
        b = a.replace(start_dir, end_dir)
        new_name = self.work_dir + "/HiseqQc/output/sickle_dir/tmp.txt"
        if self.option('qc_soft') == "fastp":
            new_name = self.work_dir + "/FastpRna/output/fastq//tmp.txt"
        with open(new_name, "w") as f:
            f.write(b)
        os.remove(old_name)
        os.rename(new_name, old_name)

        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"
        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        fq_list = origin_dir + "/QC_stat/sickle_dir/fq_list.txt"
        if self.option('qc_soft') == "fastp":
            fq_list = origin_dir + "/QC_stat/fastq/fq_list.txt"
        os.mkdir(target_dir + "/QC")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        os.link(fq_list, target_dir + "/QC/fq_list.txt")
        # Assemble
        os.mkdir(target_dir + "/Assemble")
        os.mkdir(target_dir + "/Assemble/AssembleResult")
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
            optimize_after_stat = os.path.join(origin_dir + "/filter_evaluation/trinity_stat/Trinity_stat.xls")
            os.link(optimize_after_stat, target_dir + "/Assemble/AssembleEvaluation/Optimize_after/Trinity_stat.xls")
        os.mkdir(target_dir + "/Assemble/LengthDistribution")
        if self.option("optimize") == False:
            gene_len_dis = os.path.join(origin_dir + "/assemble/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, target_dir + "/Assemble/LengthDistribution/unigene_count_stat_500.txt")
            trans_len_dis = os.path.join(origin_dir + "/assemble/trinity_stat/trans_count_stat_500.txt")
            os.link(trans_len_dis, target_dir + "/Assemble/LengthDistribution/transcript_count_stat_500.txt")
        else:
            gene_len_dis = os.path.join(origin_dir + "/filter_evaluation/trinity_stat/unigene_count_stat_500.txt")
            os.link(gene_len_dis, target_dir + "/Assemble/LengthDistribution/unigene_count_stat_500.txt")
            trans_len_dis = os.path.join(origin_dir + "/filter_evaluation/trinity_stat/trans_count_stat_500.txt")
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
        '''
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
        '''
        # AnnoQuery
        os.mkdir(target_dir + "/AnnoQuery")
        gene_anno_detail = os.path.join(origin_dir + "/annotation/anno_stat/gene_anno_detail.xls")
        os.link(gene_anno_detail, target_dir + "/AnnoQuery/unigene_anno_detail.xls")
        trans_anno_detail = os.path.join(origin_dir + "/annotation/anno_stat/trans_anno_detail.xls")
        os.link(trans_anno_detail, target_dir + "/AnnoQuery/transcript_anno_detail.xls")
        # Express
        '''
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
            gene_count = os.path.join(self.align.output_dir + "/gene.count.matrix")
            os.link(gene_count, target_dir + "/Express/ExpAnnalysis/unigene.count.matrix.xls")
            gene_tpm = os.path.join(self.align.output_dir + "/gene.tpm.matrix")
            os.link(gene_tpm, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.xls")
            trans_count = os.path.join(self.align.output_dir + "/transcript.count.matrix")
            os.link(trans_count, target_dir + "/Express/ExpAnnalysis/transcript.count.matrix.xls")
            trans_tpm = os.path.join(self.align.output_dir + "/transcript.tpm.matrix")
            os.link(trans_tpm, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.xls")
            gene_tpm_anno = os.path.join(self.align.output_dir + "/gene.tpm.matrix.annot.xls")
            os.link(gene_tpm_anno, target_dir + "/Express/ExpAnnalysis/unigene.tpm.matrix.annot.xls")
            trans_tpm_anno = os.path.join(self.align.output_dir + "/transcript.tpm.matrix.annot.xls")
            os.link(trans_tpm_anno, target_dir + "/Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls")

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
        '''
        repaths = [
            [".", "", "流程分析结果目录", 0, '201001'],
            ["Sequence_database", "", "序列文件数据库", 0, '201002'],
            ["seq_db.sqlite3 ", "", "", 1, '201003'],
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
            ["Annotation/Swiss-Prot/transcript_swissprot_anno_detail.xls", "", "transcript SwissProt注释详情表", 0, '201040'],
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
            #["DiffExpress/DiffExpress_transcript", "", "transcript差异表达分析结果", 0, '201093'],  20190517取消
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
        #greenlets_list_third = []
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.denovoass_task_info')
        task_info.add_task_info()
        self.logger.info("进行第一阶段导表")
        self.stop_timeout_check()
        greenlets_list_first.append(gevent.spawn(self.export_qc))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_assembly))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_align))
        gevent.joinall(greenlets_list_first)
        self.logger.info("进行第二阶段导表")
        greenlets_list_sec.append(gevent.spawn(self.export_denovo_annotation))
        gevent.joinall(greenlets_list_sec)
        self.logger.info("进行第三阶段导表")
        #greenlets_list_third.append(gevent.spawn(self.export_expression))
        #gevent.joinall(greenlets_list_third)
        self.logger.info("导表完成")


    def export_test(self):
        gevent.sleep()
        self.api_qc = self.api.ref_rna_qc
        from bson import ObjectId
        self.group_id = ObjectId("59ae0a75a4e1af55d523f91a")
        # self.api_qc.add_control_group(self.option("control").prop["path"], self.group_id)

    @time_count
    def export_qc(self):
        gevent.sleep()
        self.api_qc = self.api.api("denovo_rna_v2.denovoass_rna_qc")
        qc_stat = self.qc_stat_before.output_dir
        fq_type = self.option("fq_type").lower()
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before")
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat" #质控数据结果统计
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat" #原始数据结果统计
        self.api_qc.add_gragh_info(quality_stat_before, "before")
        qc_stat = self.qc_stat_after.output_dir
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="after")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        if self.option("group_table").is_set:
            self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(self.option("group_table").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
        if self.option("control").is_set:
            self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control").prop["path"], self.group_id)
            self.compare_detail = compare_detail

    @time_count
    def export_denovo_assembly(self):
        '''
        导入组装结果表格 liubinxu
        '''
        gevent.sleep()
        self.api_denovo_assembly = self.api.api("denovo_rna_v2.denovoass_assemble")
        assemble_filter_dir = self.assemble_filter.output_dir
        filter_evolution = self.filter_evaluation.output_dir

        if self.option("assemble_soft").lower() == "trinity" and self.option("assemble_method") == "total":
            unigene_evaluation = self.unigene_evaluation.output_dir
            filter_unigene_evaluation = self.filter_unigene_evaluation.output_dir
        else:
            unigene_evaluation = assemble_filter_dir
            filter_unigene_evaluation = filter_evolution
  
        result_dir = self.assemble.output_dir
        self.logger.info("{}".format(self.option("optimize")))
        if self.option("optimize") == "True" or  self.option("optimize") ==True:
            self.api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation, filter_evolution, filter_unigene_evaluation)
        else:
            self.api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation)

    @time_count
    def export_denovo_align(self):
        '''
        导入mapping率统计表格 liubinxu
        '''
        gevent.sleep()
        import pandas as pd
        self.api_denovo_align = self.api.api("denovo_rna_v2.denovoass_align")
        result_dir = self.align.output_dir
        if os.path.exists(result_dir + "/alignment_rate.txt"):
            align = pd.read_table(result_dir + "/alignment_rate.txt", sep= "\t", header=0, index_col=0)
            if self.option("fq_type") == "PE":
                align["aligned_reads"] = align["aligned_reads"] * 2
                align["total_reads"] = align["total_reads"] * 2
            align.rename(columns={"aligned_rate": "align_rate"}, inplace=True)
            align.to_csv(result_dir + "/alignment_rate.txt", sep="\t", columns=['total_reads', 'aligned_reads', 'align_rate'])
        else:
            count = pd.read_table(result_dir + '/transcript.count.matrix', sep="\t", header=0, index_col=0)
            align = pd.read_table(self.qc_stat_after.output_dir + "/fastq_stat.xls", sep= "\t", header=0, index_col=0)
            if self.option("fq_type") == "PE":
                align["aligned_reads"] = (count.sum() + 0.5).map(int) * 2
            else:
                align["aligned_reads"] = (count.sum() + 0.5).map(int)

            align['align_rate'] = align["aligned_reads"]/align["Total_Reads"]
            align = align.rename(str.lower, axis='columns')
            align.to_csv(result_dir + "/alignment_rate.txt", sep="\t", columns=['total_reads', 'aligned_reads', 'align_rate'])
        self.api_denovo_align.run(result_dir)

    @time_count
    def export_denovo_annotation(self):
        '''
        导入注释结果表格 liubinxu
        '''
        gevent.sleep()
        self.api_denovo_annotation = self.api.api("denovo_rna_v2.denovoass_annotation")
        result_dir = self.annotation.output_dir
        
        trans2gene = os.path.join(self.assemble_filter.output_dir,"Trinity.filter.gene_trans_map")
        # glob.glob(os.path.join(self.assemble.output_dir,"*.gene_trans_map"))[0]
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
        samples = list()
        for g,s in self.option("group_table").get_group_spname().items():
            samples.extend(s)
        self.api_denovo_annotation.run(result_dir, trans2gene, params,
                                       gene_exp = self.align.output_dir + '/gene.count.matrix',
                                       trans_exp = self.align.output_dir + '/transcript.count.matrix',
                                       gene_species_count = self.work_dir + '/gene_species_count.xls',
                                       trans_species_count = self.work_dir + '/transcript_species_count.xls',
                                       samples = samples,
                                       taxon=self.option("kegg_database"))



