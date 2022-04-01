# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""miRNA一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import json
import shutil
import re
import time
import gevent
import functools
from biocluster.config import Config
import pandas as pd
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
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


class SmallrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        workflow option参数设置
        """
        self._sheet = wsheet_object
        super(SmallrnaWorkflow, self).__init__(wsheet_object)
        options = [
            # ------ 基础参数设置 ------
            # --- 样本选择 ---
            {'name': 'sample_num', 'type': 'string', 'default': 'multiple'},  # 样本个数 ['multiple', 'single']
            # --- 常用参数 ---
            {'name': 'datatype', 'type': 'string', 'default': 'rawdata'},  # 数据类型 ['rawdata', 'cleandata']
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE OR other，为other时，跳过质控步骤
            {"name": "read_length", "type": "int", "default": "150"},  # 测序读长，50,75,150
            {"name": "lib_type", "type": "string", "default": "other"},  # 测序类型，typeI OR type II OR type III or other
            {"name": "cut_5", "type": "int", "default": "4"},
            # 当libtype为type I时，默认去除5'端4bp的碱基，后面仅需去除接头，,然后截取75bp的序列，最后再去除3'端4bp
            {"name": "extract_length", "type": "int", "default": "75"},
            # 当libtype为type II时，默认去除5'端3bp的碱基，然后去除接头，最后截取75bp的序列
            {"name": "cut_3", "type": "int", "default": "4"},  # 当libtype为type III时，直接截取75bp的序列，最后再去除3'端15bp
            {"name": "adapter", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # adapter序列
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {'name': 'qc_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 质控序列文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'}, # 上机名称
            {"name": "control_file", "type": "infile", "format": "sample.control_table"},  # 对照表
            {"name": "taxonmy", "type": "string", "default": ""},  # 物种类别，Animal OR Plant
            {"name": "ref_genome", "type": "string", "default": ""},  # 参考基因组，具体物种名称
            {'name': 'genome_id', 'type': 'string', 'default': None},  # 基因组编号 sg_genome_db.genome_id
            # --- sRNA鉴定 ---
            {"name": "input_type", "type": "string", "default": "clean"},  # srna模块分析的输入序列，raw为质控后，clean为mapping后
            {"name": "mirna_database", "type": "string", "default": "mirbase"},
            # 已知miRNA鉴定的参考数据库, mirbase or pmiren or none
            {"name": "mirna_specie", "type": "string", "default": ""},  # mirna 参考物种名，多选时分号分隔
            {"name": "organism_list", "type": "infile", 'format': "small_rna.common"},  # known mirna鉴定的物种列表
            # ------ 高级参数设置 ------
            # --- miRNA鉴定 ---
            {"name": "mismatch", "type": "int", "default": 1},  # 鉴定已知miRNA时允许的错配个数
            {"name": "mirna_soft", "type": "string", "default": "miRDeep2"},  # mirdeep2,mireap,mirdp2,mir_prefer
            {'name': 'rpm', 'type': 'float', 'default': 1},  # mirdeep-p2参数，用于测序数据过滤
            # --- miRNA碱基编辑分析 ---
            {"name": "is_miRNA_edit", "type": "bool", "default": True},  # 是否分析miRNA碱基编辑
            # --- 靶基因预测 ---
            {"name": "target_predict", "type": "bool", "default": False},  # 是否进行靶基因预测
            {"name": "analysis_strategy", "type": "string", "default": "ref"},  # 分析策略, ref, denovo
            {"name": "target_source", "type": "string", "default": "ref"},  # 靶基因序列来源, mongo, auto, ref
            {"name": "target_seq", "type": "infile", 'format': "sequence.fasta"},  # 预测靶基因所用序列
            {"name": "relate_task_id", "type": "string", "default": "ref"},  # 有参或无参转录组task_id
            {"name": "target_animal_predict_method", "type": "string", "default": ""},
            # 动物预测靶基因所用方法，多选时用分号分隔, miRanda|RNAhybrid|TargetScan
            {"name": "target_plant_predict_method", "type": "string", "default": ""},
            # 植物预测靶基因所用方法，多选时用分号分隔, psRobot|RNAhybrid|Targetfinder
            {"name": "min_support", "type": "int", "default": 1},  # 候选靶基因设定标准
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"},  # trinity组装结果文件
            {"name": "gene_to_trans", "type": "infile", 'format': "denovo_rna_v2.common"},  # 基因和转录本对应关系文件
            # --- 表达差异分析 ---
            {"name": "diff_method", "type": "string", "default": "DESeq2"},  # 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "diff_fdr_ci", "type": "float", "default": 0.05},  # 显著性水平
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  # 选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  # Bonferroni,Holm,BH,BY
            # --- 靶基因转录组功能注释 ---
            {"name": "nr_database", "type": "string", "default": ""},  # nr库一级分类
            {"name": "nr_sub_database", "type": "string", "default": ""},  # nr库二级分类
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
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},  # string比对使用的e值
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg,swissprot,pfam'},
            # --- 分析流程版本 ---
            {"name": "version", "type": "string", "default": "v2"},
            {'name': 'get_run_log', 'type': 'bool', 'default': False},
        ]
        # 获取输出目录
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False


        self.logger.info("options is  {}".format(self._sheet.options()))
        if self.option("version") == "v2":
            if self.option("target_seq").is_set:
                self.option("assembly_file", self.option("target_seq").path)
                self.logger.info("序列文件 {}".format(self.option("assembly_file").path))

        # 添加tool/module
        self.filecheck = self.add_tool("small_rna.file_check")
        self.download = self.add_tool("small_rna_v2.transfer_annot")

        self.gunzip = self.add_tool("small_rna.gzfastq2fastq")
        self.qc = self.add_module("small_rna.mirna_qc")
        self.uniq = self.add_tool("small_rna.fasta_uniq")
        self.qc_stat_before = self.add_module("small_rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("small_rna.hiseq_reads_stat")
        self.mapping = self.add_tool("small_rna.mapper_and_stat")
        self.annotation = self.add_module("denovo_rna_v2.denovo_annotation")
        self.cds_predict = self.add_module("denovo_rna_v2.cds_tf")
        self.srna = self.add_module("small_rna_v2.srna.srna")
        self.mirna_target = self.add_module("small_rna.target_predict")
        if self.option("sample_num") == "multiple":
            self.diffexpress = self.add_tool("small_rna_v2.diff.diffexp_batch")
            if self.option("group_table").prop["sample_number"] > 2:
                self.exp_pca = self.add_tool("small_rna.exp_pca")
                self.exp_corr = self.add_tool("small_rna.exp_corr")
            group_spname = self.option("group_table").get_group_spname()
            group_snum = [len(group_spname[g]) for g in group_spname]
            min_group_num = min(group_snum)
            if min_group_num >= 3:
                # self.ellipse = self.add_tool('graph.ellipse')
                self.ellipse = self.add_tool('graph.ellipse')

        if self.option("is_miRNA_edit"):
            self.mirna_edit = self.add_module("small_rna.mirna_edit")
        self.mirna_family = self.add_tool("small_rna.smallrna_family_analyse")
        self.mirna_atcg_bias = self.add_tool("small_rna.atcg_bias")
        self.gene_fa = self.add_tool("ref_rna_v2.gene_fa")
        self.transcripts = self.add_tool("ref_rna_v2.transcript_abstract")
        self.extract_utr3 = self.add_tool("small_rna.extract_utr3")
        self.sequence_detail = self.add_tool("small_rna_v2.sequence_detail")
        # 判断流程结束tool/module list
        self.final_tools = list()
        self.final_tools.append(self.qc_stat_after)
        self.final_tools.append(self.qc_stat_before)
        if self.option("sample_num") == "multiple":
            self.final_tools.append(self.diffexpress)
            if self.option("group_table").prop["sample_number"] > 2:
                self.final_tools.append(self.exp_pca)
                self.final_tools.append(self.exp_corr)
            if hasattr(self, 'ellipse'):
                self.final_tools.append(self.ellipse)
        if self.option("target_predict"):
            if self.option("analysis_strategy") != "denovo":
                self.final_tools.append(self.gene_fa)
            self.final_tools.append(self.mirna_target)
            if not self.option("assembly_file").is_set and self.option("target_source") == "ref":
                self.final_tools.append(self.transcripts)
        if self.option("assembly_file").is_set:
            self.final_tools.append(self.annotation)
        if self.option("mirna_database").lower() != "none":
            if self.option("is_miRNA_edit"):
                self.final_tools.append(self.mirna_edit)
            self.final_tools.append(self.mirna_atcg_bias)
            self.final_tools.append(self.mirna_family)
        self.logger.info(self.final_tools)

        # 添加step，显示在页面进度条
        self.step.add_steps("filecheck", "rna_qc", "gunzip", "uniq", "mapping", "srna",
                            "gene_fa", "annotation", "cds_predict", "exp_corr")
        if self.option("sample_num") == "multiple":
            self.step.add_steps("diffexpress")
            if self.option("group_table").prop["sample_number"] > 2:
                self.step.add_steps("exp_pca")
            if hasattr(self, 'ellipse'):
                self.step.add_steps("ellipse")
        if self.option("target_predict"):
            self.step.add_steps("mirna_target")
        if self.option("is_miRNA_edit"):
            self.step.add_steps("mirna_edit")
        if self.option("mirna_database").lower() != "none":
            self.step.add_steps("mirna_atcg_bias", "mirna_family")
        if self.option("target_source") == "mongo":
            self.step.add_steps("download")

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

        # data = os.path.join(self.work_dir, 'data.json')
        # if os.path.exists(data):
        #     with open(data, 'r') as load_f:
        #         load_dict = json.load(load_f)
        #         if 'rerun' in load_dict and load_dict['rerun']:
        #             self.logger.info("该项目重运行中，先删除mongo库中已有数据")
        #             self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'small_rna')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))
        delete = DeleteDemoMongo(self.task_id, 'small_rna')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'ref_rna_v2')
        # code = os.system(cmd)
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        ## 检查物种参数
        if self.option("mirna_database").lower() != "none":
            if not self.option("mirna_specie") or self.option("mirna_specie") == "null":
                raise OptionError("必须指定具体物种")
        elif self.option("mirna_specie").lower() not in ["all", "auto"]:
            species = set(self.option("mirna_specie").split(","))
            species_list = list()
            for i in species:
                i = i.strip()
                if len(i) >= 3:
                    if self.option("database").lower() == "pmiren":
                        pmiren = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')][
                            'pmiren']
                        try:
                            organism = pmiren.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                    else:
                        mirbase = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')][
                            'mirbase']
                        try:
                            organism = mirbase.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                species_list.append(i)
            self.option("mirna_specie", ",".join(species_list))
        # 检查靶基因预测参数
        if self.option("target_predict"):
            if not (self.option("target_animal_predict_method") or self.option("target_plant_predict_method")):
                raise OptionError('靶基因预测：软件参数必须传递且不能为空')
            if self.option("taxonmy").lower() == "animal":
                if len(self.option("target_animal_predict_method").split(",")) < self.option("min_support"):
                    raise OptionError('靶基因预测：选择的软件数量应该>=支持的软件数量')
            elif self.option("taxonmy").lower() == "plant":
                if len(self.option("target_plant_predict_method").split(",")) < self.option("min_support"):
                    raise OptionError('靶基因预测：选择的软件数量应该>=支持的软件数量')
        # 当上传的数据为cleandata时，重新定义参数，使之与原先工作流兼容
        if self.option("datatype") == "cleandata":
            self.option("fq_type", "other")
            self.option("fastq_dir", self.option("qc_dir").path)
        # miRNA edit
        if self.option("mirna_database").lower() == "none":
            self.option("is_miRNA_edit", False)
        # 判断样本个数选择、分组方案以及is_duplicate是否正确
        if self.option('group_table').is_set:
            dup = False
            group_dict = self.option("group_table").prop['group_dict']
            for g in group_dict.keys():
                if len(group_dict[g]) > 1:
                    dup = True
            if self.option('is_duplicate') != dup:
                self.logger.info(self.option('is_duplicate'))
                self.logger.info(dup)
                raise OptionError('生物学重复参数选择和分组方案不匹配，请检查输入是否有误')
            if self.option('group_table').prop['sample_number'] >= 2:
                if self.option('sample_num') == 'single':
                    raise OptionError('分组方案表中有{}个样本，与样本个数为{}相互矛盾，请重新确认参数.'.format(
                        self.option('group_table').prop['sample_number'], self.option('sample_num')))
        else:
            samples = list()
            for line in open(os.path.join(self.option('fastq_dir').path, 'list.txt')):
                sample = line.strip().split('\t')[1]
                if sample not in samples:
                    samples.append(sample)
            group_table = os.path.join(self.work_dir, 'group.txt')
            with open(group_table, 'w') as w:
                w.write('#sample\tgroup\n')
                for sample in samples:
                    w.write('{}\t{}\n'.format(sample, sample))
            self.option("group_table", group_table)
        if self.option('sample_num') == 'single':
            self.option('is_duplicate', False)
        if self.option("taxonmy").lower() == "animal" and self.option("mirna_database").lower() == "pmiren":
            raise OptionError("PmiREN数据库目前只支持植物物种")
        if not self.option("fq_type") in ["PE", "SE", "other"]:
            raise OptionError("fq序列类型应为PE或SE或other")
        if not re.match(r'^[ATCGN]*$', self.option("adapter")):
            raise OptionError("接头序列输入错误")
        if not self.option("lib_type") in ["type I", "type II", "type III", "other"]:
            raise OptionError("建库类型必须为type I或type II或type III或other")
        if self.option("lib_type") == "type I":
            if not (int(self.option("cut_5")) == 4 and int(self.option("cut_3")) == 4 and int(
                    self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type I时，5'端和3'端的截取序列长度必须为4，序列提取长度必须为75")
        elif self.option("lib_type") == "type II":
            if not (int(self.option("cut_5")) == 3 and int(self.option("cut_3")) == 0 and int(
                    self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type II时，5'端截取序列的长度必须为3, 3'端截取序列的长度必须为0，序列提取长度必须为75")
        elif self.option("lib_type") == "type III":
            if not (int(self.option("cut_5")) == 3 and int(self.option("cut_3")) == 15 and int(
                    self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type III时，5'端截取序列的长度必须为0, 3'端截取序列的长度必须为15，序列提取长度必须为75")
        if not self.option("quality_score_system").lower() in ["phred+33", "phred 33", "phred+64", "phred 64"]:
            raise OptionError("测序质量参数输入有误")
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        col = db["sg_genome_db"]
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        if self.option("genome_id"):
            genome_info = col.find_one({'genome_id': self.option('genome_id')})
        else:
            self.set_error('缺少参数 genome_id')
        if not genome_info:
            self.set_error('数据库中不存在该物种的注释信息，程序退出')
        try:
            self.ref_annot_dir = os.path.join(db_path, genome_info["anno_path_v2"])
            self.annot_class = os.path.join(db_path, genome_info["anno_path_v2"], "annot_class")
            self.g2t2p = os.path.join(db_path, genome_info["anno_path_v2"], "annot_class", "tran2gene.txt")
            self.ref_genome = os.path.join(db_path, genome_info["dna_fa"])
            self.ref_index = os.path.join(db_path, genome_info["dna_index"])
            self.ref_gtf = os.path.join(db_path, genome_info["gtf"])
            self.genome_version = genome_info["accession"]
            self.genome_id = genome_info["genome_id"]
            self.des = os.path.join(db_path, genome_info["bio_mart_annot"])
            self.des_type = genome_info["biomart_gene_annotype"]
            self.known_cds = os.path.join(db_path, genome_info["cds"])
            self.known_pep = os.path.join(db_path, genome_info["pep"])
            self.trans_fa = os.path.join(db_path, genome_info["transcript"])
            self.entrez = os.path.join(db_path, genome_info["ensemble2entrez"])
            self.genome_stat = os.path.join(db_path, genome_info["gene_stat"])
            self.species_name = genome_info["name"]
            self.species = genome_info["taxon_id"]
            self.ref_anno_version = genome_info["assembly"]
            self.hyperlink = genome_info["ensemble_web"]
            self.anno_detail = os.path.join(db_path, genome_info["anno_path_v2"],
                                            "annot_class/anno_stat/all_anno_detail.xls")
        except:
            self.set_error("数据库中该物种信息不全，程序退出")
        if os.path.exists(os.path.join(self.ref_annot_dir, "repeatmasker")):
            self.repeatmasker = True
        else:
            self.repeatmasker = False
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
        elif self.option("taxonmy") == "Animal":
            self.option("nr_database", "metazoa")
        elif self.option("taxonmy") == "Plant":
            self.option("nr_database", "viridiplantae")
        else:
            # nr = self.option("nr_database").lower()
            self.option("nr_database", "nr")
        if self.option("target_source") == "mongo":
            if not self.option("relate_task_id"):
                raise OptionError("选择转录组数据需要填写task_id")
        elif self.option("target_source") == "auto":
            if not self.option("target_seq").is_set:
                raise OptionError("自定义上传需要填写靶基因序列文件")
            if not self.option("gene_to_trans").is_set:
                raise OptionError("自定义上传需要填写基因与转录本对应关系文件")

        if self.option('sample_num') == 'multiple':
            group_size = list()
            group_dict = self.option("group_table").prop['group_dict']
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
                    self.logger.info("该项目有生物学重复,不可以使用DEGseq,DESeq2,阈值设置为0.05")
            else:
                pass
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
        mirna workflow run方法
        """
        self.filecheck.on('end', self.run_gunzip)
        self.gunzip.on('end', self.run_qc)
        self.qc.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc.on('end', self.run_uniq)
        self.uniq.on('end', self.run_mapping)
        self.mapping.on('end', self.run_srna)
        if self.option("mirna_database").lower() != "none":
            if self.option("is_miRNA_edit"):
                self.srna.on('end', self.run_mirna_edit)
            self.srna.on('end', self.run_mirna_atcg_bias)
            self.srna.on('end', self.run_mirna_family)
        if self.option("sample_num") == "multiple":
            self.srna.on('end', self.run_diffexpress)
            if self.option("group_table").prop["sample_number"] > 2:
                self.srna.on('end', self.run_exp_pca)
                self.srna.on('end', self.run_exp_corr)

        # 靶基因预测逻辑
        if self.option("target_predict"):
            if self.option("target_source") != "mongo":
                if self.option("analysis_strategy") != "denovo":
                    self.filecheck.on('end', self.run_gene_fa)
            else:
                if self.option("analysis_strategy") != "denovo":
                    self.download.on('end', self.run_gene_fa)
            if self.option("target_source") == "mongo":
                self.filecheck.on('end', self.run_download)
            if self.option("target_source") == "ref":
                self.filecheck.on('end', self.run_transcripts)
            if self.option("assembly_file").is_set and self.option("target_source") == "auto":
                self.filecheck.on('end', self.run_diamond)
                self.filecheck.on('end', self.run_cds_predict)
            if self.option("taxonmy").lower() == "animal":
                if self.option("target_source") == "mongo":
                    self.download.on('end', self.run_extract_utr3)
                    self.on_rely([self.download, self.extract_utr3, self.srna], self.run_mirna_target)
                elif self.option("target_source") == "auto":
                    self.filecheck.on('end', self.run_extract_utr3)
                    self.on_rely([self.annotation, self.extract_utr3, self.srna], self.run_mirna_target)
                else:
                    self.filecheck.on('end', self.run_extract_utr3)
                    self.on_rely([self.extract_utr3, self.srna], self.run_mirna_target)
            else:
                if self.option("target_source") == "mongo":
                    self.on_rely([self.download, self.srna], self.run_mirna_target)
                elif self.option("target_source") == "auto":
                    self.on_rely([self.annotation, self.srna], self.run_mirna_target)
                else:
                    self.srna.on('end', self.run_mirna_target)
        self.on_rely(self.final_tools, self.run_seq_detail)
        self.run_filecheck()
        super(SmallrnaWorkflow, self).run()

    def run_transcripts(self):
        opts = {
            "ref_genome": self.ref_genome,
            "ref_genome_gtf": self.ref_gtf
        }
        self.transcripts.set_options(opts)
        self.transcripts.on("end", self.set_output, "Transcripts")
        self.transcripts.run()

    def run_download(self):

        if self.option("analysis_strategy") == "ref":
            task_type = "ref_rna_v2"
        elif self.option("analysis_strategy") == "denovo":
            task_type = "denovo_rna_v2"
        else:
            self.set_error("不支持该类型的项目")
        opts = {
            'task_id': self.option("relate_task_id"),
            'task_type': task_type
        }
        self.download.set_options(opts)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.run()

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'sample_num': self.option('sample_num'),
            'is_duplicate': self.option('is_duplicate')
        }
        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        if self.option("control_file").is_set:
            opts.update({'control_file': self.option('control_file')})
        if self.option("assembly_file").is_set:
            opts.update({'assembly_file': self.option('assembly_file')})
        if self.option("gene_to_trans").is_set:
            opts.update({'gene_to_trans': self.option('gene_to_trans').prop["path"]})
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_diamond(self):
        self.logger.info("开始运行blast注释")
        self.blast_modules = list()
        blast_opts = {
            'query': self.option("assembly_file"),
            'query_type': 'nucl',
            'database': None,
            'blast': 'blastx',
            'evalue': None,
            'outfmt': 5,
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.diamond_nr = self.add_module("denovo_rna_v2.diamond")
            blast_opts.update(
                {
                    'database': self.option("nr_database"),
                    'evalue': self.option('nr_blast_evalue')
                }
            )
            self.diamond_nr.set_options(blast_opts)
            self.blast_modules.append(self.diamond_nr)
            self.diamond_nr.on('end', self.set_output, 'diamond_nr')
            self.diamond_nr.run()
        if 'cog' in self.option('database'):
            self.diamond_string = self.add_module("denovo_rna_v2.diamond")
            blast_opts.update(
                {'database': 'string', 'evalue': self.option('string_blast_evalue')}
            )
            self.diamond_string.set_options(blast_opts)
            self.blast_modules.append(self.diamond_string)
            self.diamond_string.on('end', self.set_output, 'diamond_string')
            self.diamond_string.run()
        if 'kegg' in self.option('database'):
            self.diamond_kegg = self.add_module("denovo_rna_v2.diamond")
            blast_opts.update(
                {'database': 'kegg', 'evalue': self.option('kegg_blast_evalue')}
            )
            self.diamond_kegg.set_options(blast_opts)
            self.blast_modules.append(self.diamond_kegg)
            self.diamond_kegg.on('end', self.set_output, 'diamond_kegg')
            self.diamond_kegg.run()
        if 'swissprot' in self.option('database'):
            self.blast_swissprot = self.add_module('denovo_rna_v2.blast')
            blast_opts.update(
                {'database': 'swissprot', 'evalue': self.option('swissprot_blast_evalue')}
            )
            self.blast_swissprot.set_options(blast_opts)
            self.blast_modules.append(self.blast_swissprot)
            self.blast_swissprot.on('end', self.set_output, 'blast_swissprot')
            self.blast_swissprot.run()
        self.on_rely([self.diamond_nr, self.diamond_kegg, self.diamond_string, self.blast_swissprot, self.cds_predict],
                     self.run_annotation)

    def run_cds_predict(self):
        self.logger.info("开始运行cds预测")
        opts = {
            "fasta": self.option("assembly_file"),
            "e_value": self.option('pfam_blast_evalue'),
            "species_type": self.option("taxonmy"),
            "isoform_unigene": os.path.join(self.filecheck.work_dir, "Trinity_t2g2u"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.set_output, "cds_predict")
        self.cds_predict.on('start', self.set_step, {'start': self.step.cds_predict})
        self.cds_predict.on('end', self.set_step, {'end': self.step.cds_predict})
        self.cds_predict.run()

    def run_annotation(self):
        self.logger.info("开始运行注释统计")
        anno_opts = {
            "gene2trans": os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map"),
            "go_annot": True,
            "nr_annot": True,
            "taxonomy": self.option("kegg_database"),
            "blast_nr_xml": self.diamond_nr.option('outxml'),
            "blast_kegg_xml": self.diamond_kegg.option('outxml'),
            "blast_string_xml": self.diamond_string.option('outxml'),
            "blast_swissprot_xml": self.blast_swissprot.option('outxml'),
            "pfam_domain": self.cds_predict.output_dir + "/pfam_domain"
        }
        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        self.annotation.run()

    def run_gunzip(self):
        opts = {
            'fastq_path': self.option('fastq_dir').prop["path"]
        }
        self.gunzip.set_options(opts)
        self.gunzip.on('start', self.set_step, {'start': self.step.gunzip})
        self.gunzip.on('end', self.set_step, {'end': self.step.gunzip})
        self.gunzip.run()

    def run_qc(self):
        list_path = self.option("fastq_dir").prop["path"] + "/list_re.txt"
        if self.option('quality_score_system').lower() == "phred+33" or self.option(
                'quality_score_system').lower() == "phred 33":
            fastq_format = "Q33"
        elif self.option('quality_score_system').lower() == "phred+64" or self.option(
                'quality_score_system').lower() == "phred 64":
            fastq_format = "Q64"
        if self.option("fq_type") == "other":
            skip_qc = "yes"
        else:
            skip_qc = "no"
        self.qc.set_options({
            'list_file': list_path,
            'adapter': self.option("adapter"),
            'fastq_format': fastq_format,
            'cut_left': int(self.option('cut_5')),
            'cut_tail': int(self.option('cut_3')),
            'extract_length': int(self.option("extract_length")),
            'skip_qc': skip_qc
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.rna_qc})
        self.qc.on('end', self.set_step, {'end': self.step.rna_qc})
        self.qc.run()

    def run_uniq(self):
        opts = {
            "config": self.qc.option("config_file"),
        }
        self.uniq.set_options(opts)
        self.uniq.on("end", self.set_output, "uniq")
        self.uniq.on("start", self.set_step, {"start": self.step.uniq})
        self.uniq.on("end", self.set_step, {"end": self.step.uniq})
        self.uniq.run()

    def run_qc_stat(self, event):
        if self.option('quality_score_system').lower() == "phred+33" or self.option(
                'quality_score_system').lower() == "phred 33":
            quality = 33
        else:
            quality = 64
        if event['data']:
            self.qc_stat_after.set_options({
                'fastq_dir': self.qc.option('cleandata'),
                'fq_type': "SE",
                'quality': quality,
                'dup': True
            })
        else:
            self.qc_stat_before.set_options({
                'fastq_dir': self.qc.option('rawdata'),
                'fq_type': "SE",
                'quality': quality,
                'dup': False
            })
        if event['data']:
            self.qc_stat_after.on('end', self.set_output, 'qc_stat_after')
            self.qc_stat_after.run()
        else:
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
            self.qc_stat_before.run()

    def run_mapping(self):
        opts = {
            "config": self.qc.option("config_file"),
            "ref": self.ref_genome,
            "fasta": os.path.join(self.uniq.output_dir, "uniq.fasta"),
            "index": self.ref_index,
            "gtf": self.ref_gtf
        }
        self.mapping.set_options(opts)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_srna(self):
        self.logger.info("开始运行srna注释")
        opts = {
            "category": self.option("taxonmy"),
            "reference": self.ref_genome,
            "config": self.qc.option("config_file"),
            "arf": os.path.join(self.mapping.output_dir, "reads_vs_genome.arf"),
            "repeat": self.repeatmasker,
            "gtf": self.ref_gtf,
            "list": self.option("fastq_dir").prop["path"] + "/list_re.txt",
            "qc_output": self.qc.output_dir,
            "mismatch": self.option("mismatch"),
            "input_type": self.option("input_type"),
            "database": self.option("mirna_database").lower(),
            "index": self.ref_index,
            "method": self.option("mirna_soft"),
            "rpm": self.option("rpm"),
        }
        if self.option("mirna_database").lower() != "none":
            opts.update({'species': self.option("mirna_specie")})
        if self.option("organism_list").is_set:
            opts.update({'organism_list': self.option("organism_list")})
        if self.option("input_type") == "raw":
            opts.update({'clean_fa': os.path.join(self.uniq.output_dir, "uniq.fasta")})
        else:
            opts.update({'clean_fa': os.path.join(self.mapping.output_dir, "uniq_mapped.fasta")})
        if self.repeatmasker == True:
            repeat_fa = glob.glob(self.ref_annot_dir + "/repeatmasker/*.fa")[0]
            repeat_gff = glob.glob(self.ref_annot_dir + "/repeatmasker/*.gff")[0]
            opts.update({'repeat_fa': repeat_fa})
            opts.update({'repeat_gff': repeat_gff})
        self.srna.set_options(opts)
        self.srna.on("end", self.set_output, "srna")
        self.srna.on("start", self.set_step, {"start": self.step.srna})
        self.srna.on("end", self.set_step, {"end": self.step.srna})
        self.srna.run()

    def run_mirna_edit(self):
        specie = self.option("mirna_specie").split(",")[0]
        mature_fa = self.srna.output_dir + "/known_mirna/mature.fa"
        if self.option("mirna_database").lower() == "mirbase":
            hairpin_fa = self.config.SOFTWARE_DIR + "/database/mirbase/hairpin.fa"
        elif self.option("mirna_database").lower() == "pmiren":
            hairpin_fa = self.config.SOFTWARE_DIR + "/database/PmiREN/hairpin.fa"
        opts = {
            "list_file": os.path.join(self.qc.output_dir, "clean_data/list.txt"),
            "species": specie,
            "hairpin_fa": hairpin_fa,
            "mature_fa": mature_fa,
            "index": self.ref_index
        }
        self.mirna_edit.set_options(opts)
        self.mirna_edit.on("end", self.set_output, "mirna_edit")
        self.mirna_edit.on("start", self.set_step, {"start": self.step.mirna_edit})
        self.mirna_edit.on("end", self.set_step, {"end": self.step.mirna_edit})
        self.mirna_edit.run()

    def run_mirna_family(self):
        opts = {
            "mir": self.srna.output_dir + "/known_mirna_count.xls",
            "matfa": self.srna.output_dir + "/known_mirna/mature.fa",
            "novofa": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
        }
        if self.option("mirna_database").lower() == "mirbase":
            opts.update({'family': self.config.SOFTWARE_DIR + '/database/mirbase/family_analyse/miRfamily.dat'})
            opts.update({'pre': self.config.SOFTWARE_DIR + '/database/mirbase/family_analyse/miR_pre_mature.dat'})
        elif self.option("mirna_database").lower() == "pmiren":
            opts.update({'family': self.config.SOFTWARE_DIR + '/database/PmiREN/family_analyse/miRfamily.dat'})
            opts.update({'pre': self.config.SOFTWARE_DIR + '/database/PmiREN/family_analyse/miR_pre_mature.dat'})
        self.mirna_family.set_options(opts)
        self.mirna_family.on("end", self.set_output, "mirna_family")
        self.mirna_family.on("start", self.set_step, {"start": self.step.mirna_family})
        self.mirna_family.on("end", self.set_step, {"end": self.step.mirna_family})
        self.mirna_family.run()

    def run_mirna_atcg_bias(self):
        opts = {
            "known": self.srna.output_dir + "/known_mirna/mature.fa",
            "novel": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
        }
        self.mirna_atcg_bias.set_options(opts)
        self.mirna_atcg_bias.on("end", self.set_output, "mirna_atcg_bias")
        self.mirna_atcg_bias.on("start", self.set_step, {"start": self.step.mirna_atcg_bias})
        self.mirna_atcg_bias.on("end", self.set_step, {"end": self.step.mirna_atcg_bias})
        self.mirna_atcg_bias.run()

    def run_extract_utr3(self):
        if self.option("analysis_strategy") == "ref":
            if self.option("target_source") == "ref":
                opts = {
                    "ref": self.ref_genome,
                    "gtf": self.ref_gtf
                }
            elif self.option("target_source") == "mongo":
                opts = {
                    "ref": self.download.option("genom/mnt/ilustre/users/isanger/wpm2/workspace/20211124/Smallrna_ljrb_r243gdacdbf0run02iv7g5/DenovoAnnotation/output/annot_stat/trans_anno_detail.xls").prop["path"],
                    "gtf": self.download.option("assembly_gtf").prop["path"]
                }
            else:
                opts = {
                    "transcript": self.option("assembly_fa")
                }
        elif self.option("analysis_strategy") == "denovo":
            if self.option("target_source") == "mongo":
                opts = {
                    "transcript": self.download.option("assembly_fa").prop["path"]
                }
            else:
                opts = {
                    "transcript": self.option("assembly_file").prop['path']
                }

        self.extract_utr3.set_options(opts)
        self.extract_utr3.run()

    def run_mirna_target(self):
        if self.option("target_seq").is_set:
            self.anno_detail = os.path.join(self.annotation.output_dir, "anno_stat/trans_anno_detail.xls")
            opts = {
                "known": self.srna.output_dir + "/known_mirna/mature.fa",
                "novol": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
                "ref": self.option("target_seq").prop['path'],
                "type": self.option("taxonmy").lower(),
                "species": self.option("ref_genome"),
                "version": "v1.1",
                "anno_detail": self.anno_detail,
                "min_support": int(self.option("min_support"))
            }
            if self.option("taxonmy").lower() == "plant":
                opts.update({'method': self.option("target_plant_predict_method")})
            else:
                opts.update({'method': self.option("target_animal_predict_method")})
        else:
            if self.option("target_source") == "mongo":
                self.anno_detail = self.download.option("anno_detail").prop["path"]
            elif self.option("target_source") == "auto":
                self.anno_detail = os.path.join(self.annotation.output_dir, "anno_stat/trans_anno_detail.xls")
            else:
                pass

            if self.option("taxonmy").lower() == "plant":
                opts = {
                    "known": self.srna.output_dir + "/known_mirna/mature.fa",
                    "novol": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
                    "ref": self.trans_fa,
                    "method": self.option("target_plant_predict_method"),
                    "type": self.option("taxonmy").lower(),
                    "species": self.option("ref_genome"),
                    "version": "v1.1",
                    "anno_detail": self.anno_detail,
                    "min_support": int(self.option("min_support"))
                }
            else:
                opts = {
                    "known": self.srna.output_dir + "/known_mirna/mature.fa",
                    "novol": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
                    "ref": os.path.join(self.extract_utr3.output_dir, "utr3.fa"),
                    "method": self.option("target_animal_predict_method"),
                    "version": "v1.1",
                    "type": self.option("taxonmy").lower(),
                    "species": self.option("ref_genome"),
                    "anno_detail": self.anno_detail,
                    "min_support": int(self.option("min_support"))
                }
        if self.option("mirna_database").lower() == "none":
            os.system('touch {}'.format(self.work_dir + '/mature.fa'))
            opts.update({'known': self.work_dir + '/mature.fa'})
        self.mirna_target.set_options(opts)
        self.mirna_target.on('end', self.set_output, 'mirna_target')
        self.mirna_target.on('start', self.set_step, {'start': self.step.mirna_target})
        self.mirna_target.on('end', self.set_step, {'end': self.step.mirna_target})
        self.mirna_target.run()

    def run_gene_fa(self):
        if self.option("target_source") == "mongo":
            opts = {
                "ref_new_gtf": os.path.join(self.download.work_dir, "all_transcripts.gtf"),
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

    def run_exp_pca(self):
        self.logger.info("开始运行pca")
        opts = {
            "express_matrix": self.srna.option("total_mirna_norm").prop["path"]
        }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        if hasattr(self, 'ellipse'):
            self.exp_pca.on('end', self.run_ellipse)
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
        self.logger.info("开始运行聚类分析")
        opts = {
            "express_matrix": self.srna.option("total_mirna_norm").prop["path"]
        }
        self.exp_corr.set_options(opts)
        self.exp_corr.on("end", self.set_output, "exp_corr")
        self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.run()

    def run_diffexpress(self):
        self.logger.info("开始运行基因差异表达分析")
        count_file = os.path.join(self.srna.output_dir, "total_mirna_count.xls")
        tpm_file = os.path.join(self.srna.output_dir, "total_mirna_norm.xls")
        opts = {
            "count": count_file,
            "exp": tpm_file,
            "group": self.option("group_table").prop["path"],
            "cmp": self.option("control_file").prop["path"],
            "pvalue_padjust": self.option("pvalue_padjust").lower(),
            "pvalue": float(self.option("diff_fdr_ci")),
            "fc": float(self.option("fc")),
            "padjust_way": self.option("padjust_way"),
            "method": self.option("diff_method"),
        }
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

    def run_seq_detail(self):
        if self.option("target_source") == "mongo":
            transcript_fasta = os.path.join(self.download.work_dir, "all_transcripts.fa")
            gene_fasta = transcript_fasta
            t2g = os.path.join(self.download.work_dir, "trans2gene.txt")
        elif self.option("target_source") == "auto":
            transcript_fasta = self.option("assembly_file").prop["path"]
            gene_fasta = transcript_fasta
            t2g = self.option("gene_to_trans").prop["path"]
            cds_fasta = self.cds_predict.option("cds").prop["path"]
            pep_fasta = self.cds_predict.option("pep").prop["path"]
            biomart_file = self.annotation.output_dir + "/anno_stat/trans_anno_detail.xls"
            biomart_type = "denovo_annotation"
        else:
            transcript_fasta = self.transcripts.option("trans_fa").prop["path"]
            gene_fasta = self.gene_fa.option("gene_fa").prop["path"]
            t2g = self.transcripts.option("trans2gene").prop["path"]
            cds_fasta = self.known_cds
            pep_fasta = self.known_pep
            biomart_file = self.des
            biomart_type = self.des_type
        opts = {
            'cds_seq': cds_fasta,
            'pep_seq': pep_fasta,
            'txpt_seq': transcript_fasta,
            'gene_seq': gene_fasta,
            'trans2gene': t2g,
            'biomart_file': biomart_file,
            'biomart_type': biomart_type
        }
        self.sequence_detail.set_options(opts)
        self.sequence_detail.on('end', self.end)
        self.sequence_detail.run()


    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="13700319")
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
            self.move2outputdir(obj.output_dir, 'QC')
        if event['data'] == 'qc_stat_before':
            self.move2outputdir(obj.output_dir, 'QC/stat_rawdata')
        if event['data'] == 'qc_stat_after':
            self.move2outputdir(obj.output_dir, 'QC/stat_cleandata')
        if event['data'] == 'uniq':
            self.move2outputdir(obj.output_dir, 'Uniq')
        if event['data'] == 'mapping':
            self.move2outputdir(obj.output_dir, 'Mapping')
        if event['data'] == 'diffexpress':
            self.move2outputdir(obj.output_dir, 'DiffExpress')
        if event['data'] == 'exp_corr':
            self.move2outputdir(obj.output_dir, 'ExpCorr')
        if event['data'] == 'exp_pca':
            self.move2outputdir(obj.output_dir, 'ExpPca')
        if event['data'] == 'gene_fa':
            self.move2outputdir(obj.output_dir, 'GeneFa')
        if event['data'] == 'mirna_target':
            self.move2outputdir(obj.output_dir, 'MirnaTarget')
        if event['data'] == 'mirna_atcg_bias':
            self.move2outputdir(obj.output_dir, 'MirnaBias')
        if event['data'] == 'mirna_family':
            self.move2outputdir(obj.output_dir, 'MirnaFamily')
        if event['data'] == 'mirna_edit':
            self.move2outputdir(obj.output_dir, 'MirnaEdit')
        if event['data'] == 'srna':
            self.move2outputdir(obj.output_dir, 'Srna')
        if event['data'] == 'mapping':
            self.move2outputdir(obj.output_dir, 'Mapping')
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'Annotation')

    def end(self):
        self.run_api()
        ## 更新一系列主表的字段，用于页面交互分析
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        col = db["sg_task"]
        col.update({"task_id": self.task_id}, {"$set": {"ref_gtf": self.ref_gtf}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"ref_genome": self.ref_genome}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"genome_id": self.genome_id}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"organism_name": self.option("ref_genome")}}, upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"assembly": self.genome_version}}, upsert=True)
        if not self.option("assembly_file").is_set and self.option("target_predict"):
            col.update({'task_id': self.task_id},
                        {'$set': {'seqdetail': os.path.join(self.workflow_output, "09Other/SequenceDetail/")}},
                        upsert=True)
        if self.option("mirna_database").lower() != "none":
            col.update({"task_id": self.task_id},
                       {"$set": {
                           "known_pre": os.path.join(self.workflow_output, "04sRNA_Analysis/01Known_miRNA/known_miRNA_detail.xls")}},
                       upsert=True)
            col.update({"task_id": self.task_id},
                       {"$set": {
                           "novol_pre": os.path.join(self.workflow_output, "04sRNA_Analysis/02Novel_miRNA/novel_miRNA_detail.xls")}},
                       upsert=True)
            col.update({"task_id": self.task_id},
                       {"$set": {"novol": os.path.join(self.workflow_output, "04sRNA_Analysis/02Novel_miRNA/novel_mature.fa")}},
                       upsert=True)
            col.update({"task_id": self.task_id},
                       {"$set": {"known": os.path.join(self.workflow_output, "04sRNA_Analysis/01Known_miRNA/known_mature.fa")}},
                       upsert=True)
        else:
            col.update({"task_id": self.task_id},
                       {"$set": {
                           "novol_pre": os.path.join(self.workflow_output, "04miRNA_Analysis/miRNA_detail.xls")}},
                       upsert=True)
            col.update({"task_id": self.task_id},
                       {"$set": {"novol": os.path.join(self.workflow_output, "04miRNA_Analysis/miRNA_mature.fa")}},
                       upsert=True)
        col.update({"task_id": self.task_id}, {"$set": {"version": self.option("version")}}, upsert=True)
        if self.genome_version == "v2":
            col.update({'task_id': self.task_id},
                       {'$set': {'database_version': {"kegg": "202003"}}}, upsert=True)
        db2 = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname(mtype='ref_rna_v2', dydb_forbid=True)]
        col_genome_db = db2["sg_genome_db"]
        genome_doc = col_genome_db.find_one({'genome_id': self.option('genome_id')})

        if "database_version" in genome_doc:
            col.update({'task_id': self.task_id},
                    {'$set': {
                        "database_version": genome_doc["database_version"]
                    }}, upsert=True)

        ## 用于靶基因交互分析
        if self.option("target_predict"):
            col.update({"task_id": self.task_id},
                       {"$set": {"refrna_seqdb": os.path.join(self.workflow_output, "09Other/refrna_seqs.db")}},
                       upsert=True)
            col2 = db["sg_target"]
            col2.update({"task_id": self.task_id},
                        {"$set": {"result_dir": os.path.join(self.workflow_output, "07miRNA_Target")}},
                        upsert=True)
            if self.option("assembly_file").is_set:
                col.update({"task_id": self.task_id},
                           {"$set": {"assemble_fa": os.path.join(self.workflow_output, "09Other/upload.assembly.fa")}}, upsert=True)
                col5 = db["sg_annotation_stat"]
                col5.update({"task_id": self.task_id},
                            {"$set": {"result_dir": os.path.join(self.workflow_output, "09Other/Annotation/")}},
                            upsert=True)
            else:
                col.update({"task_id": self.task_id}, {"$set": {"assemble_fa": self.trans_fa}}, upsert=True)
            if self.option("gene_to_trans").is_set:
                col.update({"task_id": self.task_id},
                           {"$set": {"gene_to_trans": os.path.join(self.workflow_output, "09Other/upload.gene2trans.txt")}},
                           upsert=True)
            else:
                col.update({"task_id": self.task_id}, {"$set": {"gene_to_trans": self.g2t2p}}, upsert=True)
        ## 用于页面展示pdf图片
        if self.option("mirna_database").lower() != "none":
            col1 = db["sg_novel_mirna"]
            col1.update({"task_id": self.task_id},
                        {"$set": {
                            "pdf_dir": os.path.join(self.workflow_output, "04sRNA_Analysis/02Novel_miRNA/novel_pre_structure")}},
                        upsert=True)
            col2 = db["sg_known_mirna"]
            col2.update({"task_id": self.task_id},
                        {"$set": {
                            "pdf_dir": os.path.join(self.workflow_output, "04sRNA_Analysis/01Known_miRNA/known_pre_structure")}},
                        upsert=True)
        else:
            col1 = db["sg_novel_mirna"]
            col1.update({"task_id": self.task_id},
                        {"$set": {"pdf_dir": os.path.join(self.workflow_output, "04miRNA_Analysis/structure_pdf")}},
                        upsert=True)
        col4 = db["sg_mapping"]
        col4.update({"task_id": self.task_id}, {"$set": {"result_dir": os.path.join(self.workflow_output, "09Other/Mapping")}},
                    upsert=True)
        col4 = db["sg_circos"]
        col4.update({"task_id": self.task_id}, {"$set": {"result_dir": os.path.join(self.workflow_output, "09Other/Mapping")}},
                    upsert=True)
        self.modify_output()
        super(SmallrnaWorkflow, self).end()

    def modify_output(self):
        # 项目结果目录中呈现文件
        self.target_dir = os.path.join(self.work_dir, 'upload')
        self.db = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')]
        if os.path.isdir(self.target_dir):
            shutil.rmtree(self.target_dir)
        os.mkdir(self.target_dir)
        # 01Background
        os.mkdir(os.path.join(self.target_dir, '01Background'))
        ## get run parameter
        if self.option("get_run_log"):
            self.run_parameter(os.path.join(self.target_dir, '01Background'))
        ## genome_stat
        genome_stat = pd.read_table(self.genome_stat, header=0)
        genome_stat.rename(columns={'Chr': 'Chromosome', 'Size(Mb)': 'Length(MB)'}, inplace=True)
        genome_stat.to_csv(os.path.join(self.target_dir, '01Background', "genome_info.xls"),
                           sep="\t", header=True, index=False)
        ## software_info
        software_info = os.path.join(self.target_dir, '01Background', "software_info.xls")
        my_collection = self.db['sg_software_database']
        my_results = my_collection.find({})
        with open(software_info, "w") as w:
            w.write("\t".join(["Soft/Database", "Version", "Analysis", "Source"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["software_database"]), str(collection["version"]), str(collection["usage"]),
                     str(collection["source"])]) + "\n")
        ## sample_info
        sample_info = os.path.join(self.target_dir, '01Background', "sample_info.xls")
        if self.option('sample_num') == 'multiple':
            group_table = self.option('group_table').prop["path"]
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
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
        with open(sample_info, "w") as w, open(group_table, "r") as f:
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
                        w.write("\t".join(
                            [line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")
        # 02QC
        os.mkdir(os.path.join(self.target_dir, '02QC'))
        if self.option("datatype") == "rawdata":
            fq_stat_before = self.output_dir + "/QC/stat_rawdata/fastq_stat.xls"
            os.link(fq_stat_before, self.target_dir + "/02QC/rawdata_stat.xls")
        fq_stat_after = self.output_dir + "/QC/stat_cleandata/fastq_stat.xls"
        os.link(fq_stat_after, self.target_dir + "/02QC/cleandata_stat.xls")
        # 03Align
        os.mkdir(os.path.join(self.target_dir, '03Align'))
        main_id = self.db['sg_mapping'].find_one({"task_id": self.task_id})["main_id"]
        stat_results = self.db['sg_mapping_stat'].find({"mapping_id": main_id})
        stat_file = os.path.join(self.target_dir, '03Align', 'align_stat.xls')
        with open(stat_file, "w") as w:
            w.write("\t".join(["Sample", "Total reads", "Total mapped", "Mapped reads (+)", "Mapped reads (-)"]) + "\n")
            for result in stat_results:
                w.write("\t".join(
                    [str(result["sample_name"]), str(result["total_reads"]), str(result["mapped_reads"]),
                     str(result["mapped_reads_for"]), str(result["mapped_reads_rec"])]) + "\n")
        distribution_results = self.db['sg_mapping_detail'].find({"mapping_id": main_id})
        distribution_file = os.path.join(self.target_dir, '03Align', 'chr_distribution.xls')
        sample_file = self.qc.work_dir + "/list.txt"
        sample_list = list()
        with open(sample_file, 'r') as f:
            for line in f.readlines():
                sample = line.strip().split("\t")[1]
                if sample in sample_list:
                    pass
                else:
                    sample_list.append(sample)
        with open(distribution_file, "w") as w:
            w.write("Chromosome" + "\t".join(sample_list) + "\n")
            for result in distribution_results:
                result_list = list()
                result_list.append(result["ref"])
                for sample in sample_list:
                    result_list.append(str(result[sample + "_total"]))
                w.write("\t".join(result_list) + "\n")
        os.mkdir(os.path.join(self.target_dir, '03Align', 'chr_circos'))
        for file in os.listdir(self.output_dir + "/Mapping"):
            if "_circos.svg" in file:
                file_path = os.path.join(self.output_dir, "Mapping", file)
                os.link(file_path, self.target_dir + "/03Align/chr_circos/" + file)
        # 04sRNA_Analysis
        if self.option("mirna_database").lower() != "none":
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis'))
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA'))
            known_miRNA_detail = os.path.join(self.srna.output_dir, 'known_mirna', 'known_mirna_detail.xls')
            os.link(known_miRNA_detail,
                    os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA', 'known_miRNA_detail.xls'))
            known_mature = os.path.join(self.srna.output_dir, 'known_mirna', 'mature.fa')
            os.link(known_mature, os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA', 'known_mature.fa'))
            known_hairpin = os.path.join(self.srna.output_dir, 'known_mirna', 'hairpin.fa')
            os.link(known_hairpin,
                    os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA', 'known_hairpin.fa'))
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA', 'known_pre_structure'))
            structure_pdf = os.path.join(self.srna.output_dir, 'known_mirna', 'structure_pdf')
            for file in os.listdir(structure_pdf):
                file_path = os.path.join(self.srna.output_dir, 'known_mirna', 'structure_pdf', file)
                os.link(file_path,
                        os.path.join(self.target_dir, '04sRNA_Analysis', '01Known_miRNA', 'known_pre_structure', file))
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA'))
            novel_mirna_detail = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_mirna_detail.xls')
            os.link(novel_mirna_detail,
                    os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA', 'novel_miRNA_detail.xls'))
            novel_mature = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_mature_seq.fa')
            os.link(novel_mature, os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA', 'novel_mature.fa'))
            novel_hairpin = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_precursor_seq.fa')
            os.link(novel_hairpin,
                    os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA', 'novel_hairpin.fa'))
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA', 'novel_pre_structure'))
            structure_pdf = os.path.join(self.srna.output_dir, 'novel_mirna', 'structure_pdf')
            for file in os.listdir(structure_pdf):
                file_path = os.path.join(self.srna.output_dir, 'novel_mirna', 'structure_pdf', file)
                os.link(file_path,
                        os.path.join(self.target_dir, '04sRNA_Analysis', '02Novel_miRNA', 'novel_pre_structure', file))
            os.mkdir(os.path.join(self.target_dir, '04sRNA_Analysis', '03sRNA_stat'))
            miRNA_stat = os.path.join(self.srna.output_dir, 'srna_stat', 'mirna_stat.xls')
            os.link(miRNA_stat, os.path.join(self.target_dir, '04sRNA_Analysis', '03sRNA_stat', 'miRNA_stat.xls'))
            sRNA_stat = os.path.join(self.srna.output_dir, 'srna_stat', 'srna_stat.xls')
            os.link(sRNA_stat, os.path.join(self.target_dir, '04sRNA_Analysis', '03sRNA_stat', 'sRNA_stat.xls'))
        else:
            os.mkdir(os.path.join(self.target_dir, '04miRNA_Analysis'))
            miRNA_detail = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_mirna_detail.xls')
            os.link(miRNA_detail, os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_detail.xls'))
            miRNA_stat = os.path.join(self.srna.output_dir, 'srna_stat', 'mirna_stat.xls')
            os.link(miRNA_stat, os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_stat.xls'))
            sRNA_stat = os.path.join(self.srna.output_dir, 'srna_stat', 'srna_stat.xls')
            os.link(sRNA_stat, os.path.join(self.target_dir, '04miRNA_Analysis', 'sRNA_stat.xls'))
            miRNA_mature = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_mature_seq.fa')
            os.link(miRNA_mature, os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_mature.fa'))
            miRNA_hairpin = os.path.join(self.srna.output_dir, 'novel_mirna', 'novel_precursor_seq.fa')
            os.link(miRNA_hairpin, os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_hairpin.fa'))
            os.mkdir(os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_pre_structure'))
            structure_pdf = os.path.join(self.srna.output_dir, 'novel_mirna', 'structure_pdf')
            for file in os.listdir(structure_pdf):
                file_path = os.path.join(self.srna.output_dir, 'novel_mirna', 'structure_pdf', file)
                os.link(file_path, os.path.join(self.target_dir, '04miRNA_Analysis', 'miRNA_pre_structure', file))
        # 05Express
        os.mkdir(os.path.join(self.target_dir, '05Express'))
        if self.option("mirna_database").lower() != "none":
            os.mkdir(os.path.join(self.target_dir, '05Express', '01Exp_Annalysis'))
            All_miRNA_TPM = os.path.join(self.srna.output_dir, 'total_mirna_norm.xls')
            os.link(All_miRNA_TPM, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'All_miRNA_TPM.xls'))
            known_miRNA_TPM = os.path.join(self.srna.output_dir, 'known_mirna_norm.xls')
            os.link(known_miRNA_TPM, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'known_miRNA_TPM.xls'))
            novel_miRNA_TPM = os.path.join(self.srna.output_dir, 'novel_mirna_norm.xls')
            os.link(novel_miRNA_TPM, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'novel_miRNA_TPM.xls'))
            All_miRNA_count = os.path.join(self.srna.output_dir, 'total_mirna_count.xls')
            os.link(All_miRNA_count, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'All_miRNA_count.xls'))
            known_miRNA_count = os.path.join(self.srna.output_dir, 'known_mirna_count.xls')
            os.link(known_miRNA_count, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'known_miRNA_count.xls'))
            novel_miRNA_count = os.path.join(self.srna.output_dir, 'novel_mirna_count.xls')
            os.link(novel_miRNA_count, os.path.join(self.target_dir, '05Express', '01Exp_Annalysis', 'novel_miRNA_count.xls'))
            if self.option("sample_num") == "multiple":
                if self.option("group_table").prop["sample_number"] > 2:
                    os.mkdir(os.path.join(self.target_dir, '05Express', '02Exp_Corr'))
                    sample_correlation = os.path.join(self.output_dir, 'ExpCorr', 'sample_correlation.xls')
                    os.link(sample_correlation,
                            os.path.join(self.target_dir, '05Express', '02Exp_Corr', 'sample_correlation.xls'))
                    os.mkdir(os.path.join(self.target_dir, '05Express', '03Exp_PCA'))
                    PCA = os.path.join(self.output_dir, 'ExpPca', 'PCA.xls')
                    os.link(PCA, os.path.join(self.target_dir, '05Express', '03Exp_PCA', 'PCA.xls'))
                    Explained_variance_ratio = os.path.join(self.output_dir, 'ExpPca', 'Explained_variance_ratio.xls')
                    with open(Explained_variance_ratio, "r") as f, open(
                            os.path.join(self.target_dir, '05Express', '03Exp_PCA', 'Explained_variance_ratio.xls'),
                            "w") as w:
                        w.write("\tProportion of Variance\n")
                        for line in f:
                            w.write(line)
        else:
            miRNA_count = os.path.join(self.srna.output_dir, 'novel_mirna_count.xls')
            os.link(miRNA_count, os.path.join(self.target_dir, '05Express', 'miRNA_count.xls'))
            miRNA_tpm = os.path.join(self.srna.output_dir, 'novel_mirna_norm.xls')
            os.link(miRNA_tpm, os.path.join(self.target_dir, '05Express', 'miRNA_tpm.xls'))
        # 06Diff_Express
        if self.option("sample_num") == "multiple":
            os.mkdir(os.path.join(self.target_dir, '06Diff_Express'))
            diff_files_stat = glob.glob(
                os.path.join(self.output_dir, 'DiffExpress',
                             '*_vs_*.{}.xls'.format(self.option("diff_method").lower())))
            df1 = pd.DataFrame()
            for diff_file_stat in diff_files_stat:
                os.link(diff_file_stat,
                        os.path.join(self.target_dir, '06Diff_Express', os.path.basename(diff_file_stat)))
                ctrl, test = \
                os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("diff_method").lower()))[
                    0].split('_vs_')
                fname = os.path.basename(diff_file_stat).split(".")[0]
                df_t = pd.read_table(diff_file_stat, index_col="seq_id")
                df_t.rename(columns={"log2fc": "{}_log2fc({}/{})".format(fname, test, ctrl),
                                     "fc": "{}_fc({}/{})".format(fname, test, ctrl),
                                     "pvalue": "{}_pvalue".format(fname),
                                     "padjust": "{}_padjust".format(fname),
                                     "significant": "{}_significant".format(fname),
                                     "regulate": "{}_regulate".format(fname)},
                            inplace=True)
                df_core = pd.DataFrame(df_t, columns=["{}_fc({}/{})".format(fname, test, ctrl),
                                                      "{}_log2fc({}/{})".format(fname, test, ctrl),
                                                      "{}_pvalue".format(fname), "{}_padjust".format(fname),
                                                      "{}_significant".format(fname), "{}_regulate".format(fname)])
                df1 = pd.concat([df1, df_core], axis=1)
            df1.index.set_names("seq_id", inplace=True)
            df1.to_csv(os.path.join(self.target_dir, '06Diff_Express',
                                    "total_diff_stat.{}.xls".format(self.option("diff_method").lower())), sep="\t")
            if glob.glob(os.path.join(self.output_dir, 'DiffExpress', '*_diff_summary.xls')):
                diff_stat = glob.glob(os.path.join(self.output_dir, 'DiffExpress', '*_diff_summary.xls'))[0]
            else:
                diff_stat = glob.glob(os.path.join(self.output_dir, 'DiffExpress', 'diff_summary_*.xls'))[0]
            os.link(diff_stat,
                    os.path.join(self.target_dir, '06Diff_Express',
                                 'diff_stat_{}.xls'.format(self.option("diff_method"))))
        # 07miRNA_Target
        if self.option("target_predict"):
            os.mkdir(os.path.join(self.target_dir, '07miRNA_Target'))
            if self.option("mirna_database").lower() != "none":
                known_detail_gz = glob.glob(self.output_dir + "/MirnaTarget/known_*_detail.txt.gz")[0]
                os.link(known_detail_gz,
                        os.path.join(self.target_dir, '07miRNA_Target', os.path.basename(known_detail_gz)))
            novel_detail_gz = glob.glob(self.output_dir + "/MirnaTarget/novol_*_detail.txt.gz")[0]
            os.link(novel_detail_gz, os.path.join(self.target_dir, '07miRNA_Target',
                                                  os.path.basename(novel_detail_gz).replace("novol", "novel")))
            target_fa = os.path.join(self.output_dir, 'MirnaTarget', 'target.fa')
            os.link(target_fa, os.path.join(self.target_dir, '07miRNA_Target', os.path.basename(target_fa)))
            known_target = os.path.join(self.output_dir, 'MirnaTarget', 'known_target.xls')
            novol_target = os.path.join(self.output_dir, 'MirnaTarget', 'novol_target.xls')
            target_predict_detail = os.path.join(self.target_dir, '07miRNA_Target', 'target_predict_detail.xls')
            target_predict_stat = os.path.join(self.target_dir, '07miRNA_Target', 'target_predict_stat.xls')
            with open(target_predict_detail, "w") as w, open(known_target, "r") as f1, open(novol_target, "r") as f2:
                for line in f1:
                    w.write(line)
                head = f2.readline()
                for line in f2:
                    w.write(line)
            target_stat = self.db["sg_target_stat"]
            with open(target_predict_stat, "w") as w:
                w.write("\t".join(["Type", "All miRNA", "miRNA with Target", "Target"]))
                result = target_stat.find_one({"type": "Known miRNA"})
                w.write("\t".join(["Known miRNA", str(result['all_mirna']), str(result['mirna_target']), str(
                    result['target'])]) + "\n")
                result = target_stat.find_one({"type": "Novel miRNA"})
                w.write("\t".join(["Novel miRNA", str(result['all_mirna']), str(result['mirna_target']), str(
                    result['target'])]) + "\n")
                result = target_stat.find_one({"type": "Total"})
                w.write("\t".join(["Total", str(result['all_mirna']), str(result['mirna_target']), str(
                    result['target'])]) + "\n")
            target_anno_detail = os.path.join(self.output_dir, 'MirnaTarget', 'All_annot_target.xls')
            os.link(target_anno_detail, os.path.join(self.target_dir, '07miRNA_Target', 'target_anno_detail.xls'))
            target_anno_stat_main_id = self.db["sg_annotation_stat"].find_one({"task_id": self.task_id})["main_id"]
            target_anno_stat = os.path.join(self.output_dir, 'MirnaTarget', 'target_anno_stat.xls')
            with open(target_anno_stat, "w") as w:
                w.write("\t".join(["", "target number", "target percent"]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "NR"})
                w.write("\t".join(["NR", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "Swiss-Prot"})
                w.write("\t".join(["Swiss-Prot", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "Pfam"})
                w.write("\t".join(["Pfam", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "KEGG"})
                w.write("\t".join(["KEGG", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "GO"})
                w.write("\t".join(["GO", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "COG"})
                w.write("\t".join(["COG", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "Total_anno"})
                w.write("\t".join(["Total_anno", str(result['gene']), str(result['gene_percent'])]) + "\n")
                result = self.db["sg_annotation_stat_detail"].find_one(
                    {"stat_id": target_anno_stat_main_id, "type": "Total"})
                w.write("\t".join(["Total", str(result['gene']), str(result['gene_percent'])]) + "\n")
        # 08miRNA_Structure
        if self.option("mirna_database").lower() != "none":
            os.mkdir(os.path.join(self.target_dir, '08miRNA_Structure'))
            os.mkdir(os.path.join(self.target_dir, '08miRNA_Structure', '01miRNA_bias'))
            first_bias_per = os.path.join(self.output_dir, 'MirnaBias', 'all_first_bias_per.xls')
            os.link(first_bias_per,
                    os.path.join(self.target_dir, '08miRNA_Structure', '01miRNA_bias', 'first_bias_per.xls'))
            loc_bias_per = os.path.join(self.output_dir, 'MirnaBias', 'all_loc_bias_per.xls')
            os.link(loc_bias_per,
                    os.path.join(self.target_dir, '08miRNA_Structure', '01miRNA_bias', 'loc_bias_per.xls'))
            os.mkdir(os.path.join(self.target_dir, '08miRNA_Structure', '02miRNA_family'))
            family_species = os.path.join(self.output_dir, 'MirnaFamily', 'family.species.xls')
            os.link(family_species,
                    os.path.join(self.target_dir, '08miRNA_Structure', '02miRNA_family', 'family_species.xls'))
            if self.option("is_miRNA_edit"):
                os.mkdir(os.path.join(self.target_dir, '08miRNA_Structure', '03miRNA_edit'))
                for file in os.listdir(self.mirna_edit.output_dir):
                    file_path = os.path.join(self.mirna_edit.output_dir, file)
                    if os.path.basename(file_path) == "all.result.xls":
                        os.link(file_path,
                                os.path.join(self.target_dir, '08miRNA_Structure', '03miRNA_edit', "edit_detail.xls", ))
        # 09Other
        os.mkdir(os.path.join(self.target_dir, '09Other'))
        ## 上传的序列文件
        if self.option("assembly_file").is_set:
            os.link(self.option("assembly_file").prop["path"], self.target_dir + "/09Other/upload.assembly.fa")
        if self.option("gene_to_trans").is_set:
            os.link(self.option("gene_to_trans").prop["path"], self.target_dir + "/09Other/upload.gene2trans.txt")
        # os.mkdir(self.target_dir + "/09Other/rawdata")
        # for file in os.listdir(os.path.join(self.qc.output_dir, "raw_data")):
        #     file_path = os.path.join(self.qc.output_dir, "raw_data", file)
        #     if os.path.basename(file_path) != "list.txt":
        #         os.link(file_path, self.target_dir + "/09Other/rawdata/" + file)
        ## MirnaTarget
        if self.option("target_predict"):
            seq_db = self.output_dir + "/Sequence_database/refrna_seqs.db"
            os.link(seq_db, self.target_dir + "/09Other/refrna_seqs.db")
            if self.option("assembly_file").is_set:
                os.mkdir(self.target_dir + '/09Other/Annotation')
                CopyFile().linkdir(self.output_dir + "/Annotation", os.path.join(self.target_dir, "09Other/Annotation"))
                if os.path.exists(os.path.join(self.target_dir, "09Other/Annotation/tran2gene.txt")):
                    os.remove(os.path.join(self.target_dir, "09Other/Annotation/tran2gene.txt"))
                os.link(os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map"),
                        os.path.join(self.target_dir, "09Other/Annotation/tran2gene.txt"))
            # os.mkdir(self.target_dir + "/09Other/miRNA_Target")
            # for file in os.listdir(self.output_dir + "/MirnaTarget"):
            #     file_path = os.path.join(self.origin_dir, "MirnaTarget", file)
            #     os.link(file_path, os.path.join(self.target_dir, "09Other/miRNA_Target", file))
        ## Mapping
        os.mkdir(self.target_dir + "/09Other/Mapping")
        for file in os.listdir(self.output_dir + "/Mapping"):
            file_path = os.path.join(self.output_dir, "Mapping", file)
            os.link(file_path, self.target_dir + "/09Other/Mapping/" + file)
        ## seq_download
        os.mkdir(self.target_dir + "/09Other/SequenceDetail")
        cds_seq = os.path.join(self.sequence_detail.work_dir, 'cds_seq')
        print "cds_seq_path"
        print cds_seq
        pep_seq = os.path.join(self.sequence_detail.work_dir, 'pep_seq')
        txpt_seq = os.path.join(self.sequence_detail.work_dir, 'txpt_seq')
        gene_detail = os.path.join(self.sequence_detail.work_dir, 'gene_detail')
        trans_detail = os.path.join(self.sequence_detail.work_dir, 'tran_detail')
        gene_stat = os.path.join(self.sequence_detail.work_dir, 'gene_stat')
        # trans_stat = os.path.join(self.sequence_detail.work_dir, 'tran_stat')
        os.link(cds_seq, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(cds_seq))
        os.link(pep_seq, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(pep_seq))
        os.link(txpt_seq, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(txpt_seq))
        os.link(gene_detail, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(gene_detail))
        os.link(trans_detail, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(trans_detail))
        os.link(gene_stat, self.target_dir + '/09Other/SequenceDetail/' + os.path.basename(gene_stat))
        # os.link(trans_stat, os.path.join(self.target_dir, '/09Other/SequenceDetail', os.path.basename(trans_stat)))
        ##Assemble




        sdir = self.add_upload_dir(self.target_dir)
        sdir.add_regexp_rules([
            [r"03Align/chr_circos/.*_circos\.svg", "svg", "不同样本染色体Reads分布circos图", 0],
            [r"04sRNA_Analysis/01Known_miRNA/known_pre_structure/.*\.pdf", "pdf", "已知miRNA前体结构图", 0],
            [r"04sRNA_Analysis/02Novel_miRNA/novel_pre_structure/.*\.pdf", "pdf", "新miRNA前体结构图", 0],
            [r"04miRNA_Analysis/miRNA_pre_structure/.*\.pdf", "pdf", "miRNA前体结构图", 0],
            [r"06Diff_Express/diff_stat_.*", "xls", "差异表达miRNA统计表", 0],
            [r"06Diff_Express/total_diff_stat_.*", "xls", "差异表达miRNA详情表（所有比较组）", 0],
            [r"06Diff_Express/.*_vs_.*\.xls", "xls", "差异表达miRNA详情表（单个比较组）", 0],
            [r"07miRNA_Target/known_.*_detail\.txt\.gz", "gz", "已知miRNA靶基因比对详细信息", 0],
            [r"07miRNA_Target/novel_.*_detail\.txt\.gz", "gz", "新miRNA靶基因比对详细信息", 0],
        ])
        sdir.add_relpath_rules([
            [".", "", "流程分析结果目录", 0],
            ["01Background", "", "项目背景结果目录", 0],
            ["01Background/genome_info.xls", "xls", "基因组注释信息表", 0],
            ["01Background/sample_info.xls", "xls", "样本信息表", 0],
            ["01Background/software_info.xls", "xls", "软件信息表", 0],
            ['01Background/run_parameter.txt', 'txt', '分析参数日志', 0],
            ["02QC", "", "测序数据质控结果目录", 0],
            ["02QC/rawdata_stat.xls", "xls", "原始数据统计表", 0],
            ["02QC/cleandata_stat.xls", "xls", "质控数据统计表", 0],
            ["03Align", "", "序列比对分析结果目录", 0],
            ["03Align/align_stat.xls", "xls", "比对结果统计表", 0],
            ["03Align/chr_distribution.xls", "xls", "染色体Reads分布统计表", 0],
            ["03Align/chr_circos", "", "染色体Reads分布circos图", 0],
            ["04sRNA_Analysis", "", "sRNA分析结果目录", 0],
            ["04sRNA_Analysis/01Known_miRNA", "", "已知miRNA分析", 0],
            ["04sRNA_Analysis/01Known_miRNA/known_miRNA_detail.xls", "xls", "已知miRNA详情表", 0],
            ["04sRNA_Analysis/01Known_miRNA/known_mature.fa", "fasta", "已知miRNA成熟体序列", 0],
            ["04sRNA_Analysis/01Known_miRNA/known_hairpin.fa", "fasta", "已知miRNA前体序列", 0],
            ["04sRNA_Analysis/01Known_miRNA/known_pre_structure", "", "已知miRNA前体结构文件", 0],
            ["04sRNA_Analysis/02Novel_miRNA", "", "新miRNA预测", 0],
            ["04sRNA_Analysis/02Novel_miRNA/novel_miRNA_detail.xls", "xls", "新miRNA详情表", 0],
            ["04sRNA_Analysis/02Novel_miRNA/novel_mature.fa", "fasta", "新miRNA成熟体序列", 0],
            ["04sRNA_Analysis/02Novel_miRNA/novel_hairpin.fa", "fasta", "新miRNA前体序列", 0],
            ["04sRNA_Analysis/02Novel_miRNA/novel_pre_structure", "", "新miRNA前体结构文件", 0],
            ["04sRNA_Analysis/03sRNA_stat", "", "sRNA统计", 0],
            ["04sRNA_Analysis/03sRNA_stat/miRNA_stat.xls", "xls", "各样本miRNA统计表", 0],
            ["04sRNA_Analysis/03sRNA_stat/sRNA_stat.xls", "xls", "各样本sRNA统计表", 0],
            ["04miRNA_Analysis", "", "miRNA分析结果目录", 0],
            ["04miRNA_Analysis/miRNA_detail.xls", "xls", "miRNA详情表", 0],
            ["04miRNA_Analysis/miRNA_stat.xls", "xls", "各样本miRNA统计表", 0],
            ["04miRNA_Analysis/sRNA_stat.xls", "xls", "各样本sRNA统计表", 0],
            ["04miRNA_Analysis/miRNA_mature.fa", "fasta", "miRNA成熟体序列", 0],
            ["04miRNA_Analysis/miRNA_hairpin.fa", "fasta", "miRNA前体序列", 0],
            ["04miRNA_Analysis/miRNA_pre_structure", "", "miRNA前体结构文件", 0],
            ["05Express", "", "miRNA表达量分析结果目录", 0],
            ["05Express/miRNA_count.xls", "", "count表达量分析结果表", 0],
            ["05Express/miRNA_tpm.xls", "", "TPM表达量分析结果表", 0],
            ["05Express/01Exp_Annalysis", "", "表达定量结果", 0],
            ["05Express/01Exp_Annalysis/All_miRNA_TPM.xls", "xls", "miRNA表达量总表（TPM）", 0],
            ["05Express/01Exp_Annalysis/known_miRNA_TPM.xls", "xls", "已知miRNA表达量表（TPM）", 0],
            ["05Express/01Exp_Annalysis/novel_miRNA_TPM.xls", "xls", "新miRNA表达量表（TPM）", 0],
            ["05Express/01Exp_Annalysis/All_miRNA_count.xls", "xls", "miRNA表达量总表（count）", 0],
            ["05Express/01Exp_Annalysis/known_miRNA_count.xls", "xls", "已知miRNA表达量表（count）", 0],
            ["05Express/01Exp_Annalysis/novel_miRNA_count.xls", "xls", "新miRNA表达量表（count）", 0],
            ["05Express/01Exp_Annalysis/miRNA_count.xls", "xls", "count表达量分析结果表", 0],
            ["05Express/01Exp_Annalysis/miRNA_tpm.xls", "xls", "TPM表达量分析结果表", 0],
            ["05Express/02Exp_Corr", "", "样本间相关性分析", 0],
            ["05Express/02Exp_Corr/sample_correlation.xls", "xls", "样本间相关系数表", 0],
            ["05Express/03Exp_PCA", "", "样本间PCA分析", 0],
            ["05Express/03Exp_PCA/PCA.xls", "xls", "PCA分析结果表", 0],
            ["05Express/03Exp_PCA/Explained_variance_ratio.xls", "xls", "主成分解释表", 0],
            ["06Diff_Express", "", "miRNA表达量差异分析结果目录", 0],
            ["07miRNA_Target", "", "miRNA靶基因预测结果目录", 0],
            ["07miRNA_Target/target_predict_detail.xls", "xls", "靶基因预测详情表", 0],
            ["07miRNA_Target/target_predict_stat.xls", "xls", "靶基因预测统计表 ", 0],
            ["07miRNA_Target/target_anno_detail.xls", "xls", "靶基因注释详情表", 0],
            ["07miRNA_Target/target_anno_stat.xls", "xls", "靶基因注释统计表 ", 0],
            ["07miRNA_Target/target.fa", "fasta", "靶基因序列", 0],
            ["08miRNA_Structure", "", "miRNA结构分析结果目录", 0],
            ["08miRNA_Structure/01miRNA_bias", "", "miRNA碱基偏好性分析", 0],
            ["08miRNA_Structure/01miRNA_bias/first_bias_per.xls", "xls", "不同长度miRNA首位碱基偏好性统计表", 0],
            ["08miRNA_Structure/01miRNA_bias/loc_bias_per.xls", "xls", "miRNA不同位点碱基的偏好性统计表", 0],
            ["08miRNA_Structure", "", "miRNA结构分析结果目录", 0],
            ["08miRNA_Structure/02miRNA_family", "", "miRNA家族分析", 0],
            ["08miRNA_Structure/02miRNA_family/family_species.xls", "xls", "miRNA的家族信息表", 0],
            ["08miRNA_Structure/03miRNA_edit", "", "miRNA碱基编辑分析", 0],
            ["08miRNA_Structure/03miRNA_edit/edit_detail.xls", "", "miRNA 碱基编辑详情表", 0],
            ["09Other", "", "页面交互分析用中间文件", 0],
        ])

    def run_parameter(self, dir_path):
        document = self.db['sg_table_relation'].find_one({})
        targets = document['target']
        parsed_tables = list()
        for target in targets:
            if target[0] in ['sg_status', 'sg_task', 'sg_software_para', 'sg_geneset', 'task_workflow_params',
                             'sg_one_touch', 'sg_file_dowload', 'sg_transcripts', 'sg_result_table_deposit',
                             'sg_splicing_rmats_count'] or target[0] in parsed_tables:
                continue
            parsed_tables.append(target[0])
            records = self.db[target[0]].find({'task_id': self.task_id})
            for record in records:
                if 'params' in record and record['params'] and record['params'] != 'none':
                    if type(record['params']) == dict:
                        params = json.loads(record['params'])
                    else:
                        params = record['params']
                    if 'submit_location' in params:
                        get_run_log = GetRunLog("small_rna", table=target[0], main_id=record['main_id'],
                                                dir_path=dir_path, append=True)
                        get_run_log.run()

    def run_api(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.smallrna_task_info')
        task_info.add_task_info()
        self.logger.info("导表开始")
        self.stop_timeout_check()
        self.export_genome_info()
        self.export_qc()
        self.export_mapping()
        self.export_srna()
        if self.option("target_predict"):
            self.export_target()
        self.export_expression()
        if not self.option("assembly_file").is_set and self.option("target_predict"):
            self.export_gene_detail()
        else:
            if self.option("target_predict"):
                self.build_seq_database()
        if self.option("mirna_database").lower() != "none":
            self.export_bias()
            self.export_family()
            if self.option("is_miRNA_edit"):
                self.export_mirna_edit()
        self.logger.info("导表完成")

    @time_count
    def export_mapping(self):
        project_sn = self.project_sn
        task_id = self.task_id
        group_id = self.group_id
        group_dict = self.option('group_table').prop['group_dict']

        self.api_mapping = self.api.api("small_rna.small_rna_mapping")
        sample_file = self.qc.work_dir + "/list.txt"
        map_dir = self.mapping.output_dir
        params = {
            'task_id': task_id,
            'submit_location': 'mapping_circs',
            'task_type': 2,
            'group_id': str(group_id),
            'group_dict': group_dict,
            'chr_num': 10
        }
        self.api_mapping.add_map_stat(sample_list=None, map_dir=map_dir, sample_file=sample_file, params=params, group=self.option("group_table").path)

    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("small_rna.genome_info")
        file_path = self.genome_stat
        species_name = self.species_name
        species = self.species
        ref_anno_version = self.ref_anno_version
        hyperlink = self.hyperlink
        self.api_geno.add_genome_info(file_path=file_path, species_name=species_name, species=species,
                                      ref_anno_version=ref_anno_version, hyperlink=hyperlink)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("small_rna.small_rna_qc")
        fq_type = "SE"
        quality_stat_after = self.qc_stat_after.output_dir
        quality_stat_before = self.qc_stat_before.output_dir
        self.api_qc.add_before_qc(quality_stat_before, fq_type=fq_type, about_qc="before", group=self.option("group_table").path)
        qc_stat = self.qc.output_dir + "/clean_data"
        self.api_qc.add_after_qc(quality_stat_after, qc_dir=qc_stat, fq_type=fq_type, about_qc="after", group=self.option("group_table").path)
        if self.option('group_table').is_set:
            if self.option('productive_table').is_set:
                self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(
                    self.option("group_table").prop["path"], productive_table=self.option('productive_table').path)
            else:
                self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(
                    self.option("group_table").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
        if self.option('control_file').is_set:
            self.control_id, self.compare_detail = self.api_qc.add_control_group(
                self.option("control_file").prop["path"],
                self.group_id)

    @time_count
    def export_srna(self):
        self.srna_api = self.api.api("small_rna_v2.srna")
        novel_mirna = os.path.join(self.srna.output_dir, "novel_mirna/novel_mirna_detail.xls")
        pdfs_novel = os.path.join(self.srna.output_dir, "novel_mirna/structure_pdf")
        ncrna_stat = os.path.join(self.srna.output_dir, "srna_stat/ncrna_stat.xls")
        srna_stat = os.path.join(self.srna.output_dir, "srna_stat/srna_stat.xls")
        mirna_stat = os.path.join(self.srna.output_dir, "srna_stat/mirna_stat.xls")
        srna_stat_for_graph = os.path.join(self.srna.output_dir, "srna_stat/srna_stat_for_graph.xls")
        category = self.option("taxonmy").lower()
        method = self.option("mirna_soft")
        if self.option("mirna_database").lower() != "none":
            known_mirna = os.path.join(self.srna.output_dir, "known_mirna/known_mirna_detail.xls")
            pdfs_known = os.path.join(self.srna.output_dir, "known_mirna/structure_pdf")
            self.srna_api.add_known_mirna(known_mirna, project_sn=self.project_sn, task_id=self.task_id,
                                          params={"method": "quantifier"}, pdfs=pdfs_known)
        # self.srna_api.add_ncrna_stat(ncrna_stat, project_sn=self.project_sn, task_id=self.task_id,
        #                              params={"method": "statistics"})
        self.srna_api.add_novel_mirna(novel_mirna, project_sn=self.project_sn, task_id=self.task_id,
                                      params={"method": method}, method=method, pdfs=pdfs_novel)
        self.srna_api.add_mirna_stat(mirna_stat, project_sn=self.project_sn, task_id=self.task_id, group=self.option("group_table").path,
                                     params={"method": method})
        self.srna_api.add_srna_stat(srna_stat, srna_stat_for_graph, project_sn=self.project_sn, task_id=self.task_id,
                                    params={"method": method}, group=self.option("group_table").path)

    @time_count
    def export_target(self):
        self.target = self.api.api("small_rna.target_annotation")
        if self.option("target_source") == "mongo":
            result_dir = self.download.option("anno_class").prop['path']
            g2t2p = self.download.option("assembly_gene2trans").prop['path']
        elif self.option("target_source") == "ref":
            result_dir = self.annot_class
            g2t2p = self.g2t2p
        else:
            result_dir = self.annotation.output_dir
            g2t2p = os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map")
        target = self.mirna_target.output_dir
        species_name = self.species_name
        new_target_file = target + '/novol_target.xls'
        known_target_file = target + '/known_target.xls'
        if self.option("mirna_database").lower() == "none":
            os.system('touch {}'.format(self.work_dir + '/mature.fa'))
            known_seq = self.work_dir + '/mature.fa'
        else:
            known_seq = self.srna.output_dir + "/known_mirna/mature.fa"
        new_seq = self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa"
        params = {
            "nr_evalue": 1e-5,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue": 1e-5,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-5,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-5,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-5,
        }
        if self.option("taxonmy").lower() == "animal":
            if "miRanda" in self.option("target_animal_predict_method"):
                miranda = "yes"
            else:
                miranda = "no"
            if "TargetScan" in self.option("target_animal_predict_method"):
                targetscan = "yes"
            else:
                targetscan = "no"
            if "RNAhybrid" in self.option("target_animal_predict_method"):
                rnahybrid = "yes"
            else:
                rnahybrid = "no"
            #modify bwy fwy 20210113
            params_target = {
                "miranda": miranda,
                "targetscan": targetscan,
                "rnahybrid": rnahybrid,
                "min_support": str(self.option("min_support"))
            }
            if miranda == "yes":
                params_target.update({'miranda_score': "160.0",
                'miranda_energy': "-20",
                'miranda_strict': "on"})
            if rnahybrid == "yes":
                params_target.update({
                    'rnahybrid_num': "100",
                    'rnahybrid_energy': "-20",
                    'rnahybrid_pvalue': "0.01"
                })
            # params_target = {
            #     "miranda": miranda,
            #     "targetscan": targetscan,
            #     "rnahybrid": rnahybrid,
            #     "min_support": str(self.option("min_support")),
            #     'miranda_score': "160.0",
            #     'miranda_energy': "-20",
            #     'miranda_strict': "on",
            #     'rnahybrid_num': "100",
            #     'rnahybrid_energy': "-20",
            #     'rnahybrid_pvalue': "0.01",
            #     # 'ps_robot_score': "2.5",
            #     # 'targetfinder_score': "4"
            # }
        else:
            if "psRobot" in self.option("target_plant_predict_method"):
                psRobot = "yes"
            else:
                psRobot = "no"
            if "targetfinder" in self.option("target_plant_predict_method"):
                targetfinder = "yes"
            else:
                targetfinder = "no"
            if "RNAhybrid" in self.option("target_plant_predict_method"):
                RNAhybrid = "yes"
            else:
                RNAhybrid = "no"
            params_target = {
                "psrobot": psRobot,
                "targetfinder": targetfinder,
                "rnahybrid": RNAhybrid,
                "min_support": str(self.option("min_support"))
            }
            if psRobot ==  "yes":
                params_target.update({
                    'ps_robot_score': "2.5"
                })
            if targetfinder == "yes":
                params_target.update({
                    'targetfinder_score': "4"
                })
            if RNAhybrid == "yes":
                params_target.update({
                    'rnahybrid_num': "100",
                    'rnahybrid_energy': "-20",
                    'rnahybrid_pvalue': "0.01"
                })
            # params_target = {
            #     "psrobot": psRobot,
            #     "targetfinder": targetfinder,
            #     "rnahybrid": RNAhybrid,
            #     "min_support": str(self.option("min_support")),
            #     # 'miranda_score': "160.0",
            #     # 'miranda_energy': "-20",
            #     # 'miranda_strict': "on",
            #     'rnahybrid_num': "100",
            #     'rnahybrid_energy': "-20",
            #     'rnahybrid_pvalue': "0.01",
            #     'ps_robot_score': "2.5",
            #     'targetfinder_score': "4"

            # }

        self.target.run(known_target_file, new_target_file, result_dir, g2t2p, params, taxon=self.option("taxonmy"),
                        exp_level='transcript', version=self.option("version"))
        self.target.import_target_detail(new_target_file, known_target_file, params_target, new_seq=new_seq,
                                         known_seq=known_seq, anno_type="origin", species_name=species_name,
                                         target_dir=target, version=self.option("version"))
        if self.option("sample_num") == "multiple":
            if glob.glob(os.path.join(self.output_dir, 'DiffExpress', '*_diff_summary.xls')):
                diff_summary = glob.glob(os.path.join(self.output_dir, 'DiffExpress', '*_diff_summary.xls'))[0]
            else:
                diff_summary = glob.glob(os.path.join(self.output_dir, 'DiffExpress', 'diff_summary_*.xls'))[0]
            self.target.add_target_geneset(new_target_file, known_target_file, diff_summary)

    @time_count
    def export_expression(self):
        if self.option("mirna_database").lower() != "none":
            known_count_xls = os.path.join(self.srna.output_dir, "known_mirna_count.xls")
            known_tpm_xls = os.path.join(self.srna.output_dir, "known_mirna_norm.xls")
        novel_count_xls = os.path.join(self.srna.output_dir, "novel_mirna_count.xls")
        novel_tpm_xls = os.path.join(self.srna.output_dir, "novel_mirna_norm.xls")

        # set basic variables
        gevent.sleep()
        all_exp = self.api.api("small_rna_v2.all_exp")
        project_sn = self.project_sn
        task_id = self.task_id
        if self.option('group_table').is_set:
            group_id = self.group_id
            group_dict = self.option('group_table').prop['group_dict']
        if self.option('control_file').is_set:
            control_id = str(self.control_id)

        # create main_table in sg_exp and add detail_table in sg_exp_detail for count
        if self.option("mirna_database").lower() != "none":
            count_main_id = all_exp.add_exp(
                exp_matrix=known_count_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='count',
                is_novel=False
            )
            count_exp_id = all_exp.add_exp(
                exp_matrix=novel_count_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='count',
                is_novel=True,
                main_id=str(count_main_id)
            )
        else:
            count_exp_id = all_exp.add_exp(
                exp_matrix=novel_count_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='count',
                is_novel=True
            )

        # create main_table in sg_exp and add detail_table in sg_exp_detail for tpm
        if self.option("mirna_database").lower() != "none":
            tpm_main_id = all_exp.add_exp(
                exp_matrix=known_tpm_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='tpm',
                is_novel=False
            )
            tpm_exp_id = all_exp.add_exp(
                exp_matrix=novel_tpm_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='tpm',
                is_novel=True,
                main_id=str(tpm_main_id)
            )
        else:
            tpm_exp_id = all_exp.add_exp(
                exp_matrix=novel_tpm_xls,
                project_sn=project_sn,
                task_id=task_id,
                exp_type='tpm',
                is_novel=True
            )

        # create main_table in sg_exp_graph and add detail_table in sg_exp_graph_*
        all_tpm_xls = os.path.join(self.work_dir, 'all_mirna_norm.xls')
        if self.option("mirna_database").lower() != "none":
            all_tpm_pd = pd.concat([
                all_exp.process_exp_matrix(known_tpm_xls),
                all_exp.process_exp_matrix(novel_tpm_xls)
            ], axis=0)
        else:
            all_tpm_pd = all_exp.process_exp_matrix(novel_tpm_xls)
        all_tpm_pd.to_csv(all_tpm_xls, sep='\t', header=True, index=True)
        if self.option("mirna_database").lower() != "none":
            params = {
                'task_id': task_id,
                'submit_location': 'expgraph',
                'task_type': 2,
                'group_id': str(group_id),
                'group_dict': group_dict,
                'exp_id': str(tpm_exp_id),
                'seq_type': 'all',
            }
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            all_exp.add_distribution(
                exp_matrix=all_tpm_xls,
                group_dict=group_dict,
                project_sn=project_sn,
                task_id=task_id,
                exp_id=str(tpm_exp_id),
                seq_type='all',
                exp_type='tpm',
                params=params
            )

            # create main_table in sg_exp_corr and add detail_table in sg_exp_corr_detail
            if self.option("sample_num") == "multiple":
                if self.option("group_table").prop["sample_number"] > 2:
                    corr_output_dir = self.exp_corr.output_dir
                    params = {
                        'task_id': task_id,
                        'submit_location': 'expcorr',
                        'task_type': 2,
                        'group_id': str(group_id),
                        'group_dict': group_dict,
                        'exp_id': str(tpm_exp_id),
                        'corr_method': 'pearson',
                        'scm': 'complete',
                        'scd': 'euclidean',
                        'use_log': 'no',
                    }
                    all_exp.add_exp_corr(
                        corr_output_dir=corr_output_dir,
                        params=params,
                        task_id=task_id,
                        project_sn=project_sn
                    )

            # create main_table in sg_exp_pca and add detail_table in sg_exp_pca_detail
            if self.option("sample_num") == "multiple":
                if self.option("group_table").prop["sample_number"] > 2:
                    pca_output_dir = self.exp_pca.output_dir
                    params = {
                        'task_id': task_id,
                        'submit_location': 'exppca',
                        'task_type': 2,
                        'group_id': str(group_id),
                        'group_dict': group_dict,
                        'exp_id': str(tpm_exp_id),
                    }
                    main_id = all_exp.add_exp_pca(
                            pca_output_dir=pca_output_dir,
                            params=params,
                            task_id=task_id,
                            project_sn=project_sn
                        )
                    if hasattr(self, 'ellipse'):
                        if os.path.exists(os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')):
                            all_exp.insert_ellipse_table(os.path.join(self.ellipse.work_dir, 'ellipse_out.xls'), main_id)

        # create main_table in sg_diff and add detail_table in sg_diff_*
        if self.option("sample_num") == "multiple":
            if len(glob.glob(self.diffexpress.output_dir + "/*")) > 2:
                diff_output_dir = self.diffexpress.output_dir
                diff_method = self.option('diff_method')
                if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                    params = {
                        'task_id': task_id,
                        'submit_location': 'diff_detail',
                        'task_type': 2,
                        'group_id': str(group_id),
                        'group_dict': group_dict,
                        'control_id': str(control_id),
                        'exp_id': str(tpm_exp_id),
                        'diff_method': self.option("diff_method"),
                        'stat_type': self.option("pvalue_padjust").lower(),
                        'stat_cutoff': str(self.option("diff_fdr_ci")),
                        'fc': str(self.option("fc")),
                        'correct_method': self.option("padjust_way"),
                        'is_batch': 'False',
                    }
                    params = json.dumps(params, sort_keys=True, separators=(',', ':'))
                    all_exp.add_diffexp(
                        diff_output=diff_output_dir,
                        exp_id=str(tpm_exp_id),
                        project_sn=project_sn,
                        task_id=task_id,
                        params=params,
                        create_geneset=True,
                        diff_method=self.option("diff_method"),
                        group_dict=group_dict
                    )
                if diff_method.lower() in ['noiseq']:
                    params = {
                        'task_id': task_id,
                        'submit_location': 'diff_detail',
                        'task_type': 2,
                        'group_id': str(group_id),
                        'group_dict': group_dict,
                        'control_id': str(control_id),
                        'exp_id': str(tpm_exp_id),
                        'diff_method': self.option("diff_method"),
                        'fc': str(self.option("fc")),
                        'stat_cutoff': str(self.option("diff_fdr_ci")),
                        'prob': float(self.option('diff_fdr_ci')),
                        'is_batch': 'False',

                    }
                    params = json.dumps(params, sort_keys=True, separators=(',', ':'))
                    all_exp.add_diffexp_noiseq(
                        diff_output=diff_output_dir,
                        exp_id=str(tpm_exp_id),
                        project_sn=project_sn,
                        task_id=task_id,
                        params=params,
                        create_geneset=True,
                        diff_method=self.option("diff_method"),
                        group_dict=group_dict
                    )

    @time_count
    def export_gene_detail(self):
        """
        导入基因详情表
        :return:
        """
        gevent.sleep()
        self.api_gene_detail = self.api.api('small_rna.add_gene_detail')
        db_path = self.output_dir + "/Sequence_database/refrna_seqs.db"
        if os.path.exists(self.output_dir + "/Sequence_database/"):
            shutil.rmtree(self.output_dir + "/Sequence_database/")
        os.mkdir(self.output_dir + "/Sequence_database/")
        gene_bed = self.gene_fa.option("gene_bed").prop["path"]
        transcript_bed = self.gene_fa.option("transcript_bed").prop["path"]
        gene_fasta = self.gene_fa.option("gene_fa").prop["path"]
        if self.option("target_source") == "mongo":
            transcript_fasta = os.path.join(self.download.work_dir, "all_transcripts.fa")
            t2g = os.path.join(self.download.work_dir, "trans2gene.txt")
        else:
            transcript_fasta = self.transcripts.option("trans_fa").prop["path"]
            t2g = self.transcripts.option("trans2gene").prop["path"]
        species_urls = self.hyperlink
        biomart_file = self.des
        biomart_type = self.des_type
        cds_fasta = self.known_cds
        pep_fasta = self.known_pep
        gene_stat = os.path.join(self.sequence_detail.work_dir, 'gene_stat')
        self.api_gene_detail.add_gene_detail(db_path, gene_fasta, transcript_fasta, t2g, biomart_file, biomart_type,
                                             cds_fasta, pep_fasta, gene_bed, transcript_bed, species_urls, gene_stat)

    @time_count
    def build_seq_database(self):
        self.export_seq = self.api.api("small_rna.seq_detail")

        fasta = self.option("assembly_file").prop["path"]
        base_name = os.path.basename(fasta)
        cds = os.path.join(self.cds_predict.output_dir, '{}.transdecoder.cds'.format(base_name))
        pep = os.path.join(self.cds_predict.output_dir, '{}.transdecoder.pep'.format(base_name))
        trans2unigene = os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map")
        if os.path.exists(self.output_dir + "/Sequence_database/"):
            shutil.rmtree(self.output_dir + "/Sequence_database/")
        os.mkdir(self.output_dir + "/Sequence_database/")
        seq_db = self.output_dir + "/Sequence_database/refrna_seqs.db"
        self.export_seq.build_seq_database(seq_db, cds, pep, fasta, trans2unigene, task_id=self.task_id)

    # 添加家族分析信息
    @time_count
    def export_family(self):
        gevent.sleep()
        self.api_family = self.api.api('small_rna.smallrna_family_analyse')
        params = {
            "name": 'family_analyse',
            "database": self.option('taxonmy'),
        }
        family_path = self.mirna_family.output_dir
        self.api_family.run(family_path, params=params)

    # 添加碱基偏好性信息
    @time_count
    def export_bias(self):
        gevent.sleep()
        self.api_bias = self.api.api('small_rna.atcg_bias')
        params = {
            "name": 'atcg_analyse',
            "database": self.option('taxonmy'),
        }
        bias_path = self.mirna_atcg_bias.output_dir
        self.api_bias.run(bias_path, params=params)

    @time_count
    def export_mirna_edit(self):
        self.edit_api = self.api.api("small_rna.mirna_edit")
        self.edit_api.add_mirna_edit(self.mirna_edit.output_dir, task_id=self.task_id, project_sn=self.project_sn)
