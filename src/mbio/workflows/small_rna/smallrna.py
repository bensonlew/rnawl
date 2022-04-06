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
            ## 基础参数设置
            {"name": "fq_type", "type": "string", "default": "PE"},  # 测序类型，PE OR SE OR other，为other时，跳过质控步骤
            {"name": "read_length", "type": "int", "default": "150"},  # 测序读长，50,75,150
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # 测序质量，phred+64 OR phred+33
            {"name": "lib_type", "type": "string", "default": "other"},  # 测序类型，typeI OR type II OR type III or other
            {"name": "cut_5", "type": "int", "default": "4"},  # 当libtype为type I时，默认去除5'端4bp的碱基，后面仅需去除接头，,然后截取75bp的序列，最后再去除3'端4bp
            {"name": "extract_length", "type": "int", "default": "75"},  # 当libtype为type II时，默认去除5'端3bp的碱基，然后去除接头，最后截取75bp的序列
            {"name": "cut_3", "type": "int", "default": "4"},  # 当libtype为type III时，直接截取75bp的序列，最后再去除3'端15bp
            {"name": "adapter", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # adapter序列
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹，必须包含list.txt文件
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"}, # 对照表
            {"name": "analysis_strategy", "type": "string", "default": "type I"},  # 分析策略，分为参考基因组序列+完整GTF注释文件+功能注释结果(type I),仅参考基因组序列(type II),无参考基因组序列(type III)
            {"name": "taxonmy", "type":"string", "default": ""}, # 物种类别，Animal OR Plant
            {"name": "ref_genome", "type": "string", "default": ""},  # 参考基因组，具体物种名称
            {"name": "genome_version", "type": "string", "default": ""},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": ""},  # 参考基因组注释版本
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"},  # trinity组装结果文件
            {"name": "gene_to_trans", "type": "infile", 'format': "denovo_rna_v2.common"},  # 基因和转录本对应关系文件
            {"name": "mirbase_category", "type":"string", "default": "Animal"}, # mirbase物种类别，Animal OR Plant
            {"name": "mirbase_specie", "type": "string", "default": ""},  # mirbase参考物种名，多选时分号分隔
            {"name": "mismatch", "type": "int", "default": 1},  # 鉴定已知miRNA时允许的错配个数

            ## 高级参数设置
            # 靶基因预测
            {"name": "target_predict", "type": "bool", "default": False}, # 是否进行靶基因预测
            {"name": "target_seq", "type": "infile", 'format': "sequence.fasta"},  # 预测靶基因所用序列
            {"name": "target_animal_predict_method", "type": "string", "default": ""},  # 动物预测靶基因所用方法，多选时用分号分隔
            {"name": "target_plant_predict_method", "type": "string", "default": ""},  # 植物预测靶基因所用方法，多选时用分号分隔
            {"name": "min_support", "type": "int", "default": 1},  # 候选靶基因设定标准
            # 靶基因转录组功能注释
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
            # 表达差异分析
            {"name": "diff_method", "type": "string", "default": "DESeq2"},# 差异表达分析方法，DESeq2 or DEGseq or edgeR
            {"name": "diff_fdr_ci", "type": "float", "default": 0.05},  # 显著性水平
            {"name": "fc", "type": "float", "default": 2},
            {"name": "pvalue_padjust", "type": "string", "default": "padjust"},  #选择判断显著性水平的指标
            {"name": "diff_fdr_ci", "type": "string", "default": 0.05},  # 显著性水平
            {"name": "padjust_way", "type": "string", "default": "BH"},  #Bonferroni,Holm,BH,BY
            {"name": "is_miRNA_edit", "type": "bool", "default": True}, # 是否分析miRNA碱基编辑
            {"name": "version", "type": "string", "default": "v1.1"},
        ]
        #获取输出目录
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

        #添加tool/module
        self.filecheck = self.add_tool("small_rna.file_check")
        self.gunzip = self.add_tool("small_rna.gzfastq2fastq")
        self.qc = self.add_module("small_rna.mirna_qc")
        self.uniq = self.add_tool("small_rna.fasta_uniq")
        self.qc_stat_before = self.add_module("small_rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("small_rna.hiseq_reads_stat")
        self.mapping = self.add_tool("small_rna.mapper_and_stat")
        self.annotation = self.add_module("denovo_rna_v2.denovo_annotation")
        self.cds_predict = self.add_module("denovo_rna_v2.cds_tf")
        self.srna = self.add_module("small_rna.srna.srna")
        if self.option("target_predict") == True:
            self.mirna_target = self.add_module("small_rna.target_predict")
        self.exp_pca = self.add_tool("small_rna.exp_pca")
        self.exp_corr = self.add_tool("small_rna.exp_corr")
        self.diffexpress = self.add_tool("small_rna.diffexp")
        self.mirna_edit = self.add_module("small_rna.mirna_edit")
        self.mirna_family = self.add_tool("small_rna.smallrna_family_analyse")
        self.mirna_atcg_bias = self.add_tool("small_rna.atcg_bias")
        self.gene_fa = self.add_tool("ref_rna_v2.gene_fa")
        self.transcripts = self.add_tool("ref_rna_v2.transcript_abstract")
        self.extract_utr3 = self.add_tool("small_rna.extract_utr3")

        #判断流程结束tool/module list
        self.final_tools = [self.diffexpress, self.gene_fa, self.mirna_atcg_bias, self.mirna_family]
        if self.option("group_table").prop["sample_number"] > 2:
            self.final_tools.append(self.exp_pca)
        if self.option("target_predict") == True:
            self.final_tools.append(self.mirna_target)
            if self.option("assembly_file").is_set:
                self.final_tools.append(self.annotation)
        if self.option("is_miRNA_edit") == True:
            self.final_tools.append(self.mirna_edit)
        if not self.option("assembly_file").is_set:
            self.final_tools.append(self.transcripts)
        self.logger.info(self.final_tools)

        # 添加step，显示在页面进度条
        self.step.add_steps("filecheck", "rna_qc", "gunzip", "uniq","mapping", "srna", "diffexpress", "mirna_family",
                            "mirna_atcg_bias", "gene_fa", "annotation", "cds_predict", "exp_pca", "exp_corr", "mirna_target")

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        data = os.path.join(self.work_dir, 'data.json')
        if os.path.exists(data):
            with open(data, 'r') as load_f:
                load_dict = json.load(load_f)
                if 'rerun' in load_dict and load_dict['rerun']:
                    self.logger.info("该项目重运行中，先删除mongo库中已有数据")
                    self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
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
        if not self.option("fq_type") in ["PE", "SE", "other"]:
            raise OptionError("fq序列类型应为PE或SE或other")
        if not re.match(r'^[ATCGN]*$', self.option("adapter")):
            raise OptionError("接头序列输入错误")
        if not self.option("lib_type") in  ["type I", "type II", "type III", "other"]:
            raise OptionError("建库类型必须为type I或type II或type III或other")
        if self.option("lib_type") == "type I":
            if not (int(self.option("cut_5")) == 4 and int(self.option("cut_3")) == 4 and int(self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type I时，5'端和3'端的截取序列长度必须为4，序列提取长度必须为75")
        elif self.option("lib_type") == "type II":
            if not (int(self.option("cut_5")) == 3 and int(self.option("cut_3")) == 0 and int(self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type II时，5'端截取序列的长度必须为3, 3'端截取序列的长度必须为0，序列提取长度必须为75")
        elif self.option("lib_type") == "type III":
            if not (int(self.option("cut_5")) == 3 and int(self.option("cut_3")) == 15 and int(self.option("extract_length")) == 75):
                raise OptionError("当建库类型选择type III时，5'端截取序列的长度必须为0, 3'端截取序列的长度必须为15，序列提取长度必须为75")
        if not self.option("quality_score_system").lower() in ["phred+33", "phred 33", "phred+64", "phred 64"]:
            raise OptionError("测序质量参数输入有误")
        if not self.option("analysis_strategy") in ["type I", "type II"]:
            raise OptionError("分析策略参数输入有误")
        if not (self.option("ref_genome") and self.option("taxonmy") and self.option("genome_version") and self.option("genome_annot_version")):
            raise OptionError("参考基因组信息输入有误")
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        col = db["sg_genome_db"]
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        try:
            genome_info = col.find_one({"name" : self.option("ref_genome"), "assembly" : self.option("genome_version"), "annot_version" : self.option("genome_annot_version")})
            self.ref_annot_dir = os.path.join(db_path, genome_info["anno_path_v2"])
            self.annot_class = os.path.join(db_path, genome_info["anno_path_v2"], "annot_class")
            self.g2t2p = os.path.join(db_path, genome_info["anno_path_v2"], "annot_class", "tran2gene.txt")
            self.ref_genome = os.path.join(db_path, genome_info["dna_fa"])
            self.ref_index = os.path.join(db_path, genome_info["dna_index"])
            self.ref_gtf = os.path.join(db_path, genome_info["gtf"])
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
            self.anno_detail = os.path.join(db_path, genome_info["anno_path_v2"], "annot_class/anno_stat/all_anno_detail.xls")
            self.genome_version = genome_info.get('version', "")
        except:
            self.set_error("数据库中不存在该物种信息，程序退出")
        if os.path.exists(os.path.join(self.ref_annot_dir, "repeatmasker")):
            self.repeatmasker = True
        else:
            self.repeatmasker = False
        if self.option("analysis_strategy") == "type II" and self.option("target_predict") == True:
            if not self.option("assembly_file").is_set and not self.option("gene_to_trans").is_set:
                raise OptionError("选择该分析策略且需要预测靶基因时，必须上传组装结果文件、基因与转录本对应关系文件")
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
        mirna workflow run方法
        """
        self.filecheck.on('end', self.run_gunzip)
        self.filecheck.on('end', self.run_gene_fa)
        if not self.option("assembly_file").is_set:
            self.filecheck.on('end', self.run_transcripts)
        if self.option("assembly_file").is_set and self.option("target_predict") == True:
            self.filecheck.on('end', self.run_diamond)
            self.filecheck.on('end', self.run_cds_predict)
        self.gunzip.on('end', self.run_qc)
        self.qc.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc.on('end', self.run_uniq)
        self.uniq.on('end', self.run_mapping)
        self.mapping.on('end', self.run_srna)
        if self.option("is_miRNA_edit") == True:
            self.srna.on('end', self.run_mirna_edit)
        self.srna.on('end', self.run_mirna_family)
        self.srna.on('end', self.run_mirna_atcg_bias)
        self.srna.on('end', self.run_exp_corr)
        self.srna.on('end', self.run_diffexpress)
        if self.option("target_predict") == True:
            if self.option("target_seq").is_set:
                self.srna.on('end', self.run_mirna_target)
            else:
                if self.option("taxonmy").lower() == "animal":
                    self.filecheck.on('end', self.run_extract_utr3)
                    if self.option("assembly_file").is_set:
                        self.on_rely([self.extract_utr3, self.srna, self.annotation], self.run_mirna_target)
                    else:
                        self.on_rely([self.extract_utr3, self.srna], self.run_mirna_target)
                else:
                    self.srna.on('end', self.run_mirna_target)
        if self.option("group_table").prop["sample_number"] > 2:
            self.srna.on('end', self.run_exp_pca)
        self.on_rely(self.final_tools, self.end)
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

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
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
        self.blast_modules = []
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
        self.on_rely([self.diamond_nr, self.diamond_kegg, self.diamond_string, self.blast_swissprot, self.cds_predict], self.run_annotation)

    def run_cds_predict(self):
        self.logger.info("开始运行cds预测")
        opts = {
            "fasta" : self.option("assembly_file"),
            "e_value" : self.option('pfam_blast_evalue'),
            "species_type" : self.option("taxonmy"),
            "isoform_unigene" : os.path.join(self.filecheck.work_dir, "Trinity_t2g2u"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.set_output, "cds_predict")
        self.cds_predict.on('start', self.set_step, {'start': self.step.cds_predict})
        self.cds_predict.on('end', self.set_step, {'end': self.step.cds_predict})
        self.cds_predict.run()

    def run_annotation(self):
        self.logger.info("开始运行注释统计")
        anno_opts = {
            "gene2trans" : os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map"),
            "go_annot" : True,
            "nr_annot" : True,
            "taxonomy" : self.option("kegg_database"),
            "blast_nr_xml" : self.diamond_nr.option('outxml'),
            "blast_kegg_xml" : self.diamond_kegg.option('outxml'),
            "blast_string_xml" : self.diamond_string.option('outxml'),
            "blast_swissprot_xml" : self.blast_swissprot.option('outxml'),
            "pfam_domain" : self.cds_predict.output_dir + "/pfam_domain"
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
        if self.option('quality_score_system').lower() == "phred+33" or self.option('quality_score_system').lower() == "phred 33":
            fastq_format = "Q33"
        elif self.option('quality_score_system').lower() == "phred+64" or self.option('quality_score_system').lower() == "phred 64":
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
        if self.option('quality_score_system').lower() == "phred+33" or self.option('quality_score_system').lower() == "phred 33":
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
            "category": self.option("mirbase_category"),
            "species": self.option("mirbase_specie"),
            "reference": self.ref_genome,
            "clean_fa": os.path.join(self.uniq.output_dir, "uniq.fasta"),
            "config": self.qc.option("config_file"),
            "arf": os.path.join(self.mapping.output_dir, "reads_vs_genome.arf"),
            "repeat": self.repeatmasker,
            "gtf": self.ref_gtf,
            "list": self.option("fastq_dir").prop["path"] + "/list_re.txt",
            "qc_output": self.qc.output_dir,
            "mismatch": self.option("mismatch"),
        }
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
        specie = self.option("mirbase_specie").split(",")[0]
        mature_fa = self.srna.output_dir + "/known_mirna/mature.fa"
        hairpin_fa = self.config.SOFTWARE_DIR + "/database/mirbase/hairpin.fa"
        opts = {
            "list_file": os.path.join(self.qc.output_dir, "clean_data/list.txt"),
            "species": specie,
            "hairpin_fa": hairpin_fa,
            "mature_fa": mature_fa,
            "index": self.ref_index
        }
        self.mirna_edit.set_options(opts)
        self.mirna_edit.on("end", self.set_output, "mirna_edit")
        # self.mirna_edit.on("start", self.set_step, {"start": self.step.mirna_edit})
        # self.mirna_edit.on("end", self.set_step, {"end": self.step.mirna_edit})
        self.mirna_edit.run()

    def run_mirna_family(self):
        opts = {
            "mir": self.srna.output_dir + "/known_mirna_count.xls",
            "matfa": self.srna.output_dir + "/known_mirna/mature.fa",
            "novofa": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
        }
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
        if not self.option("assembly_file").is_set:
            if os.path.exists(self.ref_gtf):
                opts = {
                    "ref": self.ref_genome,
                    "gtf": self.ref_gtf
                }
            else:
                opts = {
                    "transcript": self.trans_fa,
                }
        else:
            opts = {
                "transcript": self.option("assembly_file").prop["path"],
            }
        self.extract_utr3.set_options(opts)
        self.extract_utr3.run()

    def run_mirna_target(self):
        if self.option("target_seq").is_set:
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
            if not self.option("assembly_file").is_set:
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
            else:
                if self.option("taxonmy").lower() == "plant":
                    opts = {
                        "known": self.srna.output_dir + "/known_mirna/mature.fa",
                        "novol": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
                        "ref": self.option("assembly_file").prop['path'],
                        "method": self.option("target_plant_predict_method"),
                        "type": self.option("taxonmy").lower(),
                        "species": self.option("ref_genome"),
                        "version": "v1.1",
                        "anno_detail": os.path.join(self.annotation.output_dir, "anno_stat/trans_anno_detail.xls"),
                        "min_support": int(self.option("min_support"))
                    }
                else:
                    opts = {
                        "known": self.srna.output_dir + "/known_mirna/mature.fa",
                        "novol": self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa",
                        "ref": os.path.join(self.extract_utr3.output_dir, "utr3.fa"),
                        "method": self.option("target_animal_predict_method"),
                        "type": self.option("taxonmy").lower(),
                        "species": self.option("ref_genome"),
                        "version": "v1.1",
                        "anno_detail": os.path.join(self.annotation.output_dir, "anno_stat/trans_anno_detail.xls"),
                        "min_support": int(self.option("min_support"))
                    }
        self.mirna_target.set_options(opts)
        self.mirna_target.on('end', self.set_output, 'mirna_target')
        self.mirna_target.on('start', self.set_step, {'start': self.step.mirna_target})
        self.mirna_target.on('end', self.set_step, {'end': self.step.mirna_target})
        self.mirna_target.run()

    def run_gene_fa(self):
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
            "express_matrix" : self.srna.option("total_mirna_norm").prop["path"]
        }
        self.exp_pca.set_options(opts)
        self.exp_pca.on("end", self.set_output, "exp_pca")
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.run()

    def run_exp_corr(self):
        self.logger.info("开始运行聚类分析")
        opts = {
            "express_matrix" : self.srna.option("total_mirna_norm").prop["path"]
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
            "count" : count_file,
            "exp" : tpm_file,
            "group" : self.option("group_table").prop["path"],
            "cmp" : self.option("control_file").prop["path"],
            "pvalue_padjust" : self.option("pvalue_padjust"),
            "pvalue" : float(self.option("diff_fdr_ci")),
            "fc" : float(self.option("fc")),
            "padjust_way" : self.option("padjust_way"),
            "method" : self.option("diff_method"),
        }
        self.diffexpress.set_options(opts)
        self.diffexpress.on("end", self.set_output, "diffexpress")
        self.diffexpress.on('start', self.set_step, {'start': self.step.diffexpress})
        self.diffexpress.on('end', self.set_step, {'end': self.step.diffexpress})
        self.diffexpress.run()

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
        if self.option("target_predict") == True:
            col.update({"task_id" : self.task_id}, {"$set": {"refrna_seqdb": self.workflow_output + "/Sequence_database/refrna_seqs.db"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"ref_gtf": self.ref_gtf}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"ref_genome": self.ref_genome}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"genome_id": self.genome_id}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"organism_name": self.option("ref_genome")}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"annot_version": self.option("genome_annot_version")}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"assembly": self.option("genome_version")}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"known_pre": self.workflow_output + "/sRNA/known_miRNA/known_miRNA_detail.xls"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"novol_pre": self.workflow_output + "/sRNA/novel_miRNA/novel_miRNA_detail.xls"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"novol": self.workflow_output + "/sRNA/novel_miRNA/novel_mature.fa"}}, upsert=True)
        col.update({"task_id" : self.task_id}, {"$set": {"known": self.workflow_output + "/sRNA/known_miRNA/known_mature.fa"}}, upsert=True)

        col.update({"task_id" : self.task_id}, {"$set": {"version": self.option("version")}}, upsert=True)
        if self.genome_version == "v2":
            col.update({'task_id': self.task_id},
                       {'$set': {'database_version': {"kegg": "202003"}}}, upsert=True)

        ## 用于靶基因交互分析
        if self.option("assembly_file").is_set:
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.workflow_output + "/Background/upload.assembly.fa"}}, upsert=True)
            col5 = db["sg_annotation_stat"]
            col5.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Annotation/"}}, upsert=True)
        else:
            col.update({"task_id" : self.task_id}, {"$set": {"assemble_fa": self.trans_fa}}, upsert=True)
        if self.option("gene_to_trans").is_set:
            col.update({"task_id" : self.task_id}, {"$set": {"gene_to_trans": self.workflow_output + "/Background/upload.gene2trans.txt"}}, upsert=True)
        else:
            col.update({"task_id" : self.task_id}, {"$set": {"gene_to_trans": self.g2t2p}}, upsert=True)
        ## 用于页面展示pdf图片
        col1 = db["sg_novel_mirna"]
        col1.update({"task_id" : self.task_id}, {"$set": {"pdf_dir": self.workflow_output + "/sRNA/novel_miRNA/novel_pre_structure_pdf"}}, upsert=True)
        col2 = db["sg_known_mirna"]
        col2.update({"task_id" : self.task_id}, {"$set": {"pdf_dir": self.workflow_output + "/sRNA/known_miRNA/known_pre_structure_pdf"}}, upsert=True)
        if self.option("target_predict") == True:
            col3 = db["sg_target"]
            col3.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/miRNA_Target"}}, upsert=True)
        col4 = db["sg_mapping"]
        col4.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Mapping"}}, upsert=True)
        col4 = db["sg_circos"]
        col4.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Mapping"}}, upsert=True)
        self.modify_output()

        super(SmallrnaWorkflow, self).end()

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
        if self.option("assembly_file").is_set:
            os.link(self.option("assembly_file").prop["path"], target_dir + "/Background/upload.assembly.fa")
        if self.option("gene_to_trans").is_set:
            os.link(self.option("gene_to_trans").prop["path"], target_dir + "/Background/upload.gene2trans.txt")
        #seq_db
        if self.option("target_predict") == True:
            os.mkdir(target_dir + "/Sequence_database")
            seq_db = origin_dir + "/Sequence_database/refrna_seqs.db"
            os.link(seq_db, target_dir + "/Sequence_database/refrna_seqs.db")
        # QC
        fq_stat_before = origin_dir + "/QC/stat_rawdata/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC/stat_cleandata/fastq_stat.xls"
        os.mkdir(target_dir + "/QC")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        os.mkdir(target_dir + "/QC/rawdata")
        for file in os.listdir(os.path.join(self.qc.output_dir, "raw_data")):
            file_path = os.path.join(self.qc.output_dir, "raw_data", file)
            if os.path.basename(file_path) != "list.txt":
                os.link(file_path, target_dir + "/QC/rawdata/" + file)
        # Uniq
        uniq_table = origin_dir + "/Uniq/table.xls"
        uniq_fa = origin_dir + "/Uniq/unique.fasta"
        os.mkdir(target_dir + "/Uniq")
        os.link(uniq_table, target_dir + "/Uniq/uniq.xls")
        os.link(uniq_fa, target_dir + "/Uniq/uniq.fa")
        # Mapping
        os.makedirs(target_dir + "/Mapping")
        for file in os.listdir(origin_dir + "/Mapping"):
            file_path = os.path.join(origin_dir, "Mapping", file)
            os.link(file_path, target_dir + "/Mapping/" + file)
        # Srna
        os.makedirs(target_dir + "/sRNA")
        ## KnownMirna
        os.makedirs(target_dir + "/sRNA/known_miRNA")
        known_mirna_detail = os.path.join(origin_dir, "Srna/known_mirna/known_mirna_detail.xls")
        os.link(known_mirna_detail, target_dir + "/sRNA/known_miRNA/known_miRNA_detail.xls")
        mature_fa = os.path.join(origin_dir, "Srna/known_mirna/mature.fa")
        os.link(mature_fa, target_dir + "/sRNA/known_miRNA/known_mature.fa")
        hairpin_fa = os.path.join(origin_dir, "Srna/known_mirna/hairpin.fa")
        os.link(hairpin_fa, target_dir + "/sRNA/known_miRNA/known_hairpin.fa")
        os.makedirs(target_dir + "/sRNA/known_miRNA/known_pre_structure_pdf")
        for file in os.listdir(origin_dir + "/Srna/known_mirna/structure_pdf"):
            file_path = os.path.join(origin_dir, "Srna/known_mirna/structure_pdf", file)
            os.link(file_path, target_dir + "/sRNA/known_miRNA/known_pre_structure_pdf/" + file)
        ## NovelMirna
        os.makedirs(target_dir + "/sRNA/novel_miRNA")
        novel_mirna_detail = os.path.join(origin_dir, "Srna/novel_mirna/novel_mirna_detail.xls")
        os.link(novel_mirna_detail, target_dir + "/sRNA/novel_miRNA/novel_miRNA_detail.xls")
        novel_mature_seq = os.path.join(origin_dir, "Srna/novel_mirna/novel_mature_seq.fa")
        os.link(novel_mature_seq, target_dir + "/sRNA/novel_miRNA/novel_mature.fa")
        novel_precursor_seq = os.path.join(origin_dir, "Srna/novel_mirna/novel_precursor_seq.fa")
        os.link(novel_precursor_seq, target_dir + "/sRNA/novel_miRNA/novel_hairpin.fa")
        os.makedirs(target_dir + "/sRNA/novel_miRNA/novel_pre_structure_pdf")
        for file in os.listdir(origin_dir + "/Srna/novel_mirna/structure_pdf"):
            file_path = os.path.join(origin_dir, "Srna/novel_mirna/structure_pdf", file)
            os.link(file_path, target_dir + "/sRNA/novel_miRNA/novel_pre_structure_pdf/" + file)
        ## SrnaStat
        os.makedirs(target_dir + "/sRNA/sRNA_stat")
        mirna_stat = os.path.join(origin_dir, "Srna/srna_stat/mirna_stat.xls")
        os.link(mirna_stat, target_dir + "/sRNA/sRNA_stat/miRNA_stat.xls")
        ncrna_stat = os.path.join(origin_dir, "Srna/srna_stat/ncrna_stat.xls")
        os.link(ncrna_stat, target_dir + "/sRNA/sRNA_stat/ncRNA_stat.xls")
        srna_stat = os.path.join(origin_dir, "Srna/srna_stat/srna_stat.xls")
        os.link(srna_stat, target_dir + "/sRNA/sRNA_stat/sRNA_stat.xls")
        # Express
        os.makedirs(target_dir + "/Express")
        os.makedirs(target_dir + "/Express/miRNA_express")
        total_mirna_count = os.path.join(origin_dir, "Srna/total_mirna_count.xls")
        os.link(total_mirna_count, target_dir + "/Express/miRNA_express/All_miRNA_count.xls")
        total_mirna_norm = os.path.join(origin_dir, "Srna/total_mirna_norm.xls")
        os.link(total_mirna_norm, target_dir + "/Express/miRNA_express/All_miRNA_tpm.xls")
        known_mirna_count = os.path.join(origin_dir, "Srna/known_mirna/known_mirna_count.xls")
        os.link(known_mirna_count, target_dir + "/Express/miRNA_express/known_miRNA_count.xls")
        novel_mirna_count = os.path.join(origin_dir, "Srna/novel_mirna/novel_mirna_count.xls")
        os.link(novel_mirna_count, target_dir + "/Express/miRNA_express/novel_miRNA_count.xls")
        known_mirna_norm = os.path.join(origin_dir, "Srna/known_mirna/known_mirna_norm.xls")
        os.link(known_mirna_norm, target_dir + "/Express/miRNA_express/known_miRNA_tmp.xls")
        novel_mirna_norm = os.path.join(origin_dir, "Srna/novel_mirna/novel_mirna_norm.xls")
        os.link(novel_mirna_norm, target_dir + "/Express/miRNA_express/novel_miRNA_tmp.xls")
        # ExpCorr
        os.makedirs(target_dir + "/Express/ExpCorr")
        file_path = os.path.join(origin_dir, "ExpCorr/sample_correlation.xls")
        os.link(file_path, target_dir + "/Express/ExpCorr/sample_correlation.xls")
        # ExpPca
        if self.option("group_table").prop["sample_number"] > 2:
            os.makedirs(target_dir + "/Express/ExpPca")
            Explained_variance_ratio = os.path.join(origin_dir, "ExpPca/Explained_variance_ratio.xls")
            os.link(Explained_variance_ratio, target_dir + "/Express/ExpPca/Explained_variance_ratio.xls")
            PCA = os.path.join(origin_dir, "ExpPca/PCA.xls")
            os.link(PCA, target_dir + "/Express/ExpPca/PCA.xls")
        # DiffExpress
        os.makedirs(target_dir + "/DiffExpress/")
        for file in os.listdir(origin_dir + "/DiffExpress"):
            file_path = os.path.join(origin_dir, "DiffExpress", file)
            os.link(file_path, os.path.join(target_dir, "DiffExpress", file))
        files = glob.glob(os.path.join(target_dir, "DiffExpress/*.DE.list"))
        for file in files:
            os.remove(file)
        # miRNA_Structure
        os.makedirs(target_dir + "/miRNA_Structure")
        ## MirnaBias
        os.makedirs(target_dir + "/miRNA_Structure/miRNA_bias")
        for file in os.listdir(origin_dir + "/MirnaBias"):
            file_path = os.path.join(origin_dir, "MirnaBias", file)
            os.link(file_path, os.path.join(target_dir, "miRNA_Structure/miRNA_bias", file))
        ## MirnaEdit
        if self.option("is_miRNA_edit") == True:
            os.makedirs(target_dir + "/miRNA_Structure/miRNA_edit")
            for file in os.listdir(self.mirna_edit.output_dir):
                file_path = os.path.join(self.mirna_edit.output_dir, file)
                if os.path.basename(file_path) == "all.result.xls":
                    os.link(file_path, os.path.join(target_dir, "miRNA_Structure/miRNA_edit/edit_stat.xls",))
                else:
                    os.link(file_path, os.path.join(target_dir, "miRNA_Structure/miRNA_edit", file.split(".")[0] + ".edit_detail.xls"))
        ## MirnaFamily
        os.makedirs(target_dir + "/miRNA_Structure/miRNA_family")
        for file in os.listdir(origin_dir + "/MirnaFamily"):
            file_path = os.path.join(origin_dir, "MirnaFamily", file)
            os.link(file_path, os.path.join(target_dir, "miRNA_Structure/miRNA_family", file))
        os.rename(os.path.join(target_dir, "miRNA_Structure/miRNA_family/family.species.xls"), os.path.join(target_dir, "miRNA_Structure/miRNA_family/family_species.xls"))
        # MirnaTarget
        if self.option("target_predict") == True:
            os.makedirs(target_dir + "/miRNA_Target")
            for file in os.listdir(origin_dir + "/MirnaTarget"):
                file_path = os.path.join(origin_dir, "MirnaTarget", file)
                os.link(file_path, os.path.join(target_dir, "miRNA_Target", file))
        if self.option("assembly_file").is_set and self.option("target_predict") == True:
            os.makedirs(target_dir + "/Annotation")
            CopyFile().linkdir(origin_dir + "/Annotation", os.path.join(target_dir, "Annotation"))
            if os.path.exists(os.path.join(target_dir, "Annotation/tran2gene.txt")):
                os.remove(os.path.join(target_dir, "Annotation/tran2gene.txt"))
            os.link(os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map"), os.path.join(target_dir, "Annotation/tran2gene.txt"))


        sdir = self.add_upload_dir(target_dir)
        sdir.add_regexp_rules([
            [r"Background/.*\.gtf\.genome_stat\.xls", "xls", "参考基因组注释信息表", 0],
            [r"Mapping/.*_map_stat\.xls", "xls", "与基因组比对结果表", 0],
            [r"Mapping/.*_circos\.svg", "", "不同染色体reads分布circos图", 0],
            [r"miRNA_Structure/miRNA_edit/.*\.edit_detail\.xls", "xls", "miRNA碱基编辑分析详情表（单个样本）", 0],
            [r"DiffExpress/.*_vs_.*\.xls", "xls", "差异分析详情表（单个比较组）", 0],
            [r"DiffExpress/.*_diff_summary\.xls", "xls", "表达量差异统计表（所有比较组）", 0],
            [r"DiffExpress/.*_vs_.*\..*normalize\.xls", "xls", "差异矩阵标准化结果表", 0],
            [r"DiffExpress/.*_vs_.*\..*sizeFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
            [r"DiffExpress/.*_vs_.*\..*normFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
        ])
        sdir.add_relpath_rules([
            [".", "", "流程分析结果目录", 0],
            ["Background", "", "项目背景目录", 0],
            ["Sequence_database", "", "序列文件数据库", 1],
            ["Sequence_database/refrna_seqs.db", "", "序列文件", 1],
            ["QC", "", "测序数据质控结果目录", 0],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表", 0],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表", 0],
            ["QC/rawdata", "", "原始测序数据", 0],
            ["Uniq", "", "Uniq序列文件目录", 0],
            ["Uniq/uniq.xls", "", "Uniq序列在各样本的分布详情表", 0],
            ["Uniq/uniq.fa", "", "Uniq序列", 0],
            ["sRNA", "", "sRNA分析结果目录", 0],
            ["sRNA/known_miRNA", "", "已知miRNA鉴定结果文件", 0],
            ["sRNA/known_miRNA/known_miRNA_detail.xls", "", "已知miRNA鉴定详情表", 0],
            ["sRNA/known_miRNA/known_mature.fa", "", "已知miRNA成熟体序列", 0],
            ["sRNA/known_miRNA/known_hairpin.fa", "", "已知miRNA前体序列", 0],
            ["sRNA/known_miRNA/known_pre_structure_pdf", "", "已知miRNA前体结构", 0],
            ["sRNA/novel_miRNA", "", "新miRNA预测结果文件", 0],
            ["sRNA/novel_miRNA/novel_miRNA_detail.xls", "", "新miRNA鉴定详情表", 0],
            ["sRNA/novel_miRNA/novel_mature.fa", "", "新miRNA成熟体序列", 0],
            ["sRNA/novel_miRNA/novel_hairpin.fa", "", "新miRNA前体序列", 0],
            ["sRNA/novel_miRNA/novel_pre_structure_pdf", "", "新miRNA前体结构", 0],
            ["sRNA/sRNA_stat", "", "sRNA统计结果文件", 0],
            ["sRNA/sRNA_stat/miRNA_stat.xls", "", "miRNA统计表", 0],
            ["sRNA/sRNA_stat/ncRNA_stat.xls", "", "ncRNA统计表", 0],
            ["sRNA/sRNA_stat/sRNA_stat.xls", "", "sRNA统计表", 0],
            ["miRNA_Structure", "", "miRNA结构分析结果目录", 0],
            ["miRNA_Structure/miRNA_bias", "", "miRNA碱基偏向性结果文件", 0],
            ["miRNA_Structure/miRNA_bias/all_first_bias.xls", "", "已知+新miRNA对应的首位碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/all_first_bias_per.xls", "", "已知+新miRNA对应的首位碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_bias/all_loc_bias.xls", "", "已知+新miRNA对应的不同长度碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/all_loc_bias_per.xls", "", "已知+新miRNA对应的不同长度碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_bias/known_first_bias.xls", "", "已知miRNA对应的首位碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/known_first_bias_per.xls", "", "已知miRNA对应的首位碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_bias/known_loc_bias.xls", "", "已知对应的不同长度碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/known_loc_bias_per.xls", "", "已知对应的不同长度碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_bias/novel_first_bias.xls", "", "已知对应的不同长度碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/novel_first_bias_per.xls", "", "已知对应的不同长度碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_bias/novel_loc_bias.xls", "", "新miRNA对应的不同长度碱基偏向性统计表", 0],
            ["miRNA_Structure/miRNA_bias/novel_loc_bias_per.xls", "", "新miRNA对应的不同长度碱基偏向性统计百分比表", 0],
            ["miRNA_Structure/miRNA_edit", "", "miRNA碱基编辑结果文件", 0],
            ["miRNA_Structure/miRNA_edit/edit_stat.xls", "", "miRNA碱基编辑分析统计表（所有样本）", 0],
            ["miRNA_Structure/miRNA_family", "", "miRNA家族分析结果文件", 0],
            ["miRNA_Structure/miRNA_family/known_miR_family.xls", "", "已知miRNA对应家族归属信息表", 0],
            ["miRNA_Structure/miRNA_family/novel_miR_family.xls", "", "新miRNA对应家族归属信息表", 0],
            ["miRNA_Structure/miRNA_family/family_species.xls", "", "新miRNA对应家族归属信息表", 0],
            ["miRNA_Target", "", "靶基因预测结果目录", 0],
            ["miRNA_Target/All_annot_target.xls", "", "已知+新miRNA对应的靶基因详情表", 0],
            ["miRNA_Target/known_target.xls", "", "已知miRNA对应的靶基因详情表", 0],
            ["miRNA_Target/novol_target.xls", "", "新miRNA对应的靶基因详情表", 0],
            ["miRNA_Target/target.fa", "", "靶基因序列信息", 0],
            [r"miRNA_Target/*detail.txt.gz", "", "靶基因比对详细信息", 0],
            ["Mapping", "", "与基因组比对结果目录", 0],
            ["Mapping/Genome_map_stat.xls", "", "与基因组比对结果表", 0],
            ["Mapping/Chro_map_stat.xls", "", "不同染色体reads分布统计表", 0],
            ["Mapping/reads_vs_genome.arf", "", "与参考基因组比对结果文件", 0],
            ["Mapping/uniq_mapped.fasta", "", "比对到基因组的序列文件", 0],
            ["Express", "", "表达量分析结果目录", 0],
            ["Express/miRNA_express", "", "表达定量结果文件", 0],
            ["Express/miRNA_express/All_miRNA_count.xls", "", "已知+新miRNA表达量count矩阵表", 0],
            ["Express/miRNA_express/All_miRNA_tpm.xls", "", "已知+新miRNA表达量tpm矩阵表", 0],
            ["Express/miRNA_express/known_miRNA_count.xls", "", "已知miRNA表达量count矩阵表", 0],
            ["Express/miRNA_express/known_miRNA_tpm.xls", "", "已知miRNA表达量tmp矩阵表", 0],
            ["Express/miRNA_express/novel_miRNA_count.xls", "", "新miRNA表达量count矩阵表", 0],
            ["Express/miRNA_express/novel_miRNA_tmp.xls", "", "新miRNA表达量tmp矩阵表", 0],
            ["Express/ExpCorr", "", "样本间相关性分析矩阵表", 0],
            ["Express/ExpCorr/sample_correlation.xls", "", "样本间相关性分析矩阵表", 0],
            ["Express/ExpPca", "", "样本间PCA结果文件", 0],
            ["Express/ExpPca/PCA.xls", "", "样本间PCA分析结果表", 0],
            ["Express/ExpPca/Explained_variance_ratio.xls", "", "样本间PCA主成分解释表", 0],
            ["DiffExpress", "", "表达量差异分析结果目录", 0],
            ["Annotation", "", "靶基因注释结果目录", 0],
            ])

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
        if self.option("target_predict") == True:
            self.export_target()
        if self.option("is_miRNA_edit") == True:
            self.export_mirna_edit()
        self.export_expression()
        if not self.option("assembly_file").is_set and self.option("target_predict") == True:
            self.export_gene_detail()
        else:
            if self.option("target_predict") == True:
                self.build_seq_database()
        self.export_family()
        self.export_bias()
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
        self.api_mapping.add_map_stat(sample_list=None, map_dir=map_dir, sample_file=sample_file, params=params)

    @time_count
    def export_genome_info(self):
        self.api_geno = self.api.api("small_rna.genome_info")
        file_path = self.genome_stat
        species_name = self.species_name
        species = self.species
        ref_anno_version = self.ref_anno_version
        hyperlink = self.hyperlink
        self.api_geno.add_genome_info(file_path=file_path,species_name=species_name, species=species, ref_anno_version=ref_anno_version,hyperlink=hyperlink)

    @time_count
    def export_qc(self):
        self.api_qc = self.api.api("small_rna.small_rna_qc")
        fq_type = "SE"
        quality_stat_after = self.qc_stat_after.output_dir
        quality_stat_before = self.qc_stat_before.output_dir
        self.api_qc.add_before_qc(quality_stat_before, fq_type=fq_type, about_qc="before")
        qc_stat = self.qc.output_dir + "/clean_data"
        self.api_qc.add_after_qc(quality_stat_after, qc_dir=qc_stat, fq_type=fq_type, about_qc="after")
        self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(self.option("group_table").prop["path"])
        self.logger.info("group_detail为：" + str(self.group_detail))
        self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control_file").prop["path"], self.group_id)
        self.compare_detail = compare_detail

    @time_count
    def export_srna(self):
        self.srna_api = self.api.api("small_rna.srna")
        known_mirna = os.path.join(self.srna.output_dir, "known_mirna/known_mirna_detail.xls")
        novel_mirna = os.path.join(self.srna.output_dir, "novel_mirna/novel_mirna_detail.xls")
        pdfs_known = os.path.join(self.srna.output_dir, "known_mirna/structure_pdf")
        pdfs_novel = os.path.join(self.srna.output_dir, "novel_mirna/structure_pdf")
        ncrna_stat = os.path.join(self.srna.output_dir, "srna_stat/ncrna_stat.xls")
        srna_stat = os.path.join(self.srna.output_dir, "srna_stat/srna_stat.xls")
        mirna_stat = os.path.join(self.srna.output_dir, "srna_stat/mirna_stat.xls")
        srna_stat_for_graph = os.path.join(self.srna.output_dir, "srna_stat/srna_stat_for_graph.xls")
        category = self.option("taxonmy").lower()
        if category == "plant":
            method = "mireap"
        else:
            method = "mirdeep2"
        self.srna_api.add_known_mirna(known_mirna, project_sn=self.project_sn, task_id=self.task_id, params={"method": "quantifier"}, pdfs=pdfs_known)
        self.srna_api.add_ncrna_stat(ncrna_stat, project_sn=self.project_sn, task_id=self.task_id, params={"method": "statistics"})
        self.srna_api.add_novel_mirna(novel_mirna, project_sn=self.project_sn, task_id=self.task_id, params={"method": method}, category=category, pdfs=pdfs_novel)
        self.srna_api.add_mirna_stat(mirna_stat, project_sn=self.project_sn, task_id=self.task_id, params={"method": method})
        self.srna_api.add_srna_stat(srna_stat, srna_stat_for_graph, project_sn=self.project_sn, task_id=self.task_id, params={"method": method})

    @time_count
    def export_target(self):
        self.target = self.api.api("small_rna.target_annotation")
        if not self.option("assembly_file").is_set:
            result_dir = self.annot_class
            g2t2p = self.g2t2p
        else:
            result_dir = self.annotation.output_dir
            g2t2p = os.path.join(self.filecheck.work_dir, "Trinity.gene_trans_map")
        target = self.mirna_target.output_dir
        species_name = self.species_name
        new_target_file = target + '/novol_target.xls'
        known_target_file = target + '/known_target.xls'
        known_seq = self.srna.output_dir + "/known_mirna/mature.fa"
        new_seq = self.srna.output_dir + "/novel_mirna/novel_mature_seq.fa"
        params = {
            "nr_evalue": 1e-5,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue":1e-5,
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
            params_target = {
                "miranda": miranda,
                "targetscan": targetscan,
                "rnahybrid": rnahybrid,
                "min_support": str(self.option("min_support")),
                'miranda_score': "160.0",
                'miranda_energy': "-20",
                'miranda_strict': "on",
                'rnahybird_num': "100",
                'rnahybird_energy': "-20",
                'rnahybird_pvalue': "0.01",
                'ps_robot_score': "2.5",
                'targetfinder_score': "4"
            }
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
                "min_support": str(self.option("min_support")),
                'miranda_score': "160.0",
                'miranda_energy': "-20",
                'miranda_strict': "on",
                'rnahybird_num': "100",
                'rnahybird_energy': "-20",
                'rnahybird_pvalue': "0.01",
                'ps_robot_score': "2.5",
                'targetfinder_score': "4"

            }

        self.target.run(known_target_file, new_target_file, result_dir, g2t2p, params, taxon=self.option("taxonmy"), exp_level='transcript', version=self.option("version"))
        self.target.import_target_detail(new_target_file, known_target_file, params_target, new_seq=new_seq, known_seq=known_seq, anno_type="origin", species_name= species_name, target_dir=target, version=self.option("version"))
        diff_summary = glob.glob(os.path.join(self.diffexpress.output_dir, "*_diff_summary.xls"))[0]
        self.target.add_target_geneset(new_target_file, known_target_file, diff_summary)

    @time_count
    def export_expression(self):
        known_count_xls = os.path.join(self.srna.output_dir, "known_mirna_count.xls")
        novel_count_xls = os.path.join(self.srna.output_dir, "novel_mirna_count.xls")
        known_tpm_xls = os.path.join(self.srna.output_dir, "known_mirna_norm.xls")
        novel_tpm_xls = os.path.join(self.srna.output_dir, "novel_mirna_norm.xls")

        # set basic variables
        gevent.sleep()
        all_exp = self.api.api("small_rna.all_exp")
        project_sn = self.project_sn
        task_id = self.task_id
        group_id = self.group_id
        group_dict = self.option('group_table').prop['group_dict']
        control_id = str(self.control_id)

        # create main_table in sg_exp and add detail_table in sg_exp_detail for count
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

        # create main_table in sg_exp and add detail_table in sg_exp_detail for tpm
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

        # create main_table in sg_exp_graph and add detail_table in sg_exp_graph_*
        all_tpm_xls = os.path.join(self.work_dir, 'all_mirna_norm.xls')
        all_tpm_pd = pd.concat([
            all_exp.process_exp_matrix(known_tpm_xls),
            all_exp.process_exp_matrix(novel_tpm_xls)
        ], axis=0)
        all_tpm_pd.to_csv(all_tpm_xls, sep='\t', header=True, index=True)
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
            all_exp.add_exp_pca(
                pca_output_dir=pca_output_dir,
                params=params,
                task_id=task_id,
                project_sn=project_sn
            )

        # create main_table in sg_diff and add detail_table in sg_diff_*
        diff_output_dir = self.diffexpress.output_dir
        params = {
            'task_id': task_id,
            'submit_location': 'diff_detail',
            'task_type': 2,
            'group_id': str(group_id),
            'group_dict': group_dict,
            'control_id': str(control_id),
            'exp_id': str(tpm_exp_id),
            'diff_method': self.option("diff_method"),
            'stat_type': self.option("pvalue_padjust"),
            'stat_cutoff': str(self.option("diff_fdr_ci")),
            'fc': str(self.option("fc")),
            'correct_method': self.option("padjust_way"),
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        all_exp.add_diffexp(
            diff_output=diff_output_dir,
            exp_id=str(tpm_exp_id),
            project_sn=project_sn,
            task_id=task_id,
            params=params,
            create_geneset=True,
            diff_method=self.option("diff_method")
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
        transcript_fasta = self.transcripts.option("trans_fa").prop["path"]
        t2g = self.transcripts.option("trans2gene").prop["path"]
        gene_fasta = self.gene_fa.option("gene_fa").prop["path"]
        species_urls = self.hyperlink
        biomart_file = self.des
        biomart_type = self.des_type
        cds_fasta = self.known_cds
        pep_fasta = self.known_pep
        self.api_gene_detail.add_gene_detail(db_path, gene_fasta, transcript_fasta, t2g, biomart_file, biomart_type,
                        cds_fasta, pep_fasta, gene_bed, transcript_bed, species_urls)

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
            "database": self.option('mirbase_category'),
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
            "database": self.option('mirbase_category'),
        }
        bias_path = self.mirna_atcg_bias.output_dir
        self.api_bias.run(bias_path, params=params)

    @time_count
    def export_mirna_edit(self):
        self.edit_api = self.api.api("small_rna.mirna_edit")
        self.edit_api.add_mirna_edit(self.mirna_edit.output_dir, task_id=self.task_id, project_sn=self.project_sn)
