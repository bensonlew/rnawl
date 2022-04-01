# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""多样性基础分析"""
from biocluster.workflow import Workflow
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
import os
import json
import shutil
from mainapp.models.mongo.submit.sequence.sample_extract import SampleExtract as SE
from mbio.packages.meta.common_function import group_file_spilt
from mbio.packages.meta.search_polution_by_list import check_pollution_pip
from mbio.packages.meta.delete_mongo import DeleteDemoMongo
import pandas as pd
import re
import glob
from functools import wraps

def tryforgood(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class MetaBaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        """
        self._sheet = wsheet_object
        super(MetaBaseWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq,sequence.fastq_dir'},  # 输入的fastq文件或fastq文件夹
            {'name': 'otu_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出的合并到一起的fasta，供后续的otu分析用
            {'name': 'identity', 'type': 'float', 'default': 0.97},  # 相似性值，范围0-1.
            {'name': 'otu_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出结果otu表
            {'name': 'otu_rep', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出结果otu代表序列
            # {'name': 'otu_seqids', 'type': 'outfile', 'format': 'meta.otu.otu_seqids'},  # 输出结果otu中包含序列列表
            {'name': 'otu_biom', 'type': 'outfile', 'format': 'meta.otu.biom'},  # 输出结果biom格式otu表
            {'name': 'revcomp', 'type': 'bool', 'default': False},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            # {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_dir'},  # 输出的otu_taxon_dir文件夹
            {"name": "estimate_indices", "type": "string", "default": "sobs,ace,chao,shannon,simpson,coverage"},
            {"name": "rarefy_indices", "type": "string", "default": "sobs,ace,chao,shannon,simpson,coverage"},  # 指数类型
            {"name": "rarefy_freq", "type": "int", "default": 100},
            {"name": "alpha_level", "type": "string", "default": "otu"},  # level水平
            {"name": "beta_analysis", "type": "string", "default": "anosim,pca,pcoa,nmds,hcluster"},
            {"name": "beta_level", "type": "string", "default": "otu"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            # {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "anosim_grouplab", "type": 'string', "default": ''},
            {"name": "plsda_grouplab", "type": 'string', "default": ''},
            {"name": "file_list", "type": "string", "default": "null"},
            {"name": "info_path", "type": "string", "default": ""},  # add by zhujuan 3018.05.03增加改参数，解决检测样品在流程中重复运行的问题
            {"name": "raw_sequence", "type": "infile", "format": "sequence.raw_sequence_txt"},
            {"name": "workdir_sample", "type": "string", "default": ""},
            {"name": "if_fungene", "type": "bool", 'default': False},
            {"name": "query_id", "type": "string","default": ""},    #add by liulinmeng @2018062
            {"name": "ref_acid", "type": "infile","format": "sequence.fasta"}, # 参考氨基酸fasta文件
            {"name": "pipeline", "type": "string", "default": ""}, # 选择常规流程还是功能基因流程
            {"name": "function_gene", "type": "string"}, # 功能基因流程功能基因名称
            {"name": "acid_length", "type": "int", "default": 80}, # 氨基酸长度阈值
            {"name": "seq_identity", "type": "float", "default": 0.4},
            {"name": "seq_style", "type": "string", "default": "nucl"},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "tax_database", "type": "string"}, # 物种分类数据库,
            {'name': 'tax_confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            #{"name": "tax_confidence", "type": "float", "default": 0.7}, # 分类置信度
            {'name': 'save_pdf', 'type': 'int', 'default': 1},
            {'name': 'otu_sub', "type": "string", 'default': "True"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.sample_rename = self.add_tool("meta.sample_rename")
        self.sample_rename = self.add_tool("metaasv.sample_rename")
        self.new_sample_extract = self.add_module("meta.sample_extract.sample_extract")  # added by shijin
        # self.sample_check = self.add_tool("meta.sample_check")
        self.sample_check = self.add_tool("meta.sample_merge")
        # self.info_abstract = self.add_tool("meta.lala.info_abstract")
        self.framebot = self.add_tool("meta.framebot")
        self.otu = self.add_tool("meta.otu.usearch_otu")
        self.phylo = self.add_tool("phylo.phylo_tree")
        self.tax = self.add_module("annotation.meta_tax")
        self.stat_depth = self.add_module("meta.otu_subsample")
        self.stat = self.add_tool("meta.otu.otu_taxon_stat")
        self.alpha = self.add_module("meta.alpha_diversity.alpha_group")
        self.beta = self.add_module("meta.beta_diversity.beta_group")
        #self.pan_core = self.add_tool("meta.otu.pan_core_otu")
        self.corr_network_analysis = self.add_module("meta.corr_network_group")
        self.composition = self.add_module("meta.composition.composition_group")
        self.diff = self.add_module("meta.diff.diff_group")
        self.env_analysis = self.add_module("meta.env_group")
        if self.option("group").is_set and self.option("envtable").is_set:
            self.step.add_steps("sample_rename", "framebot","otucluster", "taxassign", "alpha_diversity", "Beta_diversity", "Composition_analysis", "Diff_analysis", "Env_analysis","corr_network_analysis")
        elif self.option("group").is_set:
            self.step.add_steps("sample_rename", "framebot", "otucluster", "taxassign", "alpha_diversity",
                                "Beta_diversity", "Composition_analysis", "Diff_analysis","corr_network_analysis")
        elif self.option("envtable").is_set:
            self.step.add_steps("sample_rename", "framebot", "otucluster", "taxassign", "alpha_diversity",
                                "Beta_diversity", "Composition_analysis", "Env_analysis","corr_network_analysis")
        else:
            self.step.add_steps("sample_rename", "framebot", "otucluster", "taxassign", "alpha_diversity",
                                "Beta_diversity", "Composition_analysis","corr_network_analysis")
        self.spname_spid = dict()
        self.otu_id = None
        self.otu_depth_id = None
        self.env_id = None
        self.env_true_name = {}
        self.level_dict = {'Domain': 1, 'Kingdom': 2, 'Phylum': 3, 'Class': 4, 'Order': 5, 'Family': 6, 'Genus': 7, 'Species': 8, 'otu': 9}
        self.updata_status_api = self.api.meta_update_status
        self.info_path = ""
        self.work_dir_path = ""
        self.function_gene_tool = None
        self.function_gene_path = ""
        self.in_fastq_path = ""
        self.format_taxon = self.add_tool("meta.filecheck.format_taxon")  #guanqingzou 20180801
        self.pan_core_api_status = True
        self.group_dir = self.work_dir + "/all_group"
        self.basic_analysis = []

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'meta')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查参数设置
        """
        # if not self.option("fasta").is_set:
        #     raise OptionError("必须设置输入fasta文件.")
        if self.option("identity") < 0 or self.option("identity") > 1:
            raise OptionError("identity值必须在0-1范围内.", code="12700101")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置序列是否翻转", code="12700102")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件", code="12700103")
        elif self.option('pipeline') == "functional_gene":
            if  self.option('function_gene') == "custom_mode":
                if  not self.option("ref_acid").is_set or not self.option("ref_taxon").is_set:
                    raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件", code="12700103")
                if self.option("seq_style") == "nucl":
                    if not self.option("ref_fasta").is_set:
                        raise OptionError("数据库自定义模式下ASV数据类型为核苷酸必须设置参考核酸fasta序列和参考taxon文件")
            else:
                if self.option("function_gene") not in ['fgr/amoA_archaea_202012','fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012',
                                               'fgr/amoA_comammox_202012', 'fgr/nosZ_202012','fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012',
                                               'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012',
                                               'fgr/pmoA_202012', 'fgr/mmoX_202012']:
                    raise OptionError("数据库%s不被支持", variables=(self.option("function_gene")), code="12700104")
        else:
            if self.option("database") not in ['silva123/16s_bacteria', 'silva123/16s_archaea',
                                               'silva123/16s', 'silva123/18s_eukaryota', 'silva123',
                                               'silva119/16s_bacteria', 'silva119/16s_archaea',
                                               'silva119/16s', 'silva119/18s_eukaryota', 
                                               'unite8.0/its_fungi', 'unite7.2/its_fungi', 'unite7.0/its_fungi',
                                               'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS',
                                               'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA',
                                               'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                                               'maarjam081/AM', 'Human_HOMD', 'Human_HOMD_v15.2', 'Human_HPB', 'Protist_PR2_v4.5',
                                               'silva128/16s_archaea', 'silva128/16s_bacteria',
                                               'silva128/18s_eukaryota', 'silva128/16s',
                                               'silva132/16s_archaea', 'silva132/16s_bacteria',
                                               'silva132/18s_eukaryota', 'silva132/16s',
                                               'silva138/16s_archaea', 'silva138/16s_bacteria',
                                               'silva138/18s_eukaryota', 'silva138/16s',
                                               'greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria',
                                               'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea',
                                               'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s',
                                               'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',
                                               'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria','nt_v20210917/16s',
                                               'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi', "nt_v20210917",
                                               'nt', "nt_v20200604",'fgr/amoA_archaea_202012','fgr/amoA_bacteria_202012','fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012',
                                               'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012',
                                               'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012','fgr/mmoX_202012']:
                    # add by wzy 2016.11.14 silva128,2016.11.23 mrcA,2016.11.28,greengenes135,20170424,amoA(a,b)
                raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="12700104")
        # add by qindanhua 20170112 check if function gene is exist
        if self.option("if_fungene") and self.option("database").split("/")[0] != "fgr":
            raise OptionError("不支持%s功能基因", variables=(self.option("database").split("/")[0]), code="12700105")
        if self.option("group").is_set:
            with open(self.option("group").prop["path"]) as f:
                sample = []
                group1 = []
                for x in f:
                    sample.append(x.strip().split("\t")[0].strip())
                    if x.strip().split("\t")[1].strip():
                        group1.append(x.strip().split("\t")[1].strip())
                if len(sample) != len(group1):
                    self.set_error('第一个分组必须包含所有样本！', code="")
        if self.option("envtable").is_set:
            with open(self.option("envtable").prop["path"]) as f:
                data = f.readlines()
                for i in data[1:]:
                    for x in i.strip().split("\t"):
                        if x == "":
                            self.set_error('环境因子表不能有空值！', code="")
        return True

    # add by qindanhua run function gene tool
    def run_function_gene(self):
        self.step.add_steps("function_gene")
        self.function_gene_tool = self.add_tool("meta.function_gene.function_gene")
        tax_name = self.option("database").split("/")[1]
        if tax_name == 'amoA_archaea':
            tax_name = 'amoA_AOA'
        elif tax_name == 'amoA_bacteria':
            tax_name = 'amoA_AOB'
        self.function_gene_tool.set_options({
            "fasta": self.sample_check.option("otu_fasta"),
            "function_gene": tax_name,
        })
        self.function_gene_tool.on("start", self.set_step, {'start': self.step.function_gene})
        self.function_gene_tool.on("end", self.run_otu)
        self.function_gene_tool.on("end", self.set_step, {'end': self.step.function_gene})
        self.function_gene_tool.run()

    def run_pre_sample_extract(self):
        # if self.option("if_fungene"):
        #     self.in_fastq_path = self.function_gene_tool.output_dir + "/fungene_reads.fastq"
        opts = {
                "in_fastq":  self.option("in_fastq")  # modified by shijin
            }
        self.new_sample_extract.set_options(opts)
        #self.new_sample_extract.on("start", self.set_step, {'start': self.step.sample_rename})
        # self.new_sample_extract.on("end", self.set_step, {'end': self.step.pre_sample_extract})
        self.new_sample_extract.run()

    def set_run(self, opts, module, event, step, start=True):
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)

    def run_sample_rename(self):
        """
        fix by qingchen.zhang @20200922
        改成与qiime2一样的流程，主要是将这个逻辑抽取出来，只做统计
        :return:
        """
        if self.option("info_path") != "" and os.path.exists(self.option("info_path")):  #guanqing.zou 20180814
            info_txt = self.option("info_path")
        else:
            info_txt = self.new_sample_extract.work_dir + "/info.txt"
        opts = {
            # "workdir_sample": self.option("workdir_sample"),  # 从数据库中提取到的样本工作路径，以便对其进行操作
            "info_txt": info_txt,
            # "file_list": self.option("file_list")  # 对样本进行重命名
            "task_id": self._sheet.id,
            "task_type": "meta",
            "query_id": self.option("query_id")
        }
        if self.option("in_fastq").format == "sequence.fastq":
            s3_origin = self.get_origin_name()
            opts['file_origin'] = json.dumps(s3_origin)
        else:
            opts['mapping_file'] = os.path.join(self.work_dir, "remote_input/in_fastq/mapping_file.txt")
        self.sample_rename.set_options(opts)
        self.sample_rename.on("start", self.set_step, {'start': self.step.sample_rename})
        self.sample_rename.run()

    def run_samplecheck(self):
        """
        fix by qingchen.zhang @20200922
        这里只做样本合并和长度分布统计
        :return:
        """
        opts = {
            "info_path": self.sample_rename.output_dir + "/info_path.xls",
            "info_temp": self.sample_rename.output_dir + "/info_temp.xls",
            # "in_fastq": self.split_sample.output_dir,
            "sample_info": self.new_sample_extract.work_dir + "/info.txt"
        }
        if self.option("raw_sequence").is_set:
            opts.update({
                "raw_sequence": self.option("raw_sequence")
            })
        self.sample_check.set_options(opts)
        self.sample_check.on("end", self.set_output, "sample_check")
        self.sample_check.run()

    def run_framebot(self):
        opts = {
            "fasta": self.sample_check.option("otu_fasta"),
            "database": self.option("function_gene"),
            "acid_length": self.option("acid_length"),
            "seq_identity": self.option("seq_identity"),
            "ref_acid": self.option("ref_acid"),
        }
        self.framebot.set_options(opts)
        self.framebot.on("start", self.set_step, {'end': self.step.sample_rename, 'start': self.step.framebot})
        self.framebot.run()


    def run_otu(self):
        if self.option("if_fungene"):
            self.in_fasta_path = self.function_gene_tool.output_dir + "/fungene_reads.fasta"
        elif self.option("pipeline") != "functional_gene":
            self.in_fasta_path = self.sample_check.option("otu_fasta")
        else:
            if self.option("seq_style") == "nucl":
                self.in_fasta_path = self.framebot.output_dir + "/framebot_nucl.fasta"
            else:
                self.in_fasta_path = self.framebot.output_dir + "/framebot_prot.fasta"
        opts = {
            "fasta": self.in_fasta_path,
            # modified by shijin on 20170428
            "identity": self.option("identity")
        }
        self.otu.set_options(opts)
        self.otu.on("end", self.set_output, "otu")
        if self.option("pipeline") != "functional_gene":
            self.otu.on("start", self.set_step, {'end': self.step.sample_rename, 'start': self.step.otucluster})
        else:
            self.otu.on("start", self.set_step, {'start': self.step.otucluster})
        # self.otu.on("end", self.set_step, {'end':self.step.otucluster})
        self.otu.run()

    def run_phylotree(self):
        self.phylo.set_options({
            "fasta_file": self.otu.output_dir + "/otu_reps.fasta"
        })
        # self.phylo.on("start", self.set_step, {'end':self.step.otucluster, 'start':self.step.phylotree})
        self.phylo.on("end", self.set_step, {'end': self.step.otucluster})
        self.phylo.run()

    ###guanqing.zou 20180801
    def run_taxon_format(self):
        if self.option("database") == "custom_mode" or self.option("function_gene") == "custom_mode":
            opts = {"in_taxon_table":self.option('ref_taxon')}
            self.format_taxon.set_options(opts)
            self.format_taxon.on('end',self.run_taxon)
            self.format_taxon.run()
        else:
            self.run_taxon()

    def run_taxon(self):
        opts = {
            "fasta": self.otu.option("otu_rep"),
            "revcomp": self.option("revcomp"),
            "confidence": self.option("confidence"),
            "database": self.option("database")
        }
        if self.option("database") in ["nt", "nt_v20200604","nt_v20210917"]:
            opts.update({
                "query_type": "nucl",
                "blast": "blastn"
            })
        if self.option("pipeline") != "functional_gene":
            if self.option("database") == "custom_mode":
                opts.update({
                    "ref_fasta": self.option("ref_fasta"),
                    #"ref_taxon": self.option("ref_taxon")
                    "ref_taxon": self.format_taxon.option("out_taxon_table")
                })
        else:
            if self.option("function_gene") == "custom_mode":
                if self.option("seq_style") == "nucl":
                    opts.update({
                        "database": "custom_mode",
                        "confidence": self.option("tax_confidence"),
                        "ref_fasta": self.option("ref_fasta"),
                        "ref_taxon": self.format_taxon.option("out_taxon_table")
                    })
                else:
                    opts.update({
                        "database": "customer_mode",
                        "ref_fasta": self.option("ref_acid"),
                        "ref_taxon": self.format_taxon.option("out_taxon_table"),
                        "blast": "blastp",
                        "query_type": "prot",
                        "reference_type": "prot",
                        "evalue": self.option("evalue"),
                        "meta_pipeline": "functional_gene_prot",
                    })
            else:
                opts.update({
                    "confidence": self.option("tax_confidence")
                })
                if self.option("tax_database") in ["nr_v20210917"]:
                    if self.option("seq_style") == "nucl":
                        blast = "blastx"
                        query_type = "nucl"
                    else:
                        blast = "blastp"
                        query_type = "prot"
                    opts.update({
                        "database": self.option("tax_database"),
                        "blast": blast,
                        "query_type": query_type,
                        "evalue": self.option("evalue")
                    })
                elif self.option("seq_style") != "nucl":
                    opts.update({
                        "database": self.option("function_gene"),
                        "blast": "blastp",
                        "query_type": "prot",
                        "evalue": self.option("evalue"),
                        "meta_pipeline": "functional_gene_prot",
                    })

        self.tax.set_options(opts)
        self.tax.on("end", self.set_output, "tax")
        self.tax.on("start", self.set_step, {'start': self.step.taxassign})
        self.tax.on("end", self.set_step, {'end': self.step.taxassign})
        self.tax.run()

    def run_stat(self):
        self.stat.set_options({
            "in_otu_table": self.otu.option("otu_table"),
            "taxon_file": self.tax.option("taxon_file")
        })
        self.stat.on("end", self.set_output, "stat")
        # self.stat.on("end", self.set_step, {'end': self.step.taxassign})
        self.stat.run()

    # 获取抽平后的表名为 OTU_Taxon_Depth
    def run_otu_sub(self):
        self.otu_taxon_dir = self.stat_depth.output_dir
        self.group_all = self.work_dir + "/group_all.txt"
        with open(self.sample_check.output_dir + "/samples_info/samples_info.txt") as r, open(self.group_all,
                                                                                              "w") as t:
            t.write("#sample\tgroup_name\n")
            data = r.readlines()
            for i in data[1:]:
                t.write(i.strip().split("\t")[0] + "\tAll\n")
        with open(self.stat.output_dir + "/tax_summary_a/otu_taxon_otu.full.xls") as f, open(
                self.work_dir + "/otu_taxon_otu.full.xls", "w") as g:
            data1 = f.readlines()
            g.write(data1[0])
            for i in data1[1:]:
                g.write(i.replace(" ", ""))
        self.stat_depth.set_options({
            "in_otu_table": self.work_dir + "/otu_taxon_otu.full.xls",
            "group": self.group_all,
        })
        self.stat_depth.on("end", self.set_output, "stat_depth")
        self.stat_depth.run()

    def run_alpha(self):
        """
        包含alpha分析和 pan_core 分析
        :return:
        """
        opts = {
            'otu_table': self.otu_taxon_dir+"/tax_summary_a/otu_taxon_otu.full.xls",
            "level": self.option('alpha_level'),
            'estimate_indices': self.option('estimate_indices'),
            'rarefy_indices': self.option('rarefy_indices'),
            'rarefy_freq': self.option('rarefy_freq'),
            'grouptable': self.group_dir
        }
        self.alpha = self.add_batch('meta.alpha_diversity.alpha_group', ignore_error=True,
                                   batch_type="module")  # 出错也没关系
        self.set_run(opts, self.alpha, 'alpha_diversity' , self.step.alpha_diversity, False)
        #self.alpha.set_options(opts)
        self.basic_analysis.append(self.alpha)

    def run_beta(self):
        opts = {
            'analysis': 'distance,anosim,' + self.option('beta_analysis'),
            'dis_method': self.option('dis_method'),
            'otutable': self.otu_taxon_dir+"/tax_summary_a/otu_taxon_otu.xls",
            "level": self.option('beta_level'),
            'permutations': self.option('permutations'),
            "grouptable": self.group_dir
        }
        self.logger.info("beta_analysis {}".format('distance,' + self.option('beta_analysis')))
        if self.count_otus:
            opts['phy_newick'] = self.phylo.option('phylo_tre').prop['path']
        self.beta = self.add_batch('meta.beta_diversity.beta_group', ignore_error=True,
                                       batch_type="module")  # 出错也没关系
        self.set_run(opts, self.beta, 'Beta_diversity', self.step.Beta_diversity, False)
        #self.beta.set_options(opts)
        self.basic_analysis.append(self.beta)

    def run_corrnetworkcalc(self):
        """
        单因素相关性网络
        :return:
        """
        opts = {
            'otutable': self.otu_taxon_dir+"/tax_summary_a/otu_taxon_Genus.full.xls",
            'grouptable': self.group_dir,
            }
        self.corr_network_analysis = self.add_batch('meta.corr_network_group', ignore_error=True,
                                          batch_type="module")  # 出错也没关系
        self.set_run(opts, self.corr_network_analysis, 'corr_network_analysis', self.step.corr_network_analysis, False)
        #self.corr_network_analysis.set_options(opts)
        self.basic_analysis.append(self.corr_network_analysis)

    def run_composition(self):
        """
        组成分析 barpie, heatmap, venn, circos
        :return:
        """
        opts = {
            'otutable': self.otu_taxon_dir+"/tax_summary_a/",
            'grouptable': self.group_dir
        }
        self.composition = self.add_batch('meta.composition.composition_group', ignore_error=True, batch_type="module")  # 出错也没关系
        self.set_run(opts, self.composition, 'Composition_analysis', self.step.Composition_analysis, False)
        #self.composition.set_options(opts)
        self.basic_analysis.append(self.composition)

    def run_diff(self):
        """
        物种差异分析 多组, 两组, lefse
        :return:
        """
        opts = {
            'otutable': self.otu_taxon_dir+"/tax_summary_a",
            'grouptable': self.group_dir
        }
        self.diff = self.add_batch('meta.diff.diff_group', ignore_error=True, batch_type="module")  # 出错也没关系
        self.set_run(opts, self.diff, 'Diff_analysis', self.step.Diff_analysis, False)
        #self.diff.set_options(opts)
        self.basic_analysis.append(self.diff)

    def run_env_analysis(self):
        """
        物种差异分析 多组, 两组, lefse
        :return:
        """
        opts = {
            'otutable': self.otu_taxon_dir+"/tax_summary_a/",
            'grouptable': self.group_dir,
            'envtable': self.option("envtable")
        }
        self.env_analysis = self.add_batch('meta.env_group', ignore_error=True, batch_type="module")  # 出错也没关系
        self.set_run(opts, self.env_analysis, 'Env_analysis', self.step.Env_analysis, False)
        #self.env_analysis.set_options(opts)
        self.basic_analysis.append(self.env_analysis)

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def move2outputdir(self, olddir, newname, mode='link'):  # add by shenghe 20160329
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下，如果文件夹名已存在，文件夹会被完整删除。
        """
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="12700101")
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            if os.path.islink(newdir):
                os.remove(newdir)
            else:
                shutil.rmtree(newdir)  # 不可以删除一个链接
        if mode == 'link':
            # os.symlink(os.path.abspath(olddir), newdir)  # 原始路径需要时绝对路径
            shutil.copytree(olddir, newdir, symlinks=True)
        elif mode == 'copy':
            shutil.copytree(olddir, newdir)
        else:
            self.set_error('错误的移动文件方式，必须是\'copy\'或者\'link\'', code="12700102")

    def save_pdf(self):
        all_id = {}
        for i in self.group_id_dict:
            all_id[i] = self.group_id_dict[i]
        self.logger.info("save_pdf {}".format(self.option("save_pdf")))
        if isinstance(self.option("save_pdf"), int):
            self.logger.info("shuizi")
        else:
            self.logger.info("str")
        if self.option("save_pdf"):
            self.logger.info("save_pdf222 {}".format(self.option("save_pdf")))
            #task_info.update_mongo('sg_task', {"task_id": self.task_id}, {"save_pdf": 1})
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": self.sheet.id,
                "project": "meta",
                "interaction": 0,
                "group_id_detail":str(all_id),
                "otu_id":str(self.result_otu_id) if self.option('otu_sub') == "True" else ""
            })
            self.figsave.run()
        else:
            self.end()

    def set_output(self, event):
        # by houshuang 20190925 更新用户选择的数据库至sg_task表 >>>
        data_json = os.path.join(self.work_dir, 'data.json')
        file = open(data_json, 'r')
        js = file.read()
        dic = json.loads(js)
        file.close()
        task_id = str(dic['id'])
        if 'database' in dic['options']:
            database = str(dic['options']['database'])
        else:
            database = str(dic['options']['function_gene'])
        print(database)

        # <<<
        obj = event["bind_object"]
        self.logger.info("obj {}".format(obj.__dict__))
        # 设置QC报告文件
        if event['data'] == "sample_check":
            self.move2outputdir(obj.output_dir, self.output_dir + "/QC_stat")  # 代替cp
            # os.system('cp -r ' + obj.output_dir + ' ' + self.output_dir + "/QC_stat")
        # 设置OTU table文件
        if event['data'] == "otu":
            self.option("otu_table", obj.option("otu_table"))
            self.option("otu_rep", obj.option("otu_rep"))
            self.option("otu_biom", obj.option("otu_biom"))
            self.move2outputdir(obj.output_dir, self.output_dir + "/Otu")  # 代替cp
            # os.system('cp -r ' + obj.output_dir + ' ' + self.output_dir + "/Otu")
            # 设置进化树文件
        if event['data'] == "tax":
            self.option("taxon_file", obj.option("taxon_file"))
            self.move2outputdir(obj.output_dir, self.output_dir + "/Tax_assign")  # 代替cp
            # os.system('cp -r ' + obj.output_dir + ' ' + self.output_dir + "/Tax_assign")
        if event['data'] == "stat":
            self.updata_status_api.add_sample_numbers(task_id, self.count_samples)  # sg_task添加样本数量
            # self.option("otu_taxon_biom", obj.option("otu_taxon_biom"))
            # self.option("otu_taxon_table", obj.option("otu_taxon_table"))
            self.option("otu_taxon_dir", obj.option("otu_taxon_dir"))
            self.move2outputdir(obj.output_dir, self.output_dir + "/OtuTaxon_summary")  # 代替cp
            # os.system('cp -r ' + obj.output_dir + ' ' + self.output_dir + "/OtuTaxon_summary")
        if event['data'] == "stat_depth":
            self.move2outputdir(obj.output_dir, self.output_dir + "/OtuTaxon_summary_depth")  # 代替cp
        if event['data'] == "alpha_diversity":
            self.logger.info(obj.work_dir + "/AlphaGroup/output/Alpha_diversity")
            if os.path.exists(obj.work_dir+"/AlphaGroup/output/Alpha_diversity"):
                self.move2outputdir(obj.work_dir+"/AlphaGroup/output/Alpha_diversity", self.output_dir+"/Alpha_diversity")
            if os.path.exists(obj.work_dir + "/AlphaGroup/output/Pan_Core"):
                self.move2outputdir(obj.work_dir + "/AlphaGroup/output/Pan_Core",self.output_dir + "/Pan_Core")
        if event['data'] == "Beta_diversity":
            self.logger.info(obj.work_dir + "/BetaGroup/output")
            if os.path.exists(obj.work_dir+"/BetaGroup/output/Beta_diversity"):
                self.move2outputdir(obj.work_dir+"/BetaGroup/output/Beta_diversity", self.output_dir+"/Beta_diversity")
        if event['data'] == "corr_network_analysis":
            self.logger.info(obj.work_dir+"/CorrNetworkGroup/output")
            if os.path.exists(obj.work_dir+"/CorrNetworkGroup/output/CorrNetworkSpearmanGenus"):
                self.move2outputdir(obj.work_dir+"/CorrNetworkGroup/output/CorrNetworkSpearmanGenus", self.output_dir+"/CorrNetworkSpearmanGenus")
        if event['data'] == "Composition_analysis":
            self.logger.info(obj.work_dir + "/CompositionGroup/output")
            if os.path.exists(obj.work_dir+"/CompositionGroup/output"):
                self.move2outputdir(obj.work_dir+"/CompositionGroup/output", self.output_dir + "/Composition/")
        if event['data'] == "Diff_analysis":
            self.logger.info(obj.work_dir)
            if os.path.exists(obj.work_dir+"/DiffGroup/output/output/"):
                self.move2outputdir(obj.work_dir+"/DiffGroup/output/output/", self.output_dir + "/DiffGroup")
        if event['data'] == "Env_analysis":
            if os.path.exists(obj.work_dir+"/EnvGroup/output/output/"):
                self.move2outputdir(obj.work_dir+"/EnvGroup/output/output/", self.output_dir + "/Env_analysis")

    def export_base(self):
        data_json = os.path.join(self.work_dir, 'data.json')
        file = open(data_json, 'r')
        js = file.read()
        dic = json.loads(js)
        file.close()
        task_id = str(dic['id'])
        if 'database' in dic['options']:
            database = str(dic['options']['database'])
        else:
            database = str(dic['options']['function_gene'])
        print(database)
        self.updata_status_api.add_database(task_id, database)
        self.updata_status_api.add_save_pdf(task_id, self.option("save_pdf"))
        if self.option("info_path") == "" or not os.path.exists(
                self.option("info_path")):  # add by zhujuan 20180507 工作流中进行检测样品导表
            api_sample = self.api.sample_extract
            data_json = os.path.join(self.work_dir, 'data.json')
            file = open(data_json, 'r')
            js = file.read()
            dic = json.loads(js)
            file.close()
            in_fastq = str(dic['options']['in_fastq']).split("||")
            my_params = dict()
            task_id = self._sheet.id
            my_params["task_id"] = task_id
            my_params["file_info"] = {"file_list": [], "path": in_fastq[-1]}
            my_params["format"] = "sequence.fastq"  # in_fastq[0]
            my_params["query_id"] = str(dic['options']['query_id'])  # add by liulinmeng @20180627
            params = json.dumps(my_params, sort_keys=True, separators=(',', ':'))
            result = api_sample.check_sample(self._sheet.id, self.option("query_id"))
            if result:
                pass
            else:
                se = SE()
                se._config = Config()
                seq_sample_id = se.add_sg_seq_sample(my_params["task_id"], str(my_params["file_info"]), params,
                                                     my_params["query_id"])
                main_id = api_sample.add_sample_check(self._sheet.id, in_fastq[-1], params, self.option("query_id"))
                if self.option("in_fastq").format == "sequence.fastq":
                    s3_origin = self.get_origin_name()
                else:
                    s3_origin = self.get_dir_name()
                api_sample.add_seq_sample_detail(self.new_sample_extract.option("file_sample_list").prop["path"],
                                                 main_id, s3_name=s3_origin)

        api_samples = self.api.sample
        sample_info_path = self.sample_check.output_dir + "/samples_info/samples_info.txt"
        if not os.path.isfile(sample_info_path):
            self.logger.error("找不到报告文件:{}".format(sample_info_path))
            self.set_error("找不到报告文件", code="12700103")
        api_samples.add_samples_info(sample_info_path)
        self.spname_spid = api_samples.get_spname_spid()
        """
        base_info_path = ""

        with open(self.qc.output_dir + "/samples_info/samples_info.txt") as f:
            f.readline()
            for line in f:
                s = line.split('\t')[0]
                base_info_path = self.qc.output_dir + "/base_info/{}.fastq.fastxstat.txt".format(s)
                if not os.path.isfile(base_info_path):
                    raise Exception("找不到报告文件:{}".format(base_info_path))
                api_samples.add_base_info(s, base_info_path)
        """
        for step in (20, 50, 100, 200):
            reads_len_info_path = self.sample_check.output_dir + "/reads_len_info/step_{}.reads_len_info.txt".format(
                str(step))
            if not os.path.isfile(reads_len_info_path):
                self.set_error("找不到报告文件", code="12700103")
            api_samples.add_reads_len_info(step, reads_len_info_path)
        if self.option('envtable').is_set:
            api_env = self.api.env
            self.env_id,self.env_names = api_env.add_env_table(self.option('envtable').prop["path"], self.spname_spid, env_true=self.env_true_name)
        """
        原始序列信息提取新 by sj
        """
        if self.option('raw_sequence').is_set:
            raw_sequence_path = self.sample_check.output_dir + "/raw_sequence.txt"
            if not os.path.exists(raw_sequence_path):
                self.logger.error("找不到报告文件:{}".format(raw_sequence_path))
                self.set_error("找不到报告文件", code="12700103")
            with open(raw_sequence_path, 'r') as m:
                lines = m.readlines()
                if len(lines) > 1:
                    api_samples.add_raw_sequence_info(raw_sequence_path)
        valid_sequence_path = self.sample_check.output_dir + "/valid_sequence.txt"
        if not os.path.exists(valid_sequence_path):
            self.logger.error("找不到报告文件:{}".format(valid_sequence_path))
            self.set_error("找不到报告文件", code="12700103")
        api_samples.add_valid_sequence_info(valid_sequence_path)

        api_otu = self.api.meta
        otu_path = self.output_dir + "/OtuTaxon_summary/otu_taxon.xls"
        rep_path = self.output_dir + "/Otu/otu_reps.fasta"
        otu_path_phy = self.output_dir + "/OtuTaxon_summary/tax_summary_a/otu_taxon_Phylum.full.xls"
        otu_path_gen = self.output_dir + "/OtuTaxon_summary/tax_summary_a/otu_taxon_Genus.full.xls"
        if not os.path.isfile(otu_path):
            self.logger.error("找不到报告文件:{}".format(otu_path))
            self.set_error("找不到报告文件", code="12700103")
        params = {
            "group_id": 'all',
            # "size": 0,
            "size": "",  # modified by hongdongxuan 20170303
            "submit_location": 'otu_statistic',
            "filter_json": "[]",  # add by hongdongxuan 20170303
            "task_type": 'reportTask'
        }
        self.otu_id = api_otu.add_otu_table(otu_path, major=True, rep_path=rep_path, spname_spid=self.spname_spid,
                                            params=params)
        # self.updata_status_api.add_meta_status(table_id=str(self.otu_id), type_name='sg_otu')       # 不更新OTU_Taxon_Origin的信息，by zhengyuan
        api_otu_level = self.api.sub_sample
        api_otu_level.add_sg_otu_detail_level(otu_path, self.otu_id, 9)
        api_otu_level.add_otu_detail(otu_path_phy, self.otu_id, 3)
        api_otu_level.add_otu_detail(otu_path_gen, self.otu_id, 7)
        api_otu_level.add_sg_otu_seq_summary(otu_path, self.otu_id)  # 增加有效序列统计的导表，add by liulinmeng 20180611
        # self.otu_id = str(self.otu_id)
        # self.logger.info('OTU mongo ID:%s' % self.otu_id)
        if self.count_otus:
            api_tree = self.api.newicktree
            tree_path = self.phylo.option('phylo_tre').prop['path']
            if not os.path.isfile(tree_path):
                self.logger.error("找不到报告文件:{}".format(tree_path))
                self.set_error("找不到报告文件", code="12700103")
            if os.path.exists(self.output_dir + '/Otu/otu_phylo.tre'):
                os.remove(self.output_dir + '/Otu/otu_phylo.tre')
            os.link(tree_path, self.output_dir + '/Otu/otu_phylo.tre')
            api_tree.add_tree_file(tree_path, major=True, level=9, table_id=str(self.otu_id), table_type='otu', tree_type='phylo')
        if self.option('group').is_set:
            api_group = self.api.group
            self.group_id_dict,self.group_detail_dict = api_group.add_ini_group_table(self.option('group').prop["path"], self.spname_spid, sort_samples=True)
        else:
            self.group_id_dict={"All":"all"}
            all_id =[]
            for x in self.spname_spid.values():
                all_id.append(str(x))
            self.group_detail_dict = {"All":{"All":all_id}}
        self.logger.info('group_id_dict:%s' % self.group_id_dict)
        self.logger.info('group_detail_dict:%s' % self.group_detail_dict)
        if self.option('otu_sub') == "True":
            api_otu = self.api.meta
            otu_path = self.output_dir + "/OtuTaxon_summary_depth/otu_taxon.xls"
            rep_path = self.output_dir + "/Otu/otu_reps.fasta"
            if not os.path.isfile(otu_path):
                self.logger.error("找不到报告文件:{}".format(otu_path))
                self.set_error("找不到报告文件", code="12700103")
            params = {
                "group_id": 'all',
                "size": "min",
                "submit_location": 'otu_statistic',
                "filter_json": "[]",
                "task_type": ''
            }
            self.otu_depth_id = api_otu.add_otu_table(otu_path, major=True, rep_path=rep_path, from_out_table=self.otu_id, spname_spid=self.spname_spid,
                                                name="OTU_Taxon_Depth", params=params,level_id=[9])
            # self.updata_status_api.add_meta_status(table_id=str(self.otu_id), type_name='sg_otu')       # 不更新OTU_Taxon_Origin的信息，by zhengyuan
            api_otu_level = self.api.sub_sample
            api_otu_level.add_sg_otu_detail_level(otu_path, self.otu_depth_id, 9)
            api_otu_level.add_sg_otu_seq_summary(otu_path, self.otu_depth_id)
            self.result_otu_id = self.otu_depth_id
        else:
            self.result_otu_id = self.otu_id
        self.logger.info("self.result_otu_id:::{}".format(self.result_otu_id))

    def run_basic_analysis(self):
        """
        运行基础分析
        :return:
        """
        self.run_alpha()
        if self.option('group').is_set and not self.option('envtable').is_set:
            self.run_beta()
            self.run_corrnetworkcalc()
            self.run_composition()
            self.run_diff()
            #self.on_rely([self.alpha, self.beta, self.composition, self.corr_network_analysis, self.diff], self.set_output,"basic_analysis")
        elif self.option('envtable').is_set and not self.option('group').is_set:
            self.run_beta()
            self.run_corrnetworkcalc()
            self.run_composition()
            if self.count_samples >= 3:
                self.run_env_analysis()
                #self.on_rely([self.alpha, self.beta, self.composition, self.corr_network_analysis, self.env_analysis], self.set_output,"basic_analysis")
            #else:
                #self.on_rely([self.alpha, self.beta, self.composition, self.corr_network_analysis,], self.set_output, "basic_analysis")
        elif self.option('envtable').is_set and self.option('group').is_set:
            self.run_beta()
            self.run_corrnetworkcalc()
            self.run_composition()
            self.run_diff()
            if self.count_samples >= 3:
                self.run_env_analysis()
                #self.on_rely([self.alpha, self.beta, self.composition, self.corr_network_analysis, self.diff, self.env_analysis], self.set_output,"basic_analysis")
            #else:
                #self.on_rely([self.alpha, self.beta, self.composition, self.corr_network_analysis, self.diff,], self.set_output, "basic_analysis")
        else:
            self.run_beta()
            self.run_corrnetworkcalc()
            self.run_composition()
        self.on_rely(self.basic_analysis,self.run_api)
        for module in self.basic_analysis:
            module.run()
            #gevent.sleep(0)

    def check_otu_run(self):
        """
        在不同的OTU数目以及不同的样本数量下，有些分析会被跳过不做
        """
        self.count_samples = 0
        if self.option("group").is_set:
            group_file_spilt(self.option("group").path,self.group_dir)
            with open(self.sample_check.output_dir + "/samples_info/samples_info.txt") as r:
                self.count_samples = len(r.readlines()) - 1
        else:
            if not os.path.exists(self.group_dir):
                os.mkdir(self.group_dir)
            with open(self.sample_check.output_dir + "/samples_info/samples_info.txt") as r, open(self.group_dir+"/All.group.txt", "w") as t:
                t.write("#sample\tgroup_name\n")
                data = r.readlines()
                self.count_samples = len(data) - 1
                for i in data[1:]:
                    t.write(i.strip().split("\t")[0] + "\tAll\n")
        self.update_info = ""
        self.count_otus = True  # otu/代表序列数量大于3 ,change by wzy from 2 to 3
          # 样本数量是否大于等于2
        counts = 0
        for i in open(self.otu.output_dir + '/otu_reps.fasta'):
            if i[0] == '>':
                counts += 1
                if counts > 3:    # change from 2 to 3 by wzy
                    break
        else:
            self.count_otus = False
        self.logger.info(self.otu.output_dir + '/otu_reps.fasta')
        self.logger.info(counts)
        self.logger.info(self.count_otus)

        if self.count_samples < 3:  # 少于3个样本
            #self.update_info += "样本少于三个，不进行beta多样性相关分析；"
            self.update_info += "Sample size too small: SKIP Beta diversity analysis!"
            self.option('beta_analysis', '')
        if not self.count_otus:
            if 'pca' in self.option('beta_analysis'):
                self.option('beta_analysis', self.option('beta_analysis').replace('pca', '').strip(',').replace(',,', ','))
                #self.update_info += "OTU数量过少，不做PCA分析；
                self.update_info += "OTU size too small: SKIP PCA analysis!"
            if "unifrac" in self.option('dis_method'):
                self.option('dis_method', "bray_curtis")
                #self.update_info += "OTU数量过少，不使用unifrac类型距离算法，采用默认bray_curtis算法；"
                self.update_info += "OTU size too small: UniFrac distance method replaced by bray_curtis!"
            indices = self.option("estimate_indices").split(',')
            if 'pd' in indices:
                indices.remove("pd")
                indices = ','.join(indices)
                self.option('estimate_indices', indices)
                #self.update_info += "OTU数量过少没有进行进化树分析，所以移除了依赖进化树分析多样性指数 PD；"
                self.update_info += "OTU size too small: SKIP PD index analysis!"
        if self.count_otus and self.count_samples > 1:
            self.on_rely([self.tax, self.phylo], self.run_stat)
            if self.option("otu_sub") == "True":
                self.stat.on('end', self.run_otu_sub)
                self.stat_depth.on('end', self.run_basic_analysis)
            else:
                self.otu_taxon_dir = self.stat.output_dir
                self.stat.on('end', self.run_basic_analysis)
            #self.run_taxon()
            self.run_taxon_format()  #zouguanqing
            self.run_phylotree()
        elif self.count_otus and self.count_samples == 1:
            self.pan_core_api_status = False
            #self.update_info += "样本数量过少，不进行Pan Core分析；"
            self.update_info += "Sample size too small: SKIP Pan_Core analysis!"
            self.on_rely([self.tax, self.phylo], self.run_stat)
            if self.option("otu_sub") == "True":
                self.stat.on('end', self.run_otu_sub)
                self.stat_depth.on('end', self.run_basic_analysis)
            else:
                self.otu_taxon_dir = self.stat.output_dir
                self.stat.on('end', self.run_basic_analysis)
            #self.run_taxon()
            self.run_taxon_format()  #zouguanqing
            self.run_phylotree()
        else:
            #self.update_info += "OTU数量过少，不进行物种进化树分析，将不会生成物种进化树；"
            self.pan_core_api_status = False
            self.update_info += "OTU size too small: SKIP Phylogenetic analysis!"
            self.tax.on('end', self.run_stat)
            if self.option("otu_sub") == "True":
                self.stat.on('end', self.run_otu_sub)
                self.stat_depth.on('end', self.run_basic_analysis)
            else:
                self.otu_taxon_dir = self.stat.output_dir
                self.stat.on('end', self.run_basic_analysis)
            #self.run_taxon()
            self.run_taxon_format()  #zouguanqing
        if self.update_info:
            self.logger.info("分析结果异常处理：{}".format(self.update_info))
            self.step._info += self.update_info
            self.step._has_state_change = True

    def change_env_name(self):
        '''
        环境因子改名
        '''
        self.logger.info("change_env_name!!!！")
        tb = pd.read_csv(self.option('envtable').path, sep='\t')
        not_word = re.compile(r'\W')
        #if not filter(lambda x: not_word.search(x), tb.columns[1:]):
        #    return
        new_name = map(lambda x: 'env_' + str(x), range(1, len(tb.columns)))
        self.env_true_name = {new_name[i]: tb.columns[i+1] for i in range(len(new_name))}
        ret = list(tb.columns)
        ret[1:] = new_name
        tb.columns = ret
        new_table = os.path.join(self.work_dir, 'new_envtable.xls')
        self.logger.info("new_table{}!!!！".format(new_table))
        tb.to_csv(new_table, sep='\t', index=False)
        self.option('envtable').set_path(new_table)

    def run_api(self):
        """
        导表
        :return:
        """
        self.logger.info("开始运行导表！")
        if self.option('otu_sub') == "True":
            self.table_name = "Depth_"
        else:
            self.table_name = "Origin_"
        self.export_base()
        if os.path.exists(os.path.join(self.output_dir, "Alpha_diversity")):
            self.export_alpha_analysis()
        if os.path.exists(os.path.join(self.output_dir, "Pan_Core")):
            self.export_pan_core()
        if os.path.exists(os.path.join(self.output_dir, "Beta_diversity")):
            self.export_beta_analysis()
        if os.path.exists(os.path.join(self.output_dir, "CorrNetworkSpearmanGenus")):
            self.export_corr_network()
        if os.path.exists(os.path.join(self.output_dir, "Composition")):
            self.export_composition()
        if os.path.exists(os.path.join(self.output_dir, "DiffGroup")):
            self.export_diff_analysis()
        if os.path.exists(os.path.join(self.output_dir, "Env_analysis")):
            self.export_env_analysis()
        with open(self.output_dir+"/运行参数.txt","w") as t:
            if self.option("pipeline") != "functional_gene":
                t.write("OTU序列相似度: {}\n".format(self.option('identity')))
                t.write("物种分类数据库: {}\n".format(self.option('database')))
                if self.option('database') not in ['nt', "nt_v20200604"]:
                    t.write("分类置信度: {}\n".format(self.option('confidence')))
        self.save_pdf()
        #self.end()

    def export_alpha_analysis(self):
        for group_dir in os.listdir(self.output_dir + "/Alpha_diversity/"):
            api_est = self.api.estimator
            est_path = self.output_dir + "/Alpha_diversity/" + group_dir +"/Estimators/estimators.xls"
            if not os.path.isfile(est_path):
                self.logger.error("找不到报告文件:{}".format(est_path))
                self.set_error("找不到报告文件", code="12700103")
            indice = sorted(self.option("estimate_indices").split(','))
            level_id = self.level_dict[self.option('alpha_level')]
            params = {
                "level_id": level_id,
                "index_type": ','.join(indice),
                'submit_location': 'alpha_diversity_index',
                'task_type': 'reportTask',
                'group_id': str(self.group_id_dict[group_dir]),
                'group_detail': self.group_detail_dict[group_dir]
            }
            self.logger.info("params:::{}".format(params))
            est_id = api_est.add_est_table(est_path, major=True, level=level_id, otu_id=str(self.result_otu_id),name="Estimators_"+self.table_name+group_dir,
                                       params=params,index_type=','.join(indice))
            self.updata_status_api.add_meta_status(table_id=str(est_id), type_name='sg_alpha_diversity')
            # 主表写入没有加name，所以此处table_name固定


            api_rare = self.api.rarefaction
            rare_path = self.output_dir + "/Alpha_diversity/" + group_dir + "/Rarefaction/"
            # indice = sorted(self.option("rarefy_indices").split(','))
            indice = sorted(self.option("rarefy_indices").split(','))
            name = "Rarefaction_"+self.table_name+group_dir
            if self.option("rarefy_indices") != "":
                params = {
                    "freq": 100,
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': str(self.group_id_dict[group_dir]),
                    "index_type": ','.join(indice),
                    'level_id': 9,
                    'otu_id': str(self.result_otu_id),
                    'submit_location': 'alpha_rarefaction_curve',
                    "task_type": "reportTask"
                }
                rare_id = api_rare.add_rare_table(rare_path, level=9, otu_id=str(self.result_otu_id),
                                                  params=params,name=name)
                self.updata_status_api.add_meta_status(table_id=str(rare_id), type_name='sg_alpha_rarefaction_curve')
                # 主表写入没有加name，所以此处table_name固定

            if os.path.exists(self.output_dir + "/Alpha_diversity/" + group_dir + "/EstTTest/"):
                api_est_t_test = self.api.est_t_test
                name = "EstTTest_" + self.table_name+group_dir
                est_path = self.output_dir + "/Alpha_diversity/" + group_dir + "/EstTTest/"
                params = {
                    "alpha_diversity_id": str(est_id),
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': str(self.group_id_dict[group_dir]),
                    'otu_id': str(self.result_otu_id),
                    'submit_location': 'alpha_ttest',
                    'task_type': 'reportTask',
                    'test_method': 'mann'
                }
                all_key = self.group_detail_dict[group_dir].keys()
                all_compare = []
                for i in all_key:
                    for y in all_key:
                        if i!=y:
                            if i + "|" + y in all_compare or y + "|" + i in all_compare:
                                pass
                            else:
                                all_compare.append(i + "|" + y)
                compare_column= ",".join(all_compare)
                self.logger.info("compare_column:::{}".format(compare_column))
                est_t_test_id = api_est_t_test.add_est_table(est_path, level=level_id,est_id=est_id,group_name=group_dir,otu_id=str(self.result_otu_id),params=params,compare_column=compare_column,name=name)
                self.updata_status_api.add_meta_status(table_id=str(est_t_test_id), type_name='sg_alpha_ttest')

    def export_pan_core(self):
        for group_dir in os.listdir(self.output_dir + "/Pan_Core/"):
            api_pan_core = self.api.pan_core
            name = "Pan_"+self.table_name+group_dir
            params = {
                "level_id": 9,
                'group_id': self.group_id_dict[group_dir],
                'group_detail': self.group_detail_dict[group_dir],
                "submit_location": "otu_pan_core",
                "otu_id": str(self.result_otu_id),
            }
            pan_id = api_pan_core.create_pan_core_table(1, json.dumps(params), self.group_id_dict[group_dir], 9, self.result_otu_id, name,"end")
            name = "Core_"+self.table_name+group_dir
            core_id = api_pan_core.create_pan_core_table(2, json.dumps(params), self.group_id_dict[group_dir], 9, self.result_otu_id, name,"end")
            pan_path = self.output_dir + "/Pan_Core/"+group_dir+"/PanCoreOTU/pan.richness.xls"
            core_path = self.output_dir + "/Pan_Core/"+group_dir+"/PanCoreOTU/core.richness.xls"
            api_pan_core.add_pan_core_detail(pan_path, pan_id)
            api_pan_core.add_pan_core_detail(core_path, core_id)
            self.updata_status_api.add_meta_status(table_id=pan_id, type_name='sg_otu_pan_core')
            self.updata_status_api.add_meta_status(table_id=core_id, type_name='sg_otu_pan_core')

    def export_beta_analysis(self):
        for group_dir in os.listdir(self.output_dir + "/Beta_diversity/"):
            api_dist = self.api.distance
            dist_path = self.output_dir + "/Beta_diversity/" + group_dir + "/Distance/"+ os.listdir(self.output_dir + "/Beta_diversity/" + group_dir + "/Distance/")[0]
            if not os.path.isfile(dist_path):
                self.logger.error("找不到报告文件:{}".format(dist_path))
                self.set_error("找不到报告文件", code="12700103")
            level_id = self.level_dict[self.option('beta_level')]
            params = {
                # 'otu_id': str(self.otu_id),  # 在metabase中不能执行，生成self.otu_id的api可能会被截取
                'level_id': level_id,
                'distance_algorithm': self.option('dis_method'),
                'submit_location': 'beta_sample_distance_hcluster_tree',  # 为前端分析类型标识
                'task_type': 'reportTask',
                'hucluster_method': self.option('linkage'),
                'group_id': self.group_id_dict[group_dir],
                'group_detail': self.group_detail_dict[group_dir],
            }
            dist_id = api_dist.add_dist_table(dist_path, level=level_id, otu_id=self.result_otu_id, major=True, params=params)
            # self.updata_status_api.add_meta_status(table_id=str(dist_id), type_name='sg_beta_specimen_distance')  # 主表写入没有加name，所以此处table_name固定
            if 'hcluster' in self.option('beta_analysis').split(','):
                # 设置hcluster树文件
                api_hcluster = self.api.newicktree
                hcluster_path = self.output_dir + "/Beta_diversity/" + group_dir + "/Hcluster/hcluster.tre"
                if not os.path.isfile(hcluster_path):
                    self.logger.error("找不到报告文件:{}".format(hcluster_path))
                    self.set_error("找不到报告文件", code="12700103")
                name = "Tree_" +self.table_name+group_dir
                tree_id = api_hcluster.add_tree_file(hcluster_path, major=True, table_id=str(self.result_otu_id),
                                                     level=level_id,
                                                     table_type='otu', tree_type='cluster', params=params,
                                                     update_dist_id=dist_id,name=name)
                final_rank_others_file=self.output_dir + "/Beta_diversity/" + group_dir + "/Hcluster/barplot_table.xls"
                api_hcluster.update_newick(final_rank_others_file, tree_id)
                api_hcluster.add_newick_detail(final_rank_others_file, tree_id)
                self.updata_status_api.add_meta_status(table_id=str(tree_id),
                                                       type_name='sg_newick_tree')  # 主表写入没有加name，所以此处table_name固定


            beta_multi_analysis_dict = {'pca': 'beta_multi_analysis_pca', 'pcoa': 'beta_multi_analysis_pcoa',
                                        'nmds': 'beta_multi_analysis_nmds', 'dbrda': 'beta_multi_analysis_dbrda',
                                        'rda_cca': 'beta_multi_analysis_rda_cca'}  # 为前端分析类型标识
            for ana in self.option('beta_analysis').split(','):
                if ana in ['pca', 'pcoa', 'nmds', 'dbrda', 'rda_cca']:
                    api_betam = self.api.beta_multi_analysis
                    params = {
                        # 'otu_id': str(self.otu_id),  # 在metabase中不能执行，生成self.otu_id的api可能会被截取
                        'level_id': level_id,
                        'analysis_type': ana,
                        'submit_location': beta_multi_analysis_dict[ana],
                        'task_type': 'reportTask',
                        'group_id': self.group_id_dict[group_dir],
                        'group_detail': self.group_detail_dict[group_dir],
                        'diff_test_method':'anosim',
                        'change_times':'999'
                    }
                    #if self.option('envtable').is_set:
                    #    # params['env_id'] = str(self.env_id)  # 在metabase中不能执行，生成self.env_id的api可能会被截取
                    #    params['env_labs'] = ','.join(self.option('envtable').prop['group_scheme'])
                    if ana in ['pcoa', 'nmds', 'dbrda']:
                        params['distance_algorithm'] = self.option('dis_method')
                    if ana in ['pca']:
                        params['scale'] = "T"
                    main_id = api_betam.add_beta_multi_analysis_result(dir_path=self.output_dir + "/Beta_diversity/" + group_dir, analysis=ana,
                                                                       main=True, otu_id=self.result_otu_id, params=params,group_name=self.table_name+group_dir)
                    self.updata_status_api.add_meta_status(table_id=main_id,
                                                           type_name='sg_beta_multi_analysis')  # 主表写入没有加name，所以此处table_name固定
                    self.logger.info('set output beta %s over.' % ana)
            ## anosim
            if os.path.exists(self.output_dir + "/Beta_diversity/" + group_dir + "/Anosim/"):
                api_anosim = self.api.anosim
                name = "Anosim_Adonis_"+self.table_name+group_dir
                params = {
                    'distance_algorithm': "bray_curtis",
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    'level_id': 9,
                    'otu_id': str(self.result_otu_id),
                    'permutations': 999,
                    'submit_location': "beta_multi_analysis_anosim",
                    'task_type': "reportTask"
                }
                anosim_id=api_anosim.add_beta_anosim_main(name,self.result_otu_id,json.dumps(params))
                api_anosim.add_beta_anosim_result(self.output_dir + "/Beta_diversity/" + group_dir, main=False, main_id=anosim_id)


    def export_corr_network(self):
        for group_dir in os.listdir(self.output_dir + "/CorrNetworkSpearmanGenus/"):
            if os.path.isdir(self.output_dir + "/CorrNetworkSpearmanGenus/"+group_dir):
                name = "CorrNetworkSpearmanGenus_"+self.table_name + group_dir
                params = {
                    "abundance": 50,
                    "coefficient": 0.5,
                    "color_level": 3,
                    'group_id': self.group_id_dict[group_dir],
                    'group_detail': self.group_detail_dict[group_dir],
                    "lable": 0.03,
                    "level_id": 7,
                    "otu_id": str(self.result_otu_id),
                    "ratio_method": "spearman",
                    "significance": 0.05,
                    "submit_location": "corr_network_analyse",
                    "task_type": "reportTask",
                }
                api_corrnetwork = self.api.corr_network
                corr_id = api_corrnetwork.create_corrnetwork(json.dumps(params), self.group_id_dict[group_dir],
                                                             str(self.result_otu_id), name=name, level_id=7)
                node_links_path = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_by_cut.txt'
                node_abundance_path = self.output_dir + "/CorrNetworkSpearmanGenus/species_abundance.txt"
                network_clustering_path = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_clustering.txt'
                network_degree_path = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_node_degree.txt'
                network_centrality_path = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_centrality.txt'
                network_attributes_path = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_attributes.txt'
                network_degree_distribution = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/corr_network_calc/corr_network_degree_distribution.txt'
                profile1 = self.output_dir + "/CorrNetworkSpearmanGenus/newtable.txt",
                attributes_file = self.output_dir + "/CorrNetworkSpearmanGenus/" + "/corr_network_calc/corr_network_attributes.txt"
                self.logger.info('corr_id corr_id %s over.' % corr_id)
                api_corrnetwork.add_network_links_table(file_path=node_links_path, node_id_file=network_degree_path,
                                                        table_id=corr_id)
                api_corrnetwork.add_network_abundance_table(file_path=node_abundance_path,
                                                            table_id=corr_id)
                api_corrnetwork.add_network_cluster_degree(file1_path=network_degree_path,
                                                           file2_path=network_clustering_path,
                                                           table_id=corr_id)
                api_corrnetwork.add_network_centrality(file_path=network_centrality_path,
                                                       table_id=corr_id)
                api_corrnetwork.add_network_attributes(file_path=network_attributes_path,
                                                       table_id=corr_id)
                api_corrnetwork.add_network_degree_distribution(file_path=network_degree_distribution,
                                                                table_id=corr_id)
                corr_file = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/heatmap/corr.xls'
                p_file = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/heatmap/pvalue.xls'
                tree_file = self.output_dir + "/CorrNetworkSpearmanGenus/" + group_dir + '/heatmap/corr.cluster_tree.xls'
                if not os.path.exists(tree_file):
                    tree_file = None
                api_corrnetwork.add_heatmap_corr_detail(corr_id, corr_file, p_file,
                                                        tree_file=tree_file)
                self.updata_status_api.add_meta_status(table_id=corr_id, type_name='sg_corr_network')


    def export_composition(self):
        for group_dir in os.listdir(self.output_dir + "/Composition/"):
            ## barpie
            api_cluster_analysis = self.api.cluster_analysis
            name = "CommunityBarPieGenus_"+self.table_name + group_dir
            params = {
                "combine_value": 0.01,
                'group_id': self.group_id_dict[group_dir],
                'group_detail': self.group_detail_dict[group_dir],
                "group_method": "",
                "level_id": 7,
                "otu_id": str(self.result_otu_id),
                "submit_location": "otu_group_analyse",
                "task_type": "reportTask",
            }
            bar_id_genus = api_cluster_analysis.add_sg_otu(json.dumps(params),self.result_otu_id,name=name)
            api_cluster_analysis.add_sg_otu_detail(self.output_dir + "/Composition/"+group_dir+"/BarPie_Genus/taxa.percents.table.xls", bar_id_genus,self.result_otu_id)
            self.updata_status_api.add_meta_status(table_id=bar_id_genus,type_name='sg_otu')

            api_cluster_analysis = self.api.cluster_analysis
            name = "CommunityBarPiePhylum_"+self.table_name + group_dir
            params = {
                "combine_value": 0.01,
                'group_id': self.group_id_dict[group_dir],
                'group_detail': self.group_detail_dict[group_dir],
                "group_method": "",
                "level_id": 3,
                "otu_id": str(self.result_otu_id),
                "submit_location": "otu_group_analyse",
                "task_type": "reportTask",
            }
            bar_id_phy = api_cluster_analysis.add_sg_otu(json.dumps(params), self.result_otu_id, name=name)
            api_cluster_analysis.add_sg_otu_detail(self.output_dir + "/Composition/"+group_dir+"/BarPie_Phylum/taxa.percents.table.xls",
                                                   bar_id_phy, self.result_otu_id)
            self.updata_status_api.add_meta_status(table_id=bar_id_phy,type_name='sg_otu')

            ## Circos
            api_composition = self.api.composition
            name = "Circos_"+self.table_name + group_dir
            params = {
                "combine_value": 0.01,
                "graphic_type": "circos",
                'group_detail': self.group_detail_dict[group_dir],
                'group_id': self.group_id_dict[group_dir],
                "group_method": "",
                "average": "",
                "level_id": 3,
                "otu_id": str(self.result_otu_id),
                "submit_location": "circos",
                "task_type": "reportTask",
            }
            composition_id = api_composition.add_sg_composition(json.dumps(params, sort_keys=True, separators=(',', ':')), self.result_otu_id, name=name)
            api_composition.add_sg_otu_detail(
                self.output_dir + "/Composition/" + group_dir + "/Circos/taxa.percents.table.xls",
                composition_id, self.result_otu_id,"circos")
            self.updata_status_api.add_meta_status(table_id=composition_id, type_name='sg_composition',submit_location="composition")


            ## Heatmap
            api_heatmap = self.api.hierarchical_clustering_heatmap
            table_name = "CommunityHeatmap_"+self.table_name + group_dir
            params = {
                "add_Algorithm": "",
                'group_detail': self.group_detail_dict[group_dir],
                'group_id': self.group_id_dict[group_dir],
                "level_color": 3,
                "level_id": 7,
                "method": "average",
                "otu_id": str(self.result_otu_id),
                "sample_method": "average",
                "species_number": "50",
                "submit_location": "hc_heatmap",
                "task_type": "reportTask",
            }
            self.color_dict={}
            species_list=[]
            sample_list=[]
            with open(self.output_dir + "/Composition/" + group_dir + "/Heatmap/sort_sample.xls", 'r') as r, open(
                    self.work_dir + "/tmp_heatmap.xls", 'w') as w:
                for n in r:
                    n = n.strip("\n")
                    brr = n.strip().split(";")
                    if brr[0].startswith("OTU"):
                        w.write(brr[0] + "\n")
                    else:
                        w.write(brr[-1] + "\n")
                        name = re.split("\t", brr[-1])
                        self.color_dict[name[0]] = brr[2].strip()
                new_otu_file_path = self.work_dir + "/tmp_heatmap.xls"
            sample_tree=""
            species_tree=""
            species_list=""
            sample_list=""
            if os.path.exists(self.output_dir + "/Composition/" + group_dir + "/Heatmap/species_hcluster.tre"):
                species_tree_path = self.output_dir + "/Composition/" + group_dir + "/Heatmap/species_hcluster.tre"
                if os.path.exists(species_tree_path):
                    with open(species_tree_path, "r") as f:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [i[1] for i in raw_samp]
            if os.path.exists(self.output_dir + "/Composition/" + group_dir + "/Heatmap/sample_hcluster.tre"):
                sample_tree_path = self.output_dir + "/Composition/" + group_dir + "/Heatmap/sample_hcluster.tre"
                if os.path.exists(sample_tree_path):
                    with open(sample_tree_path, "r") as f:
                        sample_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', sample_tree)
                        sample_list = [i[1] for i in raw_samp]
            composition_id = api_heatmap.add_sg_hc_heatmap(json.dumps(params, sort_keys=True, separators=(',', ':')), self.result_otu_id, name=table_name)
            api_heatmap.add_sg_hc_heatmap_detail(new_otu_file_path, self.color_dict, composition_id,
                                             self.result_otu_id,
                                             sample_tree=sample_tree, sample_list=sample_list,
                                             species_tree=species_tree, species_list=species_list,
                                             otu_relative=self.output_dir + "/Composition/" + group_dir + "/Heatmap/heatmap.taxa.relative.xls")
            self.updata_status_api.add_meta_status(table_id=composition_id, type_name='sg_hc_heatmap')

            ## Venn
            if os.path.exists(self.output_dir + "/Composition/" + group_dir + "/Venn"):
                api_venn = self.api.venn
                name = "Venn_"+self.table_name + group_dir
                params = {
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    "level_id": 9,
                    "otu_id": str(self.result_otu_id),
                    "submit_location": "otu_venn",
                    "task_type": "reportTask",
                }
                venn_path = self.output_dir + "/Composition/" + group_dir + "/Venn/venn_table.xls"
                venn_graph_path = self.output_dir + "/Composition/" + group_dir + "/Venn/venn_graph.xls"
                asv_table_path=self.output_dir + "/Composition/" + group_dir + "/Venn/sort_samples.xls"
                venn_id = api_venn.create_venn_table(json.dumps(params), self.group_id_dict[group_dir], 9, self.result_otu_id,
                                                     name=name)
                api_venn.add_venn_detail(venn_path, venn_id, self.result_otu_id, 9,group_path=self.group_dir+"/"+group_dir+".group.txt")
                api_venn.add_venn_graph(venn_graph_path, venn_id)
                api_venn.add_venn_pie(venn_path, venn_id, asv_table_path, self.group_dir+"/"+group_dir+".group.txt")
                self.updata_status_api.add_meta_status(table_id=venn_id, type_name='sg_otu_venn')


    def export_diff_analysis(self):
        for group_dir in os.listdir(self.output_dir + "/DiffGroup/"):
            for analysis_dir in os.listdir(self.output_dir + "/DiffGroup/"+group_dir):
                ## Diff multi
                if "DiffStatMultiple" in analysis_dir:
                    if "Genus" in analysis_dir:
                        level_id = 7
                    else:
                        level_id = 3
                    name = analysis_dir + "_"+self.table_name+group_dir
                    api_multiple = self.api.stat_test
                    params = {
                        "correction": "fdr",
                        "coverage": 0.95,
                        'group_detail': self.group_detail_dict[group_dir],
                        'group_id': self.group_id_dict[group_dir],
                        "level_id": level_id,
                        "methor": "tukeykramer",
                        "otu_id": str(self.result_otu_id),
                        "submit_location": "species_difference_multiple",
                        "task_type": "reportTask",
                        "test": "kru_H",
                    }
                    stat_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + "/kru_H_result.xls"
                    boxfile_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + "/kru_H_boxfile.xls"
                    bar_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/kru_H_plot_group_bar.xls'
                    category_name = ",".join(self.group_detail_dict[group_dir].keys())
                    cifiles = []
                    for r, d, f in os.walk(self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir):
                        for i in f:
                            if "tukeykramer" in i:
                                ci_path = r + '/' + i
                                if not os.path.isfile(ci_path):
                                    self.logger.error("找不到报告文件:{}".format(ci_path))
                                    self.set_error("找不到报告文件", code="12702701")
                                cifiles.append(ci_path)

                    multi_id = api_multiple.creat_multi_table(self.result_otu_id,self.group_id_dict[group_dir],json.dumps(params),level_id,
                                                              category_name=category_name,name=name)
                    api_multiple.add_species_difference_check_detail(statfile=stat_path, cifiles=cifiles,
                                                                     table_id=multi_id,
                                                                     level=level_id, check_type='multiple',
                                                                     params=json.dumps(params),
                                                                     category_name=category_name,
                                                                     group_id=self.group_id_dict[group_dir],
                                                                     from_otu_table=str(self.result_otu_id), major=False,
                                                                     posthoc="tukeykramer")
                    api_multiple.add_species_difference_check_boxplot(boxfile_path, multi_id)
                    api_multiple.add_species_difference_check_barplot(bar_path, multi_id)
                    self.updata_status_api.add_meta_status(table_id=multi_id, type_name='sg_species_difference_check')

                ## two
                if "DiffStatTwoGrou" in analysis_dir:
                    api_two_group = self.api.stat_test
                    name = analysis_dir + "_"+self.table_name + group_dir
                    if "Genus" in analysis_dir:
                        level_id = 7
                    else:
                        level_id = 3
                    params = {
                        "correction": "fdr",
                        "coverage": 0.95,
                        'group_detail': self.group_detail_dict[group_dir],
                        'group_id': self.group_id_dict[group_dir],
                        "level_id": level_id,
                        "methor": "bootstrap",
                        "otu_id": str(self.result_otu_id),
                        "submit_location": "species_difference_two_group",
                        "task_type": "reportTask",
                        "test": "mann",
                        "type": "two.side",
                    }
                    stat_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/mann_result.xls'
                    boxfile_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/mann_boxfile.xls'
                    ci_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/mann_CI.xls'
                    bar_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/mann_plot_group_bar.xls'
                    if not os.path.isfile(stat_path):
                        self.logger.error("找不到报告文件:{}".format(stat_path))
                        self.set_error("找不到报告文件", code="12704101")
                    if not os.path.isfile(boxfile_path):
                        self.logger.error("找不到报告文件:{}".format(boxfile_path))
                        self.set_error("找不到报告文件", code="12704101")
                    if not os.path.isfile(ci_path):
                        self.logger.error("找不到报告文件:{}".format(ci_path))
                        self.set_error("找不到报告文件", code="12704101")
                    category_name = ",".join(self.group_detail_dict[group_dir].keys())
                    two_id = api_two_group.creat_two_table(self.result_otu_id,self.group_id_dict[group_dir],json.dumps(params),level_id,
                                                              category_name=category_name,name=name)
                    api_two_group.add_species_difference_check_detail(statfile=stat_path, cifiles=[ci_path],
                                                                      table_id=two_id,
                                                                      level=level_id,
                                                                      check_type='two_group',
                                                                      params=json.dumps(params),
                                                                      category_name=category_name,
                                                                      group_id=self.group_id_dict[group_dir],
                                                                      from_otu_table=str(self.result_otu_id),
                                                                      posthoc=None)
                    api_two_group.add_species_difference_check_boxplot(boxfile_path, two_id)
                    api_two_group.add_species_difference_check_barplot(bar_path, two_id)
                    api_two_group.update_species_difference_check(two_id, stat_path, ci_path,'twogroup')
                    self.updata_status_api.add_meta_status(table_id=two_id, type_name='sg_species_difference_check')

                ## lefse
                if analysis_dir == "LEfSe":
                    api_lefse = self.api.stat_test
                    name = analysis_dir + "_"+self.table_name + group_dir
                    params = {
                        "end_level": 7,
                        'group_detail': self.group_detail_dict[group_dir],
                        'group_id': self.group_id_dict[group_dir],
                        "lda_filter": 2,
                        "normalization": 0,
                        "otu_id": str(self.result_otu_id),
                        "second_group_detail": "",
                        "second_group_id": "all",
                        "start_level": 3,
                        "strict": 1,
                        "submit_location": "species_lefse_analyse",
                        "task_type": "",
                    }
                    lefse_path = self.output_dir + "/DiffGroup/" + group_dir + "/" + analysis_dir + '/lefse_LDA.xls'
                    lefse_id =  api_lefse.create_species_difference_lefse(json.dumps(params),group_id=self.group_id_dict[group_dir],
                                                                          from_otu_table=str(self.result_otu_id),name=name)
                    api_lefse.add_species_difference_lefse_detail(file_path=lefse_path,
                                                                  table_id=lefse_id)
                    self.updata_status_api.add_meta_status(table_id=lefse_id, type_name='sg_species_difference_lefse')


    def export_env_analysis(self):
        for group_dir in os.listdir(self.output_dir + "/Env_analysis/"):

            ## MantelTest
            if os.path.exists(self.output_dir + "/Env_analysis/"+group_dir+"/MantelTest"):
                api_mantel = self.api.meta_species_env
                name = "MantelTest_"+self.table_name + group_dir
                params = {
                    "env_id": str(self.env_id),
                    "env_labs": self.env_names,
                    "env_method": "bray_curtis",
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    "level_id": 9,
                    "otu_id": str(self.result_otu_id),
                    "otu_method": "bray_curtis",
                    "submit_location": "beta_multi_analysis_results",
                    "task_type": "reportTask",
                }
                mantel_result = self.output_dir + "/Env_analysis/"+group_dir+"/MantelTest/Discompare/mantel_results.txt"
                dis_matrix = self.output_dir + "/Env_analysis/"+group_dir+"/MantelTest/Otudistance/bray_curtis_otu_taxon_otu.full.xls"
                fac_matrix = self.output_dir + "/Env_analysis/"+group_dir+"/MantelTest/Facdistance/factor_out.xls"
                mantel_id = api_mantel.add_mantel_table(9,self.result_otu_id,self.env_id,params=params, name=name)
                api_mantel.add_mantel_detail(mantel_result, mantel_id)
                api_mantel.add_mantel_matrix(dis_matrix, "species_matrix", mantel_id)
                api_mantel.add_mantel_matrix(fac_matrix, "env_matrix", mantel_id)
                self.updata_status_api.add_meta_status(table_id=mantel_id, type_name='sg_species_mantel_check')

            ## Rda
            if os.path.exists(self.output_dir + "/Env_analysis/" + group_dir + "/Rda"):
                api_multi = self.api.beta_multi_analysis
                params = {
                    "analysis_type": "rda_cca",
                    "env_id": str(self.env_id),
                    "env_labs": self.env_names,
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    "level_id": 9,
                    "otu_id": str(self.result_otu_id),
                    "submit_location": "beta_multi_analysis_rda_cca",
                    "task_type": "reportTask",
                }
                dir_path = self.output_dir + "/Env_analysis/" + group_dir
                #rda_id = api_multi.add_sg_otu(params, self.result_otu_id, name=name)
                rda_id = api_multi.add_beta_multi_analysis_result(dir_path, "rda_cca", main=True, env_id=self.env_id,
                                                        otu_id=self.result_otu_id, params=params, group_name=self.table_name+group_dir)
                self.updata_status_api.add_meta_status(table_id=rda_id, type_name='sg_beta_multi_analysis')

            ## heatmap_genus
            if os.path.exists(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus"):
                api_correlation = self.api.meta_species_env
                name = "SpearmanCorrelationGenus_"+self.table_name + group_dir
                params = {
                    "env_cluster": "average",
                    "env_id": str(self.env_id),
                    "env_labs": self.env_names,
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    "level_id": 7,
                    "method": "spearmanr",
                    "otu_id": str(self.result_otu_id),
                    "species_cluster": "average",
                    "submit_location": "beta_multi_analysis_pearson_correlation",
                    "task_type": "reportTask",
                    "top_species": "50",
                }
                corr_path = glob.glob(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus" + "/*correlation*")
                pvalue_path = glob.glob(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus" + "/*pvalue*")
                env_tree_path = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus/final_env_tree.tre"
                species_tree_path = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus/species_tree.tre"
                sorted_otu_file = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus/sorted_otu_file.xls"
                self.name_to_name = {}
                self.env_name = {}
                with open(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus/name_to_name.xls","r") as f:
                    for line in f:
                        line = line.strip().split("\t")
                        self.name_to_name[line[0]] = line[1]
                with open(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Genus/env_name.xls", "r") as ef:
                    self.env_name = eval(ef.readline())
                env_tree = ""
                new_env_tree = ""
                env_list = []
                species_list = []
                if os.path.exists(env_tree_path):
                    with open(env_tree_path, "r") as f:
                        env_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', env_tree)
                        #env_list = [self.env_name[i[1]] for i in raw_samp]
                        env_list = [i[1] for i in raw_samp]
                        new_env_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                if os.path.exists(species_tree_path):
                    with open(species_tree_path, "r") as f:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [self.name_to_name[i[1]] for i in raw_samp]
                        new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                corr_id1 = api_correlation.add_correlation(7,str(self.result_otu_id),str(self.env_id),
                                                           species_tree=new_species_tree,env_tree=new_env_tree, env_list=env_list,
                                                           species_list=species_list,params=params, name=name)
                api_correlation.add_correlation_detail(corr_path[0], "correlation", corr_id1,sorted_otu_file=sorted_otu_file)
                api_correlation.add_correlation_detail(pvalue_path[0], "pvalue", corr_id1, species_tree=new_species_tree,
                                                       env_tree=new_env_tree, env_list=env_list,
                                                       species_list=species_list,sorted_otu_file=sorted_otu_file)
                self.updata_status_api.add_meta_status(table_id=corr_id1, type_name='sg_species_env_correlation')

            if os.path.exists(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum"):
                api_correlation = self.api.meta_species_env
                name = "SpearmanCorrelationPhylum_"+self.table_name + group_dir
                params = {
                    "env_cluster": "average",
                    "env_id": str(self.env_id),
                    "env_labs": self.env_names,
                    'group_detail': self.group_detail_dict[group_dir],
                    'group_id': self.group_id_dict[group_dir],
                    "level_id": 3,
                    "method": "spearmanr",
                    "otu_id": str(self.result_otu_id),
                    "species_cluster": "average",
                    "submit_location": "beta_multi_analysis_pearson_correlation",
                    "task_type": "reportTask",
                    "top_species": "50",
                }
                corr_path = glob.glob(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum" + "/*correlation*")
                pvalue_path = glob.glob(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum" + "/*pvalue*")
                env_tree_path = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum/final_env_tree.tre"
                species_tree_path = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum/species_tree.tre"
                sorted_otu_file = self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum/sorted_otu_file.xls"
                self.name_to_name ={}
                self.env_name ={}
                with open(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum/name_to_name.xls", "r") as f:
                    for line in f:
                        line = line.strip().split("\t")
                        self.name_to_name[line[0]] = line[1]
                with open(self.output_dir + "/Env_analysis/" + group_dir + "/SpearmanCorrelation_Phylum/env_name.xls", "r") as ef:
                    self.env_name = eval(ef.readline())
                env_tree = ""
                new_env_tree = ""
                env_list = []
                species_list = []
                if os.path.exists(env_tree_path):
                    with open(env_tree_path, "r") as f:
                        env_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', env_tree)
                        #env_list = [self.env_name[i[1]] for i in raw_samp]
                        env_list = [i[1] for i in raw_samp]
                        new_env_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                if os.path.exists(species_tree_path):
                    with open(species_tree_path, "r") as f:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [self.name_to_name[i[1]] for i in raw_samp]
                        new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                corr_id2 = api_correlation.add_correlation(3,str(self.result_otu_id),str(self.env_id),
                                                           species_tree=new_species_tree,env_tree=new_env_tree, env_list=env_list,
                                                           species_list=species_list,params=params, name=name)
                api_correlation.add_correlation_detail(corr_path[0], "correlation", corr_id2,sorted_otu_file=sorted_otu_file)
                api_correlation.add_correlation_detail(pvalue_path[0], "pvalue", corr_id2, species_tree=new_species_tree,
                                                       env_tree=new_env_tree, env_list=env_list,
                                                       species_list=species_list,sorted_otu_file=sorted_otu_file)
                self.updata_status_api.add_meta_status(table_id=corr_id2,
                                                       type_name='sg_species_env_correlation')


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        task_info = self.api.api('task_info.task_info')
        task_info.add_task_info()
        self.logger.info("<<<<<<<<<<<<<<<<<<<")
        self.logger.info(self.option("info_path"))
        if self.option('envtable').is_set:  # 环境因子统一改名
            self.logger.info("change_env_name！")
            self.change_env_name()
        if self.option("info_path") == "" or not os.path.exists(self.option("info_path")):
            self.new_sample_extract.on("end", self.run_sample_rename)
        self.sample_rename.on("end", self.run_samplecheck)
        if self.option('pipeline') != "functional_gene":
            self.sample_check.on("end", self.run_otu)
            self.otu.on('end', self.check_otu_run)
        else:
            self.sample_check.on("end", self.run_framebot)
            self.framebot.on('end', self.run_otu)
            self.otu.on('end', self.check_otu_run)
        if self.option("if_fungene"):
            self.run_function_gene()
        elif self.option("info_path") == "" or not os.path.exists(self.option("info_path")):
            self.run_pre_sample_extract()
        else:
            self.run_sample_rename()
        super(MetaBaseWorkflow, self).run()

    def send_files(self):
        repaths = [
            [".", "", "基础分析结果文件夹", 0, "110001"],
            ["QC_stat", "", "样本数据统计文件目录", 0, "110004"],
            ["QC_stat/samples_info", "", "样本信息文件目录", 0, "110007"],  # add by hongdongxuan 20170323
            ["QC_stat/samples_info/samples_info.txt", "txt", "样本信息统计文件", 0, "110008"],
            ["QC_stat/base_info", "", "单个样本碱基质量统计目录", 0, "110011"],
            ["QC_stat/reads_len_info", "", "序列长度分布统计文件目录", 0, "110005"],
            ["QC_stat/valid_sequence.txt", "txt", "优化序列信息统计表", 0, "110009"],  # add by hongdongxuan 20170323
            ["QC_stat/序列长度分布图.pdf", "txt", "序列长度分布图.pdf", 0, ""],
            ["Otu", "dir", "OTU聚类结果文件目录", 0, "110026"],
            ["Tax_assign", "", "OTU对应物种分类文件目录", 0, "110002"],
            ["Tax_assign/seqs_tax_assignments.txt", "taxon.seq_taxon", "OTU序列物种分类文件", 0, "110003"],
            ["Tax_assign/blast_table.xls", "xls", "BLAST比对结果文件", 0, "110014"],   # add by zouxuan 20180126
            ["OtuTaxon_summary", "dir", "OTU物种分类综合统计目录", 0, "110013"],
            ["OtuTaxon_summary/otu_taxon.biom", "meta.otu.biom", "biom格式的OTU物种分类统计表", 0, "110024"],
            ["OtuTaxon_summary/otu_taxon.xls", "meta.otu.otu_table", "OTU物种分类统计表", 0, "110025"],
            ["OtuTaxon_summary/otu_summary.xls", "meta.otu.otu_table", "基于 OTU 数量的统计", 0, "110023"],
            ["OtuTaxon_summary/tax_summary_a", "meta.otu.tax_summary_dir", "各分类学水平样本序列数统计表", 0, "110019"],
            ["OtuTaxon_summary/tax_summary_r", "meta.otu.tax_summary_dir", "各分类学水平样本序列数相对丰度百分比统计表", 0, "110015"],  # add by zhouxuan 20161129
            ["OtuTaxon_summary/Otu_Rank_Abundance曲线图.pdf", "pdf", "Rank_Abundance曲线图", 0, ""],
            ["OtuTaxon_summary_depth", "dir", "抽平后OTU物种分类综合统计目录", 0, "110013"],
            ["OtuTaxon_summary_depth/otu_taxon.biom", "meta.otu.biom", "biom格式的OTU物种分类统计表", 0, "110024"],
            ["OtuTaxon_summary_depth/otu_taxon.xls", "meta.otu.otu_table", "OTU物种分类统计表", 0, "110025"],
            ["OtuTaxon_summary_depth/otu_summary.xls", "meta.otu.otu_table", "基于 OTU 数量的统计", 0, "110023"],
            ["OtuTaxon_summary_depth/tax_summary_a", "meta.otu.tax_summary_dir", "各分类学水平样本序列数统计表", 0, "110019"],
            ["OtuTaxon_summary_depth/tax_summary_r", "meta.otu.tax_summary_dir", "各分类学水平样本序列数相对丰度百分比统计表", 0, "110015"],
            ["OtuTaxon_summary_depth/Otu_Rank_Abundance曲线图.pdf", "pdf", "Rank_Abundance曲线图", 0, ""],
            ["Alpha_diversity", "", "Alpha diversity文件目录", 0, "110035"],
            ["Beta_diversity", "", "Beta diversity文件目录", 0, "110046"],
            ["Composition", "dir", "组成分析结果目录", 0, ""],
            ["CorrNetworkSpearmanGenus", "", "单因素相关性网络结果目录", 0, ""],
            ["DiffGroup", "dir", "物种差异分析结果目录", 0, ""],
            ["Env_analysis", "dir", "环境因子关联分析结果目录", 0, ""],
            ["Pan_Core", "dir", "Pan/core分析结果目录", 0, "110032"],     # add 3 lines by hongdongxuan 20170323
        ]
        for group in os.listdir(self.group_dir):
            group_name = group.rstrip(".group.txt")
            if os.path.exists(self.output_dir + "/CorrNetworkSpearmanGenus/" + group_name + "/newtable.txt"):
                os.remove(self.output_dir + "/CorrNetworkSpearmanGenus/" + group_name + "/newtable.txt")
            if os.path.exists(self.output_dir + "/CorrNetworkSpearmanGenus/"+ group_name +"/species_abundance.txt"):
                os.remove(self.output_dir + "/CorrNetworkSpearmanGenus/"+ group_name +"/species_abundance.txt")
            if os.path.exists(self.output_dir + "/Env_analysis/"+ group_name +"/SpearmanCorrelation_Genus/name_to_name.xls"):
                os.remove(self.output_dir + "/Env_analysis/"+ group_name +"/SpearmanCorrelation_Genus/name_to_name.xls")
            if os.path.exists(self.output_dir + "/Env_analysis/"+ group_name +"/SpearmanCorrelation_Genus/sorted_otu_file.xls"):
                os.remove(self.output_dir + "/Env_analysis/"+ group_name +"/SpearmanCorrelation_Genus/sorted_otu_file.xls")
            repaths += [
                ["Composition/%s/Circos" % group_name, "", "样本与物种关系分析结果目录", 0, ""],
                ["Alpha_diversity/%s/Estimators" % group_name, "", "多样性指数分析文件目录", 0, "110036"],
                ["Alpha_diversity/%s/Rarefaction" % group_name, "", "稀释曲线分析文件目录", 0, "110039"],  # add by guhaidong 20171018
                ["Alpha_diversity/%s/EstTTest" % group_name, "", "alpha多样性指数检验结果目录", 0, "110236"],
                ["Beta_diversity/%s/Anosim" % group_name, "", "ANOSIM&Adonis分析结果目录", 0, "110112"],
                ["Beta_diversity/%s/AnosimBox" % group_name, "", "ANOSIM&Adonis分析结果箱图数据目录", 0],
                ["Beta_diversity/%s/Dbrda" % group_name, "", "db_RDA分析结果目录", 0, "110106"],
                ["Beta_diversity/%s/Box" % group_name, "", "距离统计和统计检验分析结果目录", 0, "110238"],
                ["Beta_diversity/%s/Distance" % group_name, "", "距离矩阵计算结果目录", 0, "110047"],
                ["Beta_diversity/%s/Hcluster" % group_name, "", "层次聚类结果目录", 0, "110049"],
                ["Beta_diversity/%s/Nmds" % group_name, "", "NMDS分析结果目录", 0, ""],
                ["Beta_diversity/%s/Pca" % group_name, "", "PCA分析结果目录", 0, "110051"],
                ["Beta_diversity/%s/Pcoa" % group_name, "", "PCoA分析结果目录", 0, "110063"],
                ["Beta_diversity/%s/Plsda" % group_name, "", "PLS_DA分析结果目录", 0, "110092"],
                ["Beta_diversity/%s/Rda" % group_name, "", "RDA_CCA分析结果目录", 0, "110097"],
                ["Composition/%s/BarPie_Genus" % group_name, "", "Genus水平BarPie分析结果目录", 0, ""],
                ["Composition/%s/BarPie_Phylum" % group_name, "", "Phylum水平BarPie分析结果目录", 0, ""],
                ["Composition/%s/Circos" % group_name, "", "样本与物种关系分析结果目录", 0, ""],
                ["Composition/%s/Heatmap" % group_name, "", "群落Heatmap图分析结果目录", 0, ""],
                ["Composition/%s/Venn" % group_name, "", "物种venn图分析结果目录", 0, ""],
                ["CorrNetworkSpearmanGenus/%s/otu_association" % group_name, "", "物种相关性计算结果输出目录", 0, "110201"],
                ["CorrNetworkSpearmanGenus/%s/corr_network_calc" % group_name, "", "物种相关性网络分析结果输出目录", 0, "110194"],
                [r"CorrNetworkSpearmanGenus/%s/heatmap/" % group_name, "", "物种相关性热图结果输出目录", 0, ""],
                ["DiffGroup/%s/DiffStatMultiple_Genus" % group_name, "", "物种差异多组比较结果目录", 0, "110135"],
                ["DiffGroup/%s/DiffStatMultiple_Phylum" % group_name, "", "物种差异多组比较结果目录", 0, "110135"],
                ["DiffGroup/%s/DiffStatTwoGroup_Genus" % group_name, "", "物种差异两组比较结果目录", 0, "110139"],
                ["DiffGroup/%s/DiffStatTwoGroup_Phylum" % group_name, "", "物种差异两组比较结果目录", 0, "110139"],
                ["Env_analysis/%s/MantelTest" % group_name, "", "MantelTest分析结果目录", 0, "110157"],
                ["Env_analysis/%s/MantelTest/Discompare" % group_name, "", "Mantel_Test分析结果目录", 0, "110160"],
                ["Env_analysis/%s/MantelTest/Facdistance" % group_name, "", "环境因子矩阵结果目录", 0, "110158"],
                ["Env_analysis/%s/MantelTest/Otudistance" % group_name, "", "群落矩阵结果目录", 0, "110162"],
                ["Env_analysis/%s/partial" % group_name, "", "限制环境因子矩阵结果目录", 0, "110164"],
                ["Env_analysis/%s/Rda" % group_name, "", "RDA_CCA分析结果目录", 0, "110097"],
                ["Env_analysis/%s/SpearmanCorrelation_Genus" % group_name, "", "相关性Heatmap分析结果目录", 0, "110154"],
                ["Env_analysis/%s/SpearmanCorrelation_Phylum" % group_name, "", "相关性Heatmap分析结果目录", 0, "110154"],
                ["Pan_Core/%s/PanCoreOTU/pan曲线图.pdf" % group_name, "pdf", "pan曲线图", 0, ""],
                ["Pan_Core/%s/PanCoreOTU/core曲线图.pdf" % group_name, "pdf", "core曲线图", 0, ""],
                ["Beta_diversity/%s/AnosimBox/ANOSIM分析箱式图.pdf" % group_name, "pdf", "ANOSIM分析箱式图", 0, ""],
                ["Beta_diversity/%s/Hcluster/样本层级聚类分析图.pdf" % group_name, "pdf", "样本层级聚类分析结果图", 0],
                ["Beta_diversity/%s/Nmds/NMDS分析散点图.pdf" % group_name, "pdf", "样本NMDS分析散点图", 0],
                ["Beta_diversity/%s/Pca/PCA分析散点图.pdf" % group_name, "pdf", "PCA分析散点图", 0, ""],
                ["Beta_diversity/%s/Pca/PCA分析箱式图.pdf" % group_name, "pdf", "PCA分析箱式图", 0, ""],
                ["Beta_diversity/%s/Pcoa/PCoA分析散点图.pdf" % group_name, "pdf", "PCoA分析散点图", 0, ""],
                ["Beta_diversity/%s/Pcoa/PCoA分析箱式图.pdf" % group_name, "pdf", "PCoA分析箱式图", 0, ""],
                ["Composition/%s/BarPie_Genus/群落柱形图.pdf" % group_name, "pdf", "群落柱形图", 0, ""],
                ["Composition/%s/BarPie_Phylum/群落柱形图.pdf" % group_name, "pdf", "群落柱形图", 0, ""],
                ["Composition/%s/Circos/Circos图.pdf" % group_name, "pdf", "样本与物种关系Circos图", 0, ""],
                ["Composition/%s/Circos/三元相图.pdf" % group_name, "pdf", "三元相图", 0, ""],
                ["Composition/%s/Heatmap/群落Heatmap图.pdf" % group_name, "pdf", "物种群落Heatmap图", 0, ""],
                ["Composition/%s/Venn/Venn图.pdf" % group_name, "pdf", "Venn图", 0],
                ["CorrNetworkSpearmanGenus/%s/heatmap/物种相关性Heatmap图.pdf" % group_name, "pdf", "物种相关性Heatmap图", 0, ""],
                ["CorrNetworkSpearmanGenus/%s/corr_network_calc/单因素相关性网络图.pdf" % group_name, "pdf", "物种与物种间的相关性网络图", 0, ""],
                ["DiffGroup/%s/*/两组比较差异检验柱形图.pdf" % group_name, "pdf", "两组比较差异检验柱形图", 0, ""],
                ["DiffGroup/%s/*/多组比较差异检验柱形图.pdf" % group_name, "pdf", "多组比较差异检验柱形图", 0, ""],
                ["DiffGroup/%s/LEfSe/LDA判别结果图.pdf" % group_name, "pdf", "LDA判别结果图", 0, ""],
                ["DiffGroup/%s/LEfSe/LEfSe多级物种层级树图.pdf" % group_name, "pdf", "LEfSe多级物种层级树图", 0, ""],
                ["Env_analysis/%s/Rda/RDA_CCA分析结果图.pdf" % group_name, "pdf", "RDA/CCA分析结果图", 0, ""],
                ["Env_analysis/%s/SpearmanCorrelation_*/相关性Heatmap图.pdf" % group_name, "pdf", "物种与环境因子相关性Heatmap图", 0, ""]
            ]
        regexps = [
            [r"Pan_Core/.+/PanCoreOTU/core\.richness\.xls", "xls", "core 表格", 0, "110034"],
            [r"Pan_Core/.+/PanCoreOTU/pan\.richness\.xls", "xls", "pan 表格", 0, "110033"],
            ["QC_stat/base_info/.*\.fastq\.fastxstat\.txt", "", "单个样本碱基质量统计文件", 0, "110012"],
            [r"QC_stat/reads_len_info/step_\d+\.reads_len_info\.txt", "", "序列长度分布统计文件", 0, "110006"],
            [r"Alpha_diversity/.+/Estimators/otu\..*\.summary", "", "单个样本多样性指数表", 0, "110037"],  # add by guhaidong 20171018
            [r'Beta_diversity/.+/Distance/%s.*\.xls$' % self.option('dis_method'), 'meta.beta_diversity.distance_matrix', '样本距离矩阵文件', 0, "110048"],
            [r'Beta_diversity/.+/Rda/.+_importance\.xls$', 'xls', '主成分变化解释度表', 0, "110104"],
            [r'Beta_diversity/.+/Rda/.+_sites\.xls$', 'xls', '样本坐标表', 0, "110100"],
            [r'Beta_diversity/.+/Rda/.+_species\.xls$', 'xls', '物种坐标表', 0, "110103"],
            [r'Beta_diversity/.+/Rda/.+_biplot\.xls$', 'xls', '数量型环境因子坐标表', 0, "110102"],
            [r'Beta_diversity/.+/Rda/.+_centroids\.xls$', 'xls', '哑变量环境因子坐标表', 0, "110105"],
            [r"Otu/otu_reps.fasta", "sequence.fasta", "OTU代表序列", 0, "110029"],
            [r"Otu/otu_seqids.txt", "txt", "每个OTU中包含的序列编号列表", 0, "110027"],
            [r"Otu/otu_table.biom", 'meta.otu.biom', "OTU表对应的Biom文件", 0, "110028"],
            [r"Otu/otu_table.xls", "meta.otu.otu_table", "各样本OTU中序列数统计表", 0, "110030"],
            [r"Otu/otu_phylo.tre", "graph.newick_tree", "OTU代表序列进化树", 0, "110031"],
            [r"QC_stat/base_info/.*\.fastq\.fastxstat\.txt", "txt", "单个样本碱基质量统计文件", 0, "110012"],
            [r"QC_stat/reads_len_info/step_\d+\.reads_len_info\.txt", "txt", "序列长度分布统计文件", 0, "110006"],
            [r"OtuTaxon_summary/tax_summary_a/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件(absolute)", 0, "110021"],
            [r"OtuTaxon_summary/tax_summary_a/.+\.xls$", "xls", "单级物种分类统计表(absolute)", 0, "110022"],
            [r"OtuTaxon_summary/tax_summary_a/.+\.full\.xls$", "xls", "多级物种分类统计表(absolute)", 0, "110020"],
            [r"OtuTaxon_summary/tax_summary_r/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件", 0, "110018"],  # add by zhouxuan (3 line) 20161129
            [r"OtuTaxon_summary/tax_summary_r/.+\.xls$", "xls", "单级物种分类统计表", 0, "110017"],
            [r"OtuTaxon_summary/tax_summary_r/.+\.full\.xls$", "xls", "多级物种分类统计表", 0, "110016"],
            [r"OtuTaxon_summary_depth/tax_summary_a/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件(absolute)", 0, "110021"],
            [r"OtuTaxon_summary_depth/tax_summary_a/.+\.xls$", "xls", "单级物种分类统计表(absolute)", 0, "110022"],
            [r"OtuTaxon_summary_depth/tax_summary_a/.+\.full\.xls$", "xls", "多级物种分类统计表(absolute)", 0, "110020"],
            [r"OtuTaxon_summary_depth/tax_summary_r/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件", 0, "110018"],
            [r"OtuTaxon_summary_depth/tax_summary_r/.+\.xls$", "xls", "单级物种分类统计表", 0, "110017"],
            [r"OtuTaxon_summary_depth/tax_summary_r/.+\.full\.xls$", "xls", "多级物种分类统计表", 0, "110016"],
            [r"Alpha_diversity/.+/Estimators/estimators.xls", "xls", "Alpha多样性指数表", 0, "110038"],
            [r"Alpha_diversity/.+/Rarefaction/.+/*\.xls", "", "alpha多样性指数T检验结果表", 0, "110237"],
            [r"Beta_diversity/.+/Anosim/anosim_results.txt", "txt", "anosim分析结果", 0, "110113"],
            [r"Beta_diversity/.+/Anosim/adonis_results.txt", "txt", "adonis分析结果", 0, "110114"],
            [r"Beta_diversity/.+/Anosim/format_results.xls", "xls", "anosim&adonis综合统计表", 0, "110116"],
            [r"Beta_diversity/.+/Dbrda/db_rda_sites.xls", "xls", "db_rda样本坐标表", 0, "110109"],
            [r"Beta_diversity/.+/Dbrda/db_rda_species.xls", "xls", "db_rda物种坐标表", 0, "110110"],
            [r"Beta_diversity/.+/Dbrda/db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表", 0, "110111"],
            [r"Beta_diversity/.+/Dbrda/db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表", 0, "110108"],
            [r"Beta_diversity/.+/Box/Stats.xls", "xls", "分组统计检验结果", 0, "110239"],
            [r"Beta_diversity/.+/Box/Distances.xls", "xls", "组内组间距离值统计结果", 0, "110240"],
            [r"Beta_diversity/.+/Hcluster/hcluster.tre", "graph.newick_tree", "层次聚类树结果表", 0, "110050"],
            [r"Beta_diversity/.+/Nmds/nmds_sites.xls", "xls", "样本各维度坐标", 0, "110061"],
            [r"Beta_diversity/.+/Nmds/nmds_stress.xls", "xls", "样本特征拟合度值", 0, "110062"],
            [r"Beta_diversity/.+/Pca/pca_importance.xls", "xls", "主成分解释度表", 0, "110055"],
            [r"Beta_diversity/.+/Pca/pca_rotation.xls", "xls", "PCA主成分贡献度表", 0, "110053"],
            [r"Beta_diversity/.+/Pca/pca_rotation_all.xls", "xls", "PCA全部主成分贡献度表", 0, "110052"],
            [r"Beta_diversity/.+/Pca/pca_sites.xls", "xls", "样本各成分轴坐标", 0, "110054"],
            [r"Beta_diversity/.+/Pca/pca_envfit_factor_scores.xls", "xls", "哑变量环境因子坐标表", 0, "110056"],
            [r"Beta_diversity/.+/Pca/pca_envfit_factor.xls", "xls", "哑变量环境因子表", 0, "110057"],
            [r"Beta_diversity/.+/Pca/pca_envfit_vector_scores.xls", "xls", "数量型环境因子坐标表", 0, "110059"],
            [r"Beta_diversity/.+/Pca/pca_envfit_vector.xls", "xls", "数量型环境因子表", 0, "110058"],
            [r"Beta_diversity/.+/Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, "110066"],
            [r"Beta_diversity/.+/Pcoa/pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比", 0, "110064"],
            [r"Beta_diversity/.+/Pcoa/pcoa_sites.xls", "xls", "样本坐标表", 0, "110065"],
            [r'Beta_diversity/.+/Rda/dca.xls', 'xls', 'DCA分析结果', 0, "110098"],
            [r"Beta_diversity/.+/Plsda/plsda_sites.xls", "xls", "样本坐标表", 0, "110096"],
            [r"Beta_diversity/.+/Plsda/plsda_rotation.xls", "xls", "物种主成分贡献度表", 0, "110095"],
            [r"Beta_diversity/.+/Plsda/plsda_importance.xls", "xls", "主成分组别特征值表", 0, "110093"],
            [r"Beta_diversity/.+/Plsda/plsda_importancepre.xls", "xls", "主成分解释度表", 0, "110094"],
            [r"Composition/.+/BarPie_Genus/taxa.table.xls", "xls", "各样本物种丰度结果表", 0, "110076"],
            [r"Composition/.+/BarPie_Genus/taxa.precents.table.xls", "xls", "各样本物种相对丰度结果表", 0, "110075"],
            [r"Composition/.+/BarPie_Phylum/taxa.table.xls", "xls", "各样本物种丰度结果表", 0, "110076"],
            [r"Composition/.+/BarPie_Phylum/taxa.precents.table.xls", "xls", "各样本物种相对丰度结果表", 0, "110075"],
            [r"Composition/.+/Circos/taxa.table.xls", "xls", "各样本物种丰度结果表", 0, "110076"],
            [r"Composition/.+/Circos/taxa.precents.table.xls", "xls", "各样本物种相对丰度结果表", 0, "110075"],
            [r"Composition/.+/Heatmap/heatmap.taxa.table.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, "110079"],
            [r"Composition/.+/Heatmap/heatmap.taxa.relative.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, ""],
            [r"Composition/.+/Heatmap/sample_hcluster.tre", "tre", "样本聚类树", 0, "110078"],
            [r"Composition/.+/Heatmap/species_hcluster.tre", "tre", "物种聚类树", 0, "110080"],
            [r"Composition/.+/Venn/venn_table.xls", "xls", "Venn表格", 0, "110073"],
            [r"CorrNetworkSpearmanGenus/.+/otu_association/shared.txt", "txt", "shared文件", 0, "110203"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_attributes.txt", "txt", "网络的单值属性表", 0, "110197"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件", 0, "110196"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_centrality.txt", "txt", "网络节点的中心系数表", 0, "110195"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_clustering.txt", "txt", "网络节点的聚类系数表", 0, "110200"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_degree_distribution.txt", "txt", "网络节点的度分布表", 0, "110199"],
            [r"CorrNetworkSpearmanGenus/.+/corr_network_calc/corr_network_node_degree.txt", "txt", "网络节点的度统计表", 0, "110198"],
            [r"DiffGroup/.+/.+/.+_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, "110137"],
            [r"DiffGroup/.+/.+/.+(-)*\.xls", "xls", "组间差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值", 0, "110136"],
            [r"DiffGroup/*/*/.+_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "110138"],
            [r"Env_analysis/.+/Discompare/mantel_results.txt", "txt", "Mantel_Test分析结果表", 0, "110161"],
            [r"Env_analysis/.+/Facdistance/factor_out.xls", "xls", "环境因子矩阵结果表", 0, "110159"],
            [r"Env_analysis/.+/Otudistance/%s_otu_taxon_otu.full.xls" % self.option("otu_method"), "xls", "群落矩阵结果表", 0, "110163"],
            [r"Env_analysis/.+/partial/factor_out.xls", "xls", "限制环境因子矩阵结果表", 0, "110165"],
            [r'Env_analysis/.+/Rda/dca.xls', 'xls', 'DCA分析结果', 0, "110098"],
            [r"Env_analysis/.+/SpearmanCorrelation_*/pearsons_correlation_at_otu_level.xls", "xls", "相关性系数表", 0, "110155"],
            [r"Env_analysis/.+/SpearmanCorrelation_*/pearsons_pvalue_at_otu_level.xls", "xls", "相关性P值", 0, "110156"],
            [r'Env_analysis/.+/Rda/.*_sites\.xls$', 'xls', '样本坐标表', 0, "110100"],
            [r'Env_analysis/.+/Rda/.*_species\.xls$', 'xls', '物种坐标表', 0, "110103"],
            [r'Env_analysis/.+/Rda/.*_biplot\.xls$', 'xls', '数量型环境因子坐标表', 0, "110102"],
            [r'Env_analysis/.+/Rda/.*_centroids\.xls$', 'xls', '哑变量环境因子坐标表', 0, "110105"],
            [r"Env_analysis/.+/.+/pearsons_correlation.xls", "xls", "相关性系数表", 0, "110155"],
            [r"Env_analysis/.+/.+/pearsons_pvalue.xls", "xls", "相关性P值", 0, "110156"],
        ]
        for group in os.listdir(self.group_dir):
            group_name = group.rstrip(".group.txt")
            for i in self.option("estimate_indices").split(","):
                repaths.append(["./Alpha_diversity/{}/Estimators/{}指数柱形图.pdf".format(group_name,i), "pdf", "各样本的{}指数柱形图".format(i), 0, ""])
            for i in self.option("rarefy_indices").split(","):
                dir_code_list = {
                    "shannon": "110040",
                    "sobs": "110042"
                }
                file_code_list = {
                    "shannon": "110041",
                    "sobs": "110043"
                }
                if i == "sobs":  # modified by hongdongxuan 20170324
                    # repaths.append(["./rarefaction", "文件夹", "{}指数结果输出目录".format(i)])
                    repaths.append(["./Alpha_diversity/{}/Rarefaction/{}".format(group_name,i), "文件夹", "{}指数结果输出目录".format(i), 0])
                    # regexps.append([r".*rarefaction\.xls", "xls", "{}指数的simpleID的稀释性曲线表".format(i)])
                    regexps.append([r".*rarefaction\.xls", "xls", "每个样本的{}指数稀释性曲线表".format(group_name,i), 0, file_code_list[i]])
                    repaths.append(["./Alpha_diversity/{}/Rarefaction/各样本的{}指数稀释性曲线图.pdf".format(group_name, i), "pdf",
                                    "各样本的{}指数稀释性曲线图".format(i), 0, ""])
                elif i == "shannon":
                    repaths.append(["./Alpha_diversity/{}/Rarefaction/{}".format(group_name,i), "文件夹", "{}指数结果输出目录".format(i), 0,
                                    dir_code_list[i]])
                    regexps.append([r".*{}\.xls".format(i), "xls", "每个样本的{}指数稀释性曲线表".format(i), 0, file_code_list[i]])
                    repaths.append(["./Alpha_diversity/{}/Rarefaction/各样本的{}指数稀释性曲线图.pdf".format(group_name, i), "pdf",
                                    "各样本的{}指数稀释性曲线图".format(i), 0, ""])
                else:
                    if i != "":
                        # repaths.append(["./{}".format(i), "文件夹", "{}指数结果输出目录".format(i)])
                        # repaths.append(["./Alpha_diversity/Rarefaction/{}".format(i), "文件夹", "{}指数结果输出目录".format(i)])
                        repaths.append(
                            ["./Alpha_diversity/{}/Rarefaction/{}".format(group_name,i), "文件夹", "稀释曲线结果输出目录", 0, "110044"])
                        repaths.append([r"./Alpha_diversity/{}/Rarefaction/各样本的{}指数稀释性曲线图.pdf".format(group_name, i), "pdf",
                                        "各样本的{}指数稀释性曲线图".format(i), 0, ""])
                        repaths.append(
                            # [r".*{}\.xls".format(i), "xls", "{}指数的simpleID的稀释性曲线表".format(i)])
                            # [r".*{}\.xls".format(i), "xls", "每个样本的{}指数稀释性曲线表".format(i)])
                            [r".*{}\.xls".format(i), "xls", "每个样本的稀释性曲线表", 0, "110045"])
                repaths.append(
                    ["./Alpha_diversity/{}/EstTTest/{}指数检验差异检验柱形图.pdf".format(group_name,i), "pdf", "{}指数的组间差异检验结果".format(i), 0,
                     ""])
                repaths.append(
                    ["./Alpha_diversity/{}/EstTTest/{}指数检验差异检验箱式图.pdf".format(group_name,i), "pdf", "{}指数的组间差异检验结果".format(i), 0,
                     ""])
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))

    def end(self):
        if self.option("save_pdf"):
            os.system("cp -r {}/* {}/".format(self.figsave.output_dir, self.output_dir))
        is_pollu=check_pollution_pip(self.output_dir+'/OtuTaxon_summary/tax_summary_r',self._sheet.id,self.work_dir,is_sanger=self._sheet.UPDATE_STATUS_API)  #做预警物种检查 201909
        #is pollu : 0,1 表示无预警，有预警.对应api 的 1,2
        self.logger.info('warning_pollu: %s'% (is_pollu+1))
        self.add_task_option('warning_pollu', is_pollu+1) #20200708
        self.logger.info(self.task_option_data)
        #self.step.update()
        self.send_files()
        super(MetaBaseWorkflow, self).end()

    def get_origin_name(self):
        """
        从data_json中获取in_fastq的s3_name
        :return:
        """
        s3_origin = {}
        data_json = os.path.join(self.work_dir, "data.json")
        with open(data_json, 'r') as f:
            js = f.read()
            json_dict = json.loads(js)
        in_fastq = str(json_dict['options']['in_fastq']).split("||")[-1]
        true_file_name = in_fastq.split("||")[-1].split("{")[0]
        new_file_name = in_fastq.split("||")[-1].split("{")[1].rstrip("}")
        if new_file_name not in s3_origin:
            s3_origin[new_file_name] = true_file_name
        # new_file_path = os.path.join(self.work_dir, "remote_input/in_fastq/", new_file_name)
        # origin_path = os.path.join(self.work_dir, true_file_name)
        # os.link(new_file_path, origin_path)
        # self.option("in_fastq").set_path(origin_path)
        return s3_origin

    def get_dir_name(self):
        """
        从mapping_file文件中获取s3_name与改名后的文件名称的对应关系
        :return:
        """
        s3_origin = {}
        mapping_file = os.path.join(self.work_dir, "remote_input/in_fastq/mapping_file.txt")
        with open(mapping_file, 'r') as mm:
            js = mm.read()
            json_dict = json.loads(js)
        fastq_list = json_dict['in_fastq']
        for fastq_dict in fastq_list:
            new_file_name = fastq_dict['alias']
            true_file_name = fastq_dict['file_path']
            if new_file_name not in s3_origin:
                s3_origin[new_file_name] = true_file_name
        return s3_origin

    def dashrepl(self, matchobj):
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        return self.env_name[matchobj.groups()[0]]