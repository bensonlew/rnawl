# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import re
import os
import json
import shutil
from bson.objectid import ObjectId
import gevent
from Bio import SeqIO
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.search_polution_by_list import check_pollution_pip
from mbio.packages.metaasv.common_function import link_dir,check_file,link_file,get_group_from_table
from mbio.packages.meta.delete_mongo import DeleteDemoMongo
import functools

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class MetaAsvWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        metaasv Qiime2工作流
        """
        self._sheet = wsheet_object
        self.DATABASE = ['unite7.2/its_fungi',"unite8.0/its_fungi",
                         'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA', 'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                         'maarjam081/AM', 'Human_HOMD', 'Human_HPB', 'Protist_PR2_v4.5',"Human_HOMD_v15.2",
                         'silva132/16s_archaea', 'silva132/16s_bacteria','silva132/18s_eukaryota', 'silva132/16s',
                         'silva138/16s_archaea', 'silva138/16s_bacteria','silva138/18s_eukaryota', 'silva138/16s',
                         'greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria',
                         'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                         'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi','nt_v20210917',
                         'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'nt', 'nt/16s','nt/18s','nt/its_fungi', 'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s','nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi', "nt_v20200604",
                         'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012','fgr/amoA_comammox_202012', 'fgr/nosZ_202012', 'fgr/nosZ_atypical_1_202012',
                         'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012','fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012', 'fgr/mmoX_202012']
        self.DISTANCE = ["abund_jaccard","binary_euclidean","binary_hamming","bray_curtis","euclidean","hellinger","unifrac","unweighted_unifrac","weighted_unifrac", ]
        self.ESTIMATORS = ["ace", "chao", "shannon", "simpson", "coverage", "sobs", "shannoneven", "simpsoneven"]
        super(MetaAsvWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fastq_type', 'type': 'string'},  # 文件类型
            {'name': 'fastq_file', 'type': 'infile', 'format': 'sequence.fastq,sequence.fastq_dir'},  # 输入的fastq文件或fastq文件夹
            {"name": "raw_sequence", "type": "infile", "format": "metaasv.raw_sequence_txt"},## 原始序列信息表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},##group分组表
            {'name': 'asv_upload', 'type': 'infile', 'format': 'sequence.fasta'},##流程2参数
            {'name': 'asv_abundance', 'type': 'infile', 'format': 'meta.otu.otu_table'},##流程2参数
            {'name': 'denoise_method', 'type': 'string'},  # 降噪方法选择
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'database_type', 'type': 'string'},  # 数据库类型
            {'name': 'anno_method', 'type': 'string'},  # 注释方法选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'identity', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'coverage', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            {"name": "estimate_indices", "type": "string", "default": "ace,chao,shannon,simpson,coverage,sobs,shannoneven,simpsoneven"},## 多样性指数类型
            {"name": "rarefy_indices", "type": "string", "default": "ace,chao,shannon,simpson,coverage,sobs,shannoneven,simpsoneven"},  # 稀释曲线指数类型
            {"name": "beta_analysis", "type": "string", "default": "pca,hcluster,pcoa,nmds"},##排序回归分析
            {"name": "composition_type", "type": "string", "default": "barpie,heatmap"},##组成分析
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},## 距离算法
            {"name": "query_id", "type": "string", "default": ""},## 文件的编号，如果工作流直接运行则必须要传此字段，否则前端无法正常显示样本检测内容
            {"name": "truc_len", "type": "int", "default": 0},## 样本的最小样本序列数
            {"name": "max_ee", "type": "int", "default": 2}, ## DADA2所用的参数
            {"name": "trunc_q", "type": "int", "default": 0}, ##判断打断序列的大小，如果为0则不打断，即不质控
            {"name": "min_size", "type": "int", "default": 1},
            {"name": "min_reads", "type": "int", "default": 1},
            {"name": "ref_acid", "type": "infile", "format": "sequence.fasta"},  # 参考氨基酸fasta文件
            {"name": "pipeline", "type": "string", "default": ""},  # 选择常规流程还是功能基因流程
            {'name': 'identity', 'type': 'float', 'default': 0.4},  # FrameBot移码校正的序列一致性阈值，范围0-1.
            {"name": "fungene_database", "type": "string"},  # 功能基因流程功能基因名称
            {"name": "acid_length", "type": "int", "default": 80},  # 氨基酸长度阈值
            {"name": "seq_identity", "type": "float", "default": 0.4},
            {"name": "seq_style", "type": "string", "default": "nucl"},
            {'name': 'fungene_anno_method', 'type': 'string'},  # 注释方法选择
            {'name': 'fungene_identity', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'fungene_coverage', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'fungene_confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "tax_database", "type": "string"},  # 物种分类数据库,
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sample_rename = self.add_tool("metaasv.sample_rename")
        self.sample_check = self.add_tool("metaasv.sample_check")
        self.split_sample = self.add_module("metaasv.sample_split")
        self.denoise = self.add_module("metaasv.qiime2_denoise")##降噪模块
        self.framebot = self.add_module("metaasv.framebot")
        self.phylo = self.add_tool("phylo.phylo_tree")
        self.anno = self.add_module("metaasv.asv_annotation")
        # self.alpha = self.add_module("meta.alpha_diversity.alpha_diversity")
        self.alpha = self.add_tool('meta.alpha_diversity.estimators')
        self.rarefaction = self.add_module("meta.alpha_diversity.rarefaction")
        self.beta = self.add_module("metaasv.beta_diversity")
        self.composition = self.add_module("metaasv.composition_analysis")
        self.pan_core = self.add_tool("metaasv.pan_core")
        self.rank_abundance = self.add_tool('metaasv.rank_abundance')
        self.step.add_steps("qc_stat", "denoise", "framebot","taxassign", "basic_analysis")
        self.asv_phylotree = []
        self.basic_analysis = []
        self.spname_spid = dict()
        self.level_dict = {'Domain': 1, 'Kingdom': 2, 'Phylum': 3, 'Class': 4, 'Order': 5, 'Family': 6, 'Genus': 7, 'Species': 8, 'asv': 9}
        self.updata_status_api = self.api.meta_update_status
        self.info_path = ""
        if self.option("asv_upload").is_set and self.option("asv_abundance").is_set:
            self.pipline = "pipeline2"
        else:
            self.pipline = "pipeline1"

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'metaasv')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        参数二次检查
        """
        self.beta_analysis = self.option("beta_analysis").split(",")
        self.estimate_indices = self.option("estimate_indices").split(",")
        if self.option("identity") < 0 or self.option("identity") > 1:
            raise OptionError("identity值必须在0-1范围内.")
        if self.option("coverage") < 0 or self.option("coverage") > 1:
            raise OptionError("coverage值必须在0-1范围内.")
        if self.option("confidence") < 0 or self.option("confidence") > 1:
            raise OptionError("confidence值必须在0-1范围内.")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件")
        elif self.option('pipeline') == "functional_gene":
            if  self.option('fungene_database') == "custom_mode":
                if  not self.option("ref_acid").is_set or not self.option("ref_taxon").is_set:
                    raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件")
                if self.option("seq_style") == "nucl":
                    if not self.option("ref_fasta").is_set:
                        raise OptionError("数据库自定义模式下ASV数据类型为核苷酸必须设置参考核酸fasta序列和参考taxon文件")
            else:
                if self.option("fungene_database") not in ['fgr/amoA_archaea_202012','fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012',
                                               'fgr/amoA_comammox_202012', 'fgr/nosZ_202012','fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012',
                                               'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012',
                                               'fgr/pmoA_202012', 'fgr/mmoX_202012']:
                    raise OptionError("数据库%s不被支持", variables=(self.option("fungene_database")))
        else:
            if self.option("database") not in self.DATABASE:
                raise OptionError("数据库%s不被支持", variables=(self.option("database")))
            if self.option("database") in ['nt', 'nt_v20200604']:
                if self.option("anno_method") not in ["multi_blast"]:
                    raise OptionError("数据库为nt时只支持Multi_Blast方法")
        if len(self.beta_analysis) == 1 and "pca" not in self.beta_analysis:
            if not self.option("dis_method"):
                raise OptionError("必须设置距离算法!")
        for indices in self.estimate_indices:
            if indices not in self.ESTIMATORS:
                raise OptionError("多样性指数%s不正确，请检查多样性指数的配置", variables=(indices))
        if self.option("dis_method") not in self.DISTANCE:
            raise OptionError("距离方法：%s不正确，请重新选择距离方法进行计算", variables=(self.option("dis_method")))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        """
        运行和设置运行逻辑
        :return:
        """
        task_info = self.api.api('task_info.metaasv_task_info')
        task_info.add_task_info()
        task_info.update_sg_task(self._sheet.id, self.option("database_type"), self.option("database"), self.option("pipeline"))
        self.logger.info("<<<<<<<<<<<<<<<<<<<")
        if self.pipline in ["pipeline1"]:
            if self.option("fastq_type") in ["fastq_dir"]:
                if self.option("pipeline") == "functional_gene":
                    self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
                    self.seq_extract.on("end", self.run_sample_rename)
                    self.seq_extract.on("end", self.split_sample_fastq)
                    self.on_rely([self.split_sample, self.sample_rename], self.run_samplecheck)
                    self.sample_check.on("end", self.run_denoise)
                    self.denoise.on('end', self.run_framebot)
                    self.framebot.on('end', self.check_otu_run)
                    self.is_zip = check_file('dir', self.option("fastq_file").prop["path"])
                    if self.is_zip:
                        self.sequence = self.add_module('metaasv.reads_unzip')
                        self.sequence.on("end", self.run_seq_extract)
                        self.run_sequence(self.option("fastq_file").prop['path'])
                    else:
                        self.run_seq_extract()
                else:
                    self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
                    self.seq_extract.on("end", self.run_sample_rename)
                    self.seq_extract.on("end", self.split_sample_fastq)
                    self.on_rely([self.split_sample, self.sample_rename], self.run_samplecheck)
                    self.sample_check.on("end", self.run_denoise)
                    self.denoise.on('end', self.check_otu_run)
                    self.is_zip = check_file('dir', self.option("fastq_file").prop["path"])
                    if self.is_zip:
                        self.sequence = self.add_module('metaasv.reads_unzip')
                        self.sequence.on("end", self.run_seq_extract)
                        self.run_sequence(self.option("fastq_file").prop['path'])
                    else:
                        self.run_seq_extract()

            elif self.option("fastq_type") in ["fastq"]:
                if self.option("pipeline") == "functional_gene":
                    self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
                    self.seq_extract.on("end", self.run_sample_rename)
                    self.seq_extract.on("end", self.split_sample_fastq)
                    self.on_rely([self.split_sample, self.sample_rename], self.run_samplecheck)
                    self.sample_check.on("end", self.run_denoise)
                    self.denoise.on('end', self.run_framebot)
                    self.framebot.on('end', self.check_otu_run)
                    self.is_zip = check_file('file', self.option("fastq_file").prop["path"])
                    if self.is_zip:
                        self.sequence = self.add_tool('sequence.fastq_ungz')
                        self.sequence.on("end", self.run_seq_extract)
                        self.run_file_sequence(self.option("fastq_file").prop['path'])
                    else:
                        self.run_seq_extract()
                else:
                    self.seq_extract = self.add_module("meta.sample_extract.sample_extract")
                    self.seq_extract.on("end", self.run_sample_rename)
                    self.seq_extract.on("end", self.split_sample_fastq)
                    self.on_rely([self.split_sample, self.sample_rename], self.run_samplecheck)
                    self.sample_check.on("end", self.run_denoise)
                    self.denoise.on('end', self.check_otu_run)
                    self.is_zip = check_file('file', self.option("fastq_file").prop["path"])
                    if self.is_zip:
                        self.sequence = self.add_tool('sequence.fastq_ungz')
                        self.sequence.on("end", self.run_seq_extract)
                        self.run_file_sequence(self.option("fastq_file").prop['path'])
                    else:
                        self.run_seq_extract()
        else:
            self.check_otu_run()
        super(MetaAsvWorkflow, self).run()

    def check_otu_run(self):
        """
        检查ASV数目和样本数目，确定下游分析内容
        """
        self.update_info = ""
        self.count_samples = 0  # 样本数量是否大于等于2
        if self.pipline in ["pipeline1"]:
            counts = SeqIO.index(os.path.join(self.denoise.output_dir, 'ASV_reps.fasta'), 'fasta')
            with open(self.sample_rename.output_dir + "/info_path.xls") as r:
                self.count_samples = len(r.readlines()) - 1
        else:
            counts = SeqIO.index(self.option("asv_upload").prop['path'], 'fasta')
            with open(self.option("asv_abundance").prop['path'], 'r') as r:
                self.count_samples = len(r.readline()) - 1
        if counts > 3:    #
            self.count_otus = True
        else:
            self.count_otus = False
        if self.count_samples < 3:  # 少于3个样本,不进行beta多样性相关分析
            self.update_info += "Sample size too small: SKIP Beta diversity analysis!"
            self.option('beta_analysis', '')
        if not self.count_otus:
            genus_path = os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a","asv_taxon_Genus.full.xls")
            with open(genus_path, 'r') as r:
                genus_counts = len(r.readline()) - 1
            if genus_counts <= 3:
                if 'pca' in self.option('beta_analysis'): # "OTU数量过少，不做PCA分析；
                    self.option('beta_analysis', self.option('beta_analysis').replace('pca', '').strip(',').replace(',,', ','))
                    self.update_info += "ASV size too small: SKIP PCA analysis!"
                if "unifrac" in self.option('dis_method'):
                    self.option('dis_method', "bray_curtis")##"OTU数量过少，不使用unifrac类型距离算法，采用默认bray_curtis算法；"
                    self.update_info += "ASV size too small: UniFrac distance method replaced by bray_curtis!"

        if self.count_otus and self.count_samples > 1:
            self.on_rely([self.phylo, self.anno], self.run_basic_analysis)
            self.run_taxon_format()
            self.run_phylotree()
        elif self.count_otus and self.count_samples == 1:
            ##与上面的区别是：样本数太少，不做pan_core分析
            self.on_rely([self.phylo, self.anno], self.run_basic_analysis2)
            self.run_taxon_format()
            self.run_phylotree()
        else:
            #"OTU数量过少，不进行物种进化树分析，将不会生成物种进化树；"
            self.update_info += "OTU size too small: SKIP Phylogenetic analysis!"
            self.anno.on('end', self.run_alpha)
            line_counts = 0 ## fix by qingchen.zhang@20200703 原因小于两条分析genus会报错
            genus_path = os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls")
            if os.path.exists(genus_path):
                with open(genus_path) as f:
                    for line in f:
                        line_counts += 1
            if line_counts > 2:
                self.anno.on('end', self.run_beta)
                self.anno.on("end", self.run_pan_core)
                self.anno.on("end", self.run_rank_abundance)
                self.anno.on('end', self.run_composition)
                if self.count_samples >= 500:###应产品线张俊彪需求大于500样本不跑稀释曲线分析，原因是太耗时
                    self.on_rely([self.alpha, self.beta, self.composition, self.pan_core, self.rank_abundance], self.set_output)
                else:
                    self.anno.on("end", self.run_rarefaction)
                    self.on_rely([self.alpha, self.beta, self.composition, self.pan_core, self.rank_abundance, self.rarefaction], self.set_output)
            else:
                if self.count_samples >= 500:
                    self.on_rely([self.alpha], self.set_output)
                else:
                    self.anno.on("end", self.run_rarefaction)
                    self.on_rely([self.alpha, self.rarefaction], self.set_output)
            self.run_taxon_format()
        if self.update_info:
            self.logger.info("分析结果异常处理：{}".format(self.update_info))
            self.step._info += self.update_info
            self.step._has_state_change = True

        task_info = self.api.api('task_info.metaasv_task_info')
        task_info.add_sample_numbers(self._sheet.id, self.count_samples)

    def run_basic_analysis(self):
        """
        运行基础分析
        :return:
        """
        self.run_alpha()
        if not self.option('group').is_set:
            group_table = os.path.join(self.work_dir, "group_file.xls")
            get_group_from_table(os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"), group_table)
        if "unifrac" in self.option('dis_method'):
            self.run_beta()
            self.run_pan_core()
            self.run_rank_abundance()
            if self.count_samples >= 500:###应产品线张俊彪需求大于500样本不跑稀释曲线分析，原因是太耗时
                self.on_rely([self.alpha, self.beta, self.pan_core, self.rank_abundance], self.set_output)
            else:
                self.run_rarefaction()
                self.on_rely([self.alpha, self.beta, self.pan_core, self.rank_abundance, self.rarefaction], self.set_output)
        else:
            line_counts = 0
            genus_path = os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls")
            if os.path.exists(genus_path):
                with open(genus_path) as f:
                    for line in f:
                        line_counts += 1
            if line_counts > 2:
                self.run_composition()
                self.run_beta()
                self.run_pan_core()
                self.run_rank_abundance()
                if self.count_samples >= 500:
                    self.on_rely([self.alpha, self.beta, self.composition, self.pan_core, self.rank_abundance], self.set_output)
                else:
                    self.run_rarefaction()
                    self.on_rely([self.alpha, self.beta, self.composition, self.pan_core, self.rank_abundance, self.rarefaction], self.set_output)
            else:
                if self.count_samples >= 500:
                    self.on_rely([self.alpha], self.set_output)
                else:
                    self.run_rarefaction()
                    self.on_rely([self.alpha, self.rarefaction], self.set_output)
        for module in self.basic_analysis:
            module.run()
            gevent.sleep(0)

    def run_basic_analysis2(self):
        """
        运行基础分析
        :return:
        """
        self.run_alpha()
        if not self.option('group').is_set:
            group_table = os.path.join(self.work_dir, "group_file.xls")
            get_group_from_table(os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"), group_table)
        if "unifrac" in self.option('dis_method'):
            self.run_beta()
            self.run_rank_abundance()
            if self.count_samples >= 500:
                self.on_rely([self.alpha, self.beta, self.rank_abundance], self.set_output)
            else:
                self.run_rarefaction()
                self.on_rely([self.alpha, self.beta, self.rank_abundance, self.rarefaction], self.set_output)
        else:
            line_counts = 0
            genus_path = os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls")
            if os.path.exists(genus_path):
                with open(genus_path) as f:
                    for line in f:
                        line_counts += 1
            if line_counts > 2:
                self.run_composition()
                self.run_beta()
                self.run_rank_abundance()
                if self.count_samples >= 500:
                    self.on_rely([self.alpha, self.beta, self.composition, self.rank_abundance], self.set_output)
                else:
                    self.run_rarefaction()
                    self.on_rely([self.alpha, self.beta, self.composition, self.rank_abundance, self.rarefaction], self.set_output)
            else:
                if self.count_samples >= 500:
                    self.on_rely([self.alpha], self.set_output)
                else:
                    self.run_rarefaction()
                    self.on_rely([self.alpha, self.rarefaction], self.set_output)
        for module in self.basic_analysis:
            module.run()
            gevent.sleep(0)

    def end(self):
        """
        结束、导表、文件连接、预警和上传结果文件
        :return:
        """
        is_pollu =check_pollution_pip(self.output_dir+'/ASVTaxon_summary/tax_summary_r',self._sheet.id,self.work_dir,is_sanger=self._sheet.UPDATE_STATUS_API)  #做预警物种检查
        #is pollu : 0,1 表示无预警，有预警.对应api 的 1,2
        self.logger.info('warning_pollu: %s' % (is_pollu + 1))
        self.add_task_option('warning_pollu', is_pollu + 1)  # 20200708
        self.logger.info(self.task_option_data)
        self.run_api()
        self.send_files()
        super(MetaAsvWorkflow, self).end()

    def run_sequence(self,dir):
        """
        对上传的文件夹进行解压，上传序列为压缩格式
        :return:
        """
        opts = {
            'fastq_dir': dir,
        }
        self.sequence.set_options(opts)
        self.sequence.run()

    def run_file_sequence(self,file):
        """
        对上传的文件进行解压，上传序列为压缩格式
        目前只支持gz和tar.gz格式的进行解压
        :return:
        """
        file_name = os.path.basename(file)
        file_name = file_name.strip(".gz") if re.search(r'\.gz', file_name) else  file_name.strip(".tar.gz")
        result_path = os.path.join(self.work_dir, "unzip_results")
        if not os.path.exists(result_path):
            os.mkdir(result_path)
        opts = {
            'fastq': file,
            'sample_name': file_name,
            'direction': "s",
            'result_path': result_path,
            "pipeline": "metaasv"
        }
        self.sequence.set_options(opts)
        self.sequence.run()

    def run_seq_extract(self):
        """
        进行样本检测fastq和fastq_dir
        :return:
        """
        if self.option("fastq_type") in ["fastq_dir"]:
            if self.is_zip:
                in_fastq = os.path.join(self.sequence.output_dir, 'data')
            else:
                in_fastq = self.option("fastq_file")
            opts = {
                "in_fastq": in_fastq
            }
        elif self.option("fastq_type") in ["fastq"]:
            if self.is_zip:
                file_path = ""
                dir = os.path.join(self.work_dir, 'unzip_results')
                for file in os.listdir(dir):
                    file_path = os.path.join(dir, file)
                    break
                in_fastq = file_path
            else:
                in_fastq = self.option("fastq_file").prop['path']
            opts = {
                "in_fastq": in_fastq,
            }
        self.seq_extract.set_options(opts)
        self.sample_rename.on("start", self.set_step, {'start': self.step.qc_stat})
        self.sample_rename.on("end", self.set_step, {'end': self.step.qc_stat})
        self.seq_extract.run()

    def split_sample_fastq(self):
        """
        fastq序列按照样本进行拆分
        :return:
        """
        if self.option("fastq_type") in ["fastq_dir"]:
            if self.is_zip:
                in_fastq = os.path.join(self.sequence.output_dir, 'data')
            else:
                in_fastq = self.option("fastq_file")
            opts = {
                "in_fastq": in_fastq
            }
        elif self.option("fastq_type") in ["fastq"]:
            if self.is_zip:
                file_path = ""
                dir = os.path.join(self.work_dir, 'unzip_results')
                for file in os.listdir(dir):
                    file_path = os.path.join(dir, file)
                    break
                in_fastq = file_path
            else:
                in_fastq = self.option("fastq_file").prop['path']
            opts = {
                "in_fastq": in_fastq,
            }
        self.split_sample.set_options(opts)
        self.split_sample.run()

    def run_sample_rename(self):
        """
        样本重命名，前端已经修改完成，直接从MongoDB中导出，但是需要根据修改后的名称进行样本的合并操作
        是为质控后的样本统计文件
        :return:
        """
        info_txt = self.seq_extract.work_dir + "/info.txt"
        opts = {
            "info_txt": info_txt,
            "task_id": self._sheet.id,
            "denoise_method": self.option("denoise_method"),
            "query_id" : str(self.option("query_id"))
        }
        if self.option("fastq_file").format == "sequence.fastq":
            s3_origin = self.get_origin_name()
            opts['file_origin'] = json.dumps(s3_origin)
        else:
            opts['mapping_file'] = os.path.join(self.work_dir, "remote_input/fastq_file/mapping_file.txt")
        self.logger.info("task_id: {}".format(self._sheet.id))
        self.sample_rename.set_options(opts)
        self.sample_rename.run()

    def run_samplecheck(self):
        """
        根据上传的原始序列信息统计表进行样本二次检查和合并序列文件
        合并fq
        :return:
        """
        if len(os.listdir(self.split_sample.output_dir)) == 0:
            gevent.sleep(5)
            self.logger.info("链接结果文件，否则会报错!")
        if self.option("raw_sequence").is_set:
            opts = {
                "info_path": self.sample_rename.output_dir + "/info_path.xls",
                "info_temp": self.sample_rename.output_dir + "/info_temp.xls",
                "raw_sequence": self.option("raw_sequence"),
                "in_fastq": self.split_sample.output_dir
            }
        else:
            opts = {
                "info_path": self.sample_rename.output_dir + "/info_path.xls",
                "info_temp": self.sample_rename.output_dir + "/info_temp.xls",
                "in_fastq": self.split_sample.output_dir
            }
        self.sample_check.set_options(opts)
        self.sample_check.run()

    def run_denoise(self):
        """
        降噪
        :return:
        """
        fastq_dir = os.path.join(self.sample_check.output_dir, "fq")
        opts = {
                "fastq_dir": fastq_dir,
                "denoise_method": self.option('denoise_method'),
                "database": self.option("database")
        }
        if self.option("pipeline") == "functional_gene":
            opts["database"] = self.option("fungene_database")
        if self.option("ref_fasta").is_set:
            opts["ref_fasta"] = self.option("ref_fasta")
        if self.option("denoise_method") in ["Deblur"]:
            if self.option("truc_len") != 0:
                opts["truc_len"] = int(float(self.option("truc_len")))
            else:
                sample_info_path = self.sample_rename.output_dir + "/info_path.xls"
                min_value, mean_value, max_value = self.get_min_seq(sample_info_path)
                opts["truc_len"] = min_value
            if self.option("min_size"):
                opts["min_size"] = int(float(self.option("min_size")))
            if self.option("min_reads"):
                opts["min_reads"] = int(float(self.option("min_reads")))
        if self.option("max_ee"):
            opts["max_ee"] = int(float(self.option("max_ee")))
        if self.option("trunc_q"):
            opts["trunc_q"] = int(float(self.option("trunc_q")))

        self.denoise.set_options(opts)
        self.denoise.on("start", self.set_step, {'start': self.step.denoise})
        self.denoise.on("end", self.set_step, {'end': self.step.denoise})
        self.denoise.run()

    def run_framebot(self):
        opts = {
            "fasta": os.path.join(self.denoise.output_dir, "ASV_reps.fasta"),
            "database": self.option("fungene_database"),
            "acid_length": self.option("acid_length"),
            "seq_identity": self.option("seq_identity"),
            "ref_acid": self.option("ref_acid"),
            "path" : self.denoise.output_dir,
        }
        self.framebot.set_options(opts)
        self.framebot.on("start", self.set_step, {'start': self.step.framebot})
        self.framebot.on("end", self.set_step, {'end': self.step.framebot})
        self.framebot.run()

    def run_phylotree(self):
        """
        构建进化树 asv水平
        :return:
        """
        if self.pipline in ["pipeline2"]:
            self.phylo.set_options({
            "fasta_file": self.option("asv_upload")##降噪得到的代表序列
            })
        else:
            if self.option("pipeline") == "functional_gene":
                self.phylo.set_options({
                    "fasta_file": os.path.join(self.framebot.output_dir, "framebot_nucl.fasta")
                })
            else:
                self.phylo.set_options({
                    "fasta_file": os.path.join(self.denoise.output_dir, "ASV_reps.fasta")
                })
        self.phylo.run()

    def run_taxon_format(self):
        """
        对上传的数据库做format的处理
        :return:
        """
        self.format_taxon = self.add_tool("meta.filecheck.format_taxon")
        if self.option("database") in ["custom_mode"] or self.option("fungene_database") in ["custom_mode"]:
            opts = {"in_taxon_table": self.option('ref_taxon')}
            self.format_taxon.set_options(opts)
            self.format_taxon.on('end',self.run_taxon)
            self.format_taxon.run()
        else:
            self.run_taxon()

    def run_taxon(self):
        """
        物种注释分类和统计
        :return:
        """
        if self.option("pipeline") != "functional_gene":
            opts = {
                "anno_method": self.option("anno_method"),
                "database": self.option("database")
            }
            if self.pipline in ["pipeline2"]:
                opts["in_otu_table"] = self.option("asv_abundance")  ##统计丰度
                opts["fasta"] = self.option("asv_upload")  ##降噪得到的代表序列
            else:
                opts["fasta"] = os.path.join(self.denoise.output_dir, "ASV_reps.fasta")  ##降噪得到的代表序列
                opts["in_otu_table"] = os.path.join(self.denoise.output_dir, "ASV_table.xls")  ##统计丰度
            if self.option("anno_method") in ["multi_blast"]:  ##自写blast比对
                opts.update({
                    "identity": 100 * self.option("identity"),
                    "coverage": 100 * self.option("coverage"),
                })
            elif self.option("anno_method") in ["vsearch", "blast"]:  ##qiime2注释
                opts.update({
                    "identity": 100 * self.option("identity"),
                    "coverage": 100 * self.option("coverage"),
                })
                if self.pipline in ["pipeline1"]:
                    opts["qza_fasta"] = os.path.join(self.denoise.output_dir, "ASV_reps.qza")  ##降噪得到的代表序列
                    opts["asv_md5"] = os.path.join(self.denoise.output_dir, "ASV_md5.xls")  ## MD5值用于替换名称
            elif self.option("anno_method") in ["bayes"]:  ##qiime2注释
                opts.update({
                    "confidence": self.option("confidence")
                })
                if self.pipline in ["pipeline1"]:
                    opts["qza_fasta"] = os.path.join(self.denoise.output_dir, "ASV_reps.qza")  ##降噪得到的代表序列
                    opts["asv_md5"] = os.path.join(self.denoise.output_dir, "ASV_md5.xls")  ## MD5值用于替换名称
            else:  ##rdp注释老流程
                opts.update({
                    "confidence": self.option("confidence")
                })
            if self.option("database") in ["custom_mode"]:
                opts.update({
                    "ref_fasta": self.option("ref_fasta"),
                    "ref_taxon": self.format_taxon.option("out_taxon_table")
                })
        else:
            opts = {
                "anno_method": self.option("fungene_anno_method"),
                "database": self.option("fungene_database")
            }
            opts["fasta"] = self.framebot.option("framebot_nucl_fasta")  ##降噪得到的代表序列
            opts["in_otu_table"] = os.path.join(self.framebot.output_dir, "ASV_table.xls")  ##统计丰度
            if self.option("fungene_anno_method") in ["blastp"] and self.option("tax_database")  not in ["nr_v20210917"]:
                opts.update({
                    "fasta": self.framebot.option("framebot_prot_fasta"),
                    "blast": "blastp",
                    "coverage": 100 * self.option("fungene_coverage"),
                    "evalue": self.option("evalue"),
                    "query_type": "prot",
                    "meta_pipeline": "metaasv",
                })
            elif self.option("fungene_anno_method") in ["diamond2"] or self.option("tax_database") in ["nr_v20210917"]:
                if self.option("seq_style") == "nucl":
                    blast_type = "blastx"
                    query_type = "nucl"
                    opts["fasta"] = self.framebot.option("framebot_nucl_fasta")
                else:
                    blast_type = "blastp"
                    query_type = "prot"
                    opts["fasta"] = self.framebot.option("framebot_prot_fasta")
                opts.update({
                    "database": self.option("tax_database"),
                    "blast": blast_type,
                    "evalue": self.option("evalue"),
                    "query_type": query_type,
                    "meta_pipeline": "metaasv",
                    "anno_method": "diamond2"
                })
            elif self.option("fungene_anno_method") in ["multi_blast"]:  ##自写blast比对
                opts.update({
                    "identity": 100 * self.option("fungene_identity"),
                    "coverage": 100 * self.option("fungene_coverage"),
                })
            elif self.option("fungene_anno_method") in ["vsearch", "blast"]:  ##qiime2注释
                opts.update({
                    "identity": 100 * self.option("fungene_identity"),
                    "coverage": 100 * self.option("fungene_coverage"),
                })
                if self.pipline in ["pipeline1"]:
                    opts["qza_fasta"] = os.path.join(self.framebot.output_dir, "ASV_reps.qza")  ##降噪得到的代表序列
                    opts["asv_md5"] = os.path.join(self.framebot.output_dir, "ASV_md5.xls")  ## MD5值用于替换名称
            elif self.option("fungene_anno_method") in ["bayes"]:  ##qiime2注释
                opts.update({
                    "confidence": self.option("fungene_confidence")
                })
                if self.pipline in ["pipeline1"]:
                    opts["qza_fasta"] = os.path.join(self.framebot.output_dir, "ASV_reps.qza")  ##降噪得到的代表序列
                    opts["asv_md5"] = os.path.join(self.framebot.output_dir, "ASV_md5.xls")  ## MD5值用于替换名称
            else:  ##rdp注释老流程
                opts.update({
                    "confidence": self.option("confidence")
                })
            if self.option("fungene_database") in ["custom_mode"]:
                if self.option("seq_style") == "nucl":
                    opts.update({
                        "ref_fasta": self.option("ref_fasta"),
                        "ref_taxon": self.format_taxon.option("out_taxon_table")
                    })
                else:
                    opts.update({
                        "ref_fasta": self.option("ref_acid"),
                        "ref_taxon": self.format_taxon.option("out_taxon_table"),
                        "blast": "blastp",
                        "query_type":"prot",
                    })

        self.anno.set_options(opts)
        self.anno.on("start", self.set_step, {'start': self.step.taxassign})
        self.anno.on("end", self.set_step, {'end': self.step.taxassign})
        self.anno.run()

    def run_alpha(self):
        """
        Alpha多样性
        :return:
        """
        options = ({
            'otu_table': os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_asv.xls"),
            'indices': self.option('estimate_indices')
        })
        self.alpha = self.add_batch('meta.alpha_diversity.estimators', ignore_error=True, batch_type="tool")  # 出错也没关系
        self.alpha.set_options(options)
        self.basic_analysis.append(self.alpha)
        # self.alpha.run()

    def run_rarefaction(self):
        """
        rarefaction 稀释曲线分析
        :return:
        """
        options = ({
            'otu_table': os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_asv.xls"),
            'indices': self.option('rarefy_indices'),
            'freq': 100
        })
        self.rarefaction = self.add_batch('meta.alpha_diversity.rarefaction', ignore_error=True, batch_type="module")  # 出错也没关系
        self.rarefaction.set_options(options)
        self.basic_analysis.append(self.rarefaction)
        # self.rarefaction.run()

    def run_pan_core(self):
        """
        Pan_core物种分析
        :return:
        """
        opts = {
            'in_otu_table': os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"),
        }
        if self.option("group").is_set:
            opts['group_table'] = self.option("group")
        self.pan_core = self.add_batch('metaasv.pan_core', ignore_error=True, batch_type="tool")  # 出错也没关系
        self.pan_core.set_options(opts)
        self.basic_analysis.append(self.pan_core)
        # self.pan_core.run()

    def run_rank_abundance(self):
        """
        rank_abundance 排序
        :return:
        """
        opts = {
            "otu_table":  os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"),
            "method": "relative",
            "step" : 1
        }
        self.rank_abundance = self.add_batch('metaasv.rank_abundance', ignore_error=True, batch_type="tool")  # 出错也没关系
        self.rank_abundance.set_options(opts)
        self.basic_analysis.append(self.rank_abundance)
        # self.rank_abundance.run()

    def run_beta(self):
        """
        排序回归分析 pca、pcoa、nmds、hcluster
        :return:
        """
        opts = {
            'analysis': 'distance' + self.option('beta_analysis'),
            'otutable': os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"),
            "level": 'asv',
        }
        if 'pca' in self.beta_analysis:
            opts["permutations"] = 999
            opts["scale"] = "T"
        if self.option('dis_method') != "":
            opts["dis_method"] = self.option('dis_method')
        if self.count_otus:
            opts['phy_newick'] = self.phylo.option('phylo_tre').prop['path']
        if self.option('group').is_set:
            opts['group'] = self.option('group')
        self.beta = self.add_batch('metaasv.beta_diversity', ignore_error=True, batch_type="module")  # 出错也没关系
        self.beta.set_options(opts)
        self.basic_analysis.append(self.beta)
        # self.beta.run()

    def run_composition(self):
        """
        组成分析 barpie和heatmap， 默认跑Genus水平
        :return:
        """
        group_table = os.path.join(self.work_dir, "group_file.xls")
        opts = {
            'analysis': self.option('composition_type'),
            'abundtable': os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"),
        }
        if self.option('group').is_set:
            opts['group'] = self.option('group')
        else:
            opts['group'] = get_group_from_table(os.path.join(self.anno.output_dir, "ASVTaxon_summary", "tax_summary_a", "asv_taxon_Genus.full.xls"), group_table)
        if "heatmap" in self.option("composition_type").split(","):
            opts["method"] = "average"
            opts["sample_method"] = "average"
            opts["sample_distance"] = self.option("dis_method")
            opts["species_distance"] = self.option("dis_method")
            opts["species_number"] = "50"
        else:
            opts["others"] = 0.01  #bu zhaozhigang 20210413 配合页面
        self.composition = self.add_batch('metaasv.composition_analysis', ignore_error=True, batch_type="module")  # 出错也没关系
        self.composition.set_options(opts)
        self.basic_analysis.append(self.composition)
        # self.composition.run()

    def run_api(self):
        """
        导表
        :return:
        """
        self.logger.info("开始运行导表！")
        if self.pipline in ["pipeline1"]:
            self.export_specimen()
            self.export_group()
            self.export_datastat()
            self.export_sample_check()
            self.export_asv()
            self.export_phylo_tree()
            if os.path.exists(os.path.join(self.output_dir, "Alpha_diversity/Estimators")):
                self.export_alpha_diversity()
            if os.path.exists(os.path.join(self.output_dir, "Alpha_diversity/Rarefaction")):
                self.export_rarefaction()
            if os.path.exists(os.path.join(self.output_dir, "Pan_Core")):
                self.export_pan_core()
            if os.path.exists(os.path.join(self.output_dir, "Rank_abundance")):
                self.export_rank_abundance()
            if os.path.exists(os.path.join(self.output_dir, "Beta_diversity")):
                self.export_beta()
            if os.path.exists(os.path.join(self.output_dir, "CompositionAnalysis")):
                self.export_composition()
        elif self.pipline in ["pipeline2"]:
            self.export_specimen2()
            self.export_group()
            self.export_asv()
            self.export_phylo_tree()
            if os.path.exists(os.path.join(self.output_dir, "Alpha_diversity/Estimators")):
                self.export_alpha_diversity()
            if os.path.exists(os.path.join(self.output_dir, "Alpha_diversity/Rarefaction")):
                self.export_rarefaction()
            if os.path.exists(os.path.join(self.output_dir, "Pan_Core")):
                self.export_pan_core()
            if os.path.exists(os.path.join(self.output_dir, "Rank_abundance")):
                self.export_rank_abundance()
            if os.path.exists(os.path.join(self.output_dir, "Beta_diversity")):
                self.export_beta()
            if os.path.exists(os.path.join(self.output_dir, "CompositionAnalysis")):
                self.export_composition()

    def export_specimen(self):
        """
        运行样本导表
        :return:
        """
        self.logger.info("开始运行specimen导表")
        sample_info_path = self.output_dir + "/QC_Stat/valid_sequence.txt"
        api_specimen = self.api.api("metaasv.specimen")
        main_id = api_specimen.add_samples_info()
        self.spname_spid = api_specimen.add_specimen_detail(main_id, sample_info_path)

    def export_specimen2(self):
        """
        流程2运行样本导表
        :return:
        """
        self.logger.info("开始运行specimen导表")
        asv_abundance = self.option("asv_abundance").prop['path']
        asv_table = os.path.join(self.work_dir, "valid_sequence.txt")
        with open(asv_abundance, 'r') as f, open(asv_table, 'w') as w:
            line_list = f.readline().strip().split("\t")
            line_list = [x.strip() for x in line_list]
            w.write("\n".join(line_list))
        api_specimen = self.api.api("metaasv.specimen")
        main_id = api_specimen.add_samples_info()
        self.spname_spid = api_specimen.add_specimen_detail(main_id, asv_table)

    def export_group(self):
        """
        运行分组导表
        :return:
        """
        self.logger.info("开始运行specimen_group导表")
        api_group = self.api.api("metaasv.group")
        if self.option("group").is_set:
            group_id_list = api_group.add_ini_group_table(self.option('group').prop["path"], spname_spid=self.spname_spid, sort_samples=False)
            self.group_id = str(group_id_list[0])
        else:
            #如果为All的话，前端直接调用，不要再导入MongoDB的分组方案 @20200529
            # group_table = os.path.join(self.work_dir, "group_file.xls")
            # api_group.add_ini_group_table(group_table, self.spname_spid, sort_samples=False)
            self.group_id = "all"

    def export_datastat(self):
        """
        运行质控导表
        :return:
        """
        self.logger.info("开始运行datastat导表")
        api_datastat = self.api.api("metaasv.data_stat")
        main_id = api_datastat.add_datastat()
        if self.option('raw_sequence').is_set:
            raw_sequence_path = self.sample_check.output_dir + "/raw_sequence.txt"
            # raw_sequence_path = self.option("raw_sequence").prop["path"]
            column_number = self.option("raw_sequence").prop["column_number"]
            api_datastat.add_datastat_detail(main_id, raw_sequence_path, "raw", column_number=column_number)
        clean_path = self.output_dir + "/QC_Stat/valid_sequence.txt"
        api_datastat.add_datastat_clean(main_id, clean_path, "clean")
        if os.path.exists(self.output_dir + "/QC_Stat/info_path.xls"):
            os.remove(self.output_dir + "/QC_Stat/info_path.xls")
        denoise_path = self.output_dir + "/Denoise_Stat/%s_sequence_info.txt"%(self.option("denoise_method"))
        if os.path.exists(denoise_path):
            api_datastat.add_datastat_denoise(main_id, denoise_path, "denoise", self.option("denoise_method"))

    def export_sample_check(self):
        """
        运行样本检测结果
        如果已经做过样本检测，不再进行导表
        :return:
        """
        self.logger.info("开始运行sample_check导表")
        api_sample_check = self.api.api("metaasv.sample_check")
        result = api_sample_check.check_sample(self._sheet.id, self.option("query_id"))
        if result:
            pass
        else:
            params = {
                'fastq_type': self.option("fastq_type"),
                "fastq_file": self.option("fastq_file").prop['path'],
                "task_id": self._sheet.id,
                "submit_location": "sample_check",
                "task_type": "1",
                "query_id" : str(self.option("query_id"))
            }
            if self.option("fastq_file").format == "sequence.fastq":
                s3_origin = self.get_origin_name()
            else:
                s3_origin = self.get_dir_name()
            main_id = api_sample_check.add_seq_sample(params, self._sheet.id, name="Sample_check_Origin", query_id=self.option("query_id"))
            info_path = self.seq_extract.work_dir + "/info.txt"
            api_sample_check.add_seq_sample_detail(info_path, main_id, s3_name=s3_origin)

    def export_asv(self):
        """
        运行asv导表
        :return:
        """
        self.logger.info("开始运行asv导表")
        self.api_common = self.api.api("metaasv.common_api")
        otu_path = self.output_dir + "/ASVTaxon_summary/asv_taxon.xls"
        otu_absolute = self.output_dir + "/ASVTaxon_summary/tax_summary_a/asv_taxon_asv.full.xls"
        otu_phy = self.output_dir + "/ASVTaxon_summary/tax_summary_a/asv_taxon_Phylum.full.xls"
        otu_gen = self.output_dir + "/ASVTaxon_summary/tax_summary_a/asv_taxon_Genus.full.xls"
        if  self.pipline in ["pipeline2"]:
            rep_path = self.option("asv_upload").prop["path"]
        else:
            rep_path = self.output_dir + "/ASV/ASV_reps.fasta"
        if not os.path.isfile(otu_path):
            self.logger.error("找不到报告文件:{}".format(otu_path))
            self.set_error("找不到报告文件")
        self.level_id = 9
        params = {
            "group_id": "all",
            "size": "",
            "submit_location": 'asv',
            "filter_json": "",
            "task_type": "2",
        }
        self.asv_id = self.api_common.add_otu_table(otu_path, major=True, rep_path=rep_path, spname_spid=self.spname_spid, params=params)
        api_otu_level = self.api.api("metaasv.sub_sample")
        api_otu_level.add_sg_otu_detail_level(otu_absolute, self.asv_id, 9)
        api_otu_level.add_otu_detail(otu_phy, self.asv_id, 3)
        api_otu_level.add_otu_detail(otu_gen, self.asv_id, 7)
        api_otu_level.add_sg_otu_seq_summary(otu_path, self.asv_id)
        # self.api_common.add_meta_status(table_id=self.asv_id, type_name="asv")

    def export_phylo_tree(self):
        """
        运行进化树导表
        :return:
        """
        api_tree = self.api.api("metaasv.phylo_tree")
        tree_path = self.output_dir + "/ASV/ASV_phylo.tre"
        if not os.path.isfile(tree_path):
            self.logger.error("找不到报告文件:{}".format(tree_path))
            self.set_error("找不到报告文件")
        if os.path.exists(tree_path):
            # params = {
            #     "asv_id": str(self.asv_id),
            #     "submit_location": "phylo_tree",
            #     "task_type": "2",
            #     "level_id": 9,
            #     "method":"NJ",
            #     "group_id": self.group_id
            # }
            params = ""
            main_id = api_tree.add_phylo_tree_info_for_meta(asv_id=self.asv_id, params=params)
            api_tree.add_phylo_tree_info_workflow(main_id, tree_path)
            # self.api_common.add_meta_status(table_id=main_id, type_name="phylo_tree")

    def export_alpha_diversity(self):
        """
        运行alpha多样性导表
        :return:
        """
        self.logger.info("开始运行alpha_diversity导表")
        api_est = self.api.api("metaasv.estimator")
        est_path = self.output_dir + "/Alpha_diversity/Estimators/estimators.xls"
        if not os.path.isfile(est_path):
            self.logger.error("找不到报告文件:{}".format(est_path))
            self.set_error("找不到报告文件")
        indice = sorted(self.option("estimate_indices").split(','))
        self.level_id = 9
        params = {
            "level_id": self.level_id,
            "index_type": ','.join(indice),
            'submit_location': 'alpha_diversity',
            'task_type': "2",
            'group_id': self.group_id
        }
        est_id = api_est.add_est_table(est_path, major=True, level=self.level_id, otu_id=str(self.asv_id),params=params, spname_spid=self.spname_spid, indices=self.option("estimate_indices"), group_id=self.group_id)
        api_est.add_est_bar(est_path, est_id,indices=self.option("estimate_indices"))
        # self.api_common.add_meta_status(table_id=str(est_id), type_name='alpha_diversity')

    def export_rarefaction(self):
        """
        运行稀释曲线分析导表
        :return:
        """
        self.logger.info("开始运行rarefaction导表")
        api_rare = self.api.api("metaasv.rarefaction")
        rare_path = self.output_dir + "/Alpha_diversity/Rarefaction/"  # 此路径比较准确
        params = {
            "level_id": self.level_id,
            "index_type": self.option("rarefy_indices"),
            'submit_location': 'rarefaction',
            'task_type': "2",
            'group_id': self.group_id
        }
        if self.option("group").is_set:
            group = self.option("group").prop["path"]
        else:
            group = os.path.join(self.work_dir, "group_file.xls")
        rare_id = api_rare.add_rare_table(rare_path, level=self.level_id, otu_id=str(self.asv_id),params=params, spname_spid=self.spname_spid,group_id=self.group_id)
        api_rare.add_rarefaction_detail(rare_id, rare_path, self.option("rarefy_indices"), group=group)
        # self.api_common.add_meta_status(table_id=str(rare_id), type_name='rarefaction')

    def export_pan_core(self):
        """
        运行pancore物种分析导表
        :return:
        """
        self.logger.info("开始运行pan_core导表")
        api_pan_core = self.api.api("metaasv.pan_core")
        params = {
            "level_id": 7,
            'submit_location': 'pan_core',
            "asv_id": str(self.asv_id),
            'task_type': "2",
            'group_id': self.group_id
        }
        main_id = api_pan_core.add_pan_core(params=params, level_id=7, asv_id=str(self.asv_id),spname_spid=self.spname_spid,group_id=self.group_id)
        pan_path = self.output_dir + "/Pan_Core/Pan.richness.xls"
        core_path = self.output_dir + "/Pan_Core/Core.richness.xls"
        api_pan_core.add_pan_core_detail(pan_path, main_id, "pan")
        api_pan_core.add_pan_core_detail(core_path, main_id, "core")
        # self.api_common.add_meta_status(table_id=str(main_id), type_name='pan_core', submit_location="pancore")

    def export_rank_abundance(self):
        """
        运行rank_abundance物种分析导表
        :return:
        """
        self.logger.info("开始运行rank_abundance导表")
        api_rank = self.api.api("metaasv.rank_abundance")
        params = {
                'level_id': 7,
                'submit_location': "rank_abundance",
                "asv_id": str(self.asv_id),
                'task_type': "2",
                "group_id": self.group_id}
        rank_file = os.path.join(self.output_dir, "Rank_abundance/Rank_abundance.xls")
        main_id = api_rank.add_rank(asv_id=self.asv_id, params=params, name="Rank_abundance_Origin", spname_spid=self.spname_spid,group_id=self.group_id)
        api_rank.add_rank_detail(rank_file, main_id)
        # self.api_common.add_meta_status(table_id=str(main_id), type_name='rank_abundance',submit_location="rankabundace")

    def export_beta(self):
        """
        运行beta多样性导表
        :return:
        """
        self.logger.info("开始运行beta导表")
        beta_diversity = self.api.api("metaasv.beta_diversity")
        if 'hcluster' in self.option('beta_analysis').split(','):
            params = {
                'level_id': 7,
                'submit_location': "hcluster",
                "asv_id": str(self.asv_id),
                'task_type': "2",
                'group_id': self.group_id,
                "hcluster_method":"average",
                "distance_algorithm": self.option("dis_method")
                    }
            hcluster_path = self.output_dir + "/Beta_diversity/Hcluster/hcluster.tre"
            if not os.path.isfile(hcluster_path):
                self.logger.error("找不到报告文件:{}".format(hcluster_path))
                self.set_error("找不到报告文件",)
            dir_path = self.output_dir + "/Beta_diversity/Hcluster"
            link_file(os.path.join(self.output_dir, "Beta_diversity/Distance/%s_asv_taxon_Genus.full.xls"%self.option("dis_method")), self.output_dir + "/Beta_diversity/Hcluster/%s_asv_taxon_asv.xls"%self.option("dis_method"))
            if self.option("group").is_set:
                group_file = self.option("group").prop["path"]
            else:
                group_file = os.path.join(self.work_dir, "group_file.xls")
            main_id = beta_diversity.add_beta_multi_analysis_result(dir_path,"hcluster", main=True, otu_id=str(self.asv_id), level=7, params=params,spname_spid=self.spname_spid, group_id=self.group_id, group_file=group_file)
            # self.api_common.add_meta_status(table_id=main_id, type_name='hcluster')  # 主表写入没有加name，所以此处table_name固定

        for ana in self.option('beta_analysis').split(','):
            if ana in ['pca', 'pcoa', 'nmds']:
                api_betam = self.api.api("metaasv.beta_diversity")
                params = {
                    'level_id': 7,
                    'submit_location': ana,
                    "asv_id": str(self.asv_id),
                    'group_id': self.group_id,
                    'task_type': "2",
                    'diff_test_method': "none",
                    "change_times" : str(999)
                }
                if ana in ['nmds','pcoa']:
                    params['distance_algorithm'] = self.option('dis_method')
                else:
                    params["scale"] = "T"
                self.beta_dict = {"pca": "Pca", "pcoa": "Pcoa", "nmds": "Nmds"}
                dir_path = self.output_dir + "/Beta_diversity/" + self.beta_dict[ana]
                main_id = api_betam.add_beta_multi_analysis_result(dir_path=dir_path, analysis=ana,
                                                                   main=True,otu_id=self.asv_id, params=params,
                                                                   spname_spid=self.spname_spid,group_id=self.group_id)
                # self.api_common.add_meta_status(table_id=main_id, type_name=ana)
                self.logger.info('set output beta %s over.' % ana)

    def export_composition(self):
        """
        运行组成分析导表
        :return:
        """
        self.logger.info("开始运行composition导表")
        for ana_type in self.option('composition_type').split(','):
            api_barpie = self.api.api("metaasv.barpie")
            params = {
                'level_id': 7,
                'submit_location': ana_type,
                "asv_id": str(self.asv_id),
                'task_type': "2",
                'group_id': self.group_id,
            }
            if ana_type in ['barpie']:
                params["combine_value"] = str(0.01)
                params["group_method"] = "none"
                composition_psth = os.path.join(self.output_dir, "CompositionAnalysis/CommunityBarPie/taxa.percents.table.xls")
                if os.path.exists(composition_psth):
                    main_id = api_barpie.add_barpie(params,from_otu_table=self.asv_id,spname_spid=self.spname_spid,group_id=self.group_id)
                    api_barpie.add_sg_otu_detail(composition_psth, main_id)
                # self.api_common.add_meta_status(table_id=main_id, type_name="barpie")
            elif ana_type in ['heatmap']:
                params["species_method"] = "average"
                params["sample_method"] = "average"
                params["top"] = "50"
                params["distance_method"] = self.option('dis_method')
                params["sample_distance_method"] = self.option('dis_method')
                api_heatmap = self.api.api("metaasv.composition_heatmap")
                species_tree_path = os.path.join(self.output_dir, "CompositionAnalysis/CommunityHeatmap/species_hcluster.tre")
                if os.path.exists(species_tree_path):
                    with open(species_tree_path, "r") as f:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [i[1].split("; ")[-1].strip() for i in raw_samp]
                sample_tree_path = os.path.join(self.output_dir, "CompositionAnalysis/CommunityHeatmap/specimen_hcluster.tre")
                if os.path.exists(sample_tree_path):
                    with open(sample_tree_path, "r") as f:
                        sample_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', sample_tree)
                        sample_list = [i[1].split("; ")[-1].strip() for i in raw_samp]
                if self.option("group").is_set:
                    sample_name = open(self.option("group").prop["path"], "r")
                    content = sample_name.readlines()
                    for f in content:
                        f = f.strip("\n")
                        arr = f.strip().split("\t")
                        if arr[0] != "#sample":
                            if arr[0] not in sample_list:
                                sample_list.append(arr[0].split("; ")[-1].strip())
                insert_asv_table = os.path.join(self.output_dir, "CompositionAnalysis/CommunityHeatmap/taxa.table.xls")
                insert_asv_percents_table = os.path.join(self.output_dir, "CompositionAnalysis/CommunityHeatmap/taxa.percents.table.xls")
                main_id = api_heatmap.add_heatmap(params,from_otu_table=self.asv_id,spname_spid=self.spname_spid, group_id=self.group_id)
                if os.path.exists(insert_asv_table):
                    api_heatmap.add_heatmap_detail(insert_asv_table, main_id, "absolute",specimen_sorts=sample_list,species_sorts=species_list)
                if os.path.exists(insert_asv_percents_table):
                    api_heatmap.add_heatmap_detail(insert_asv_percents_table, main_id, "relative",specimen_sorts=sample_list,species_sorts=species_list)
                if os.path.exists(species_tree_path):
                    api_heatmap.insert_tree_table(species_tree_path, main_id, "species")
                if os.path.exists(sample_tree_path):
                    api_heatmap.insert_tree_table(sample_tree_path, main_id, "specimen")
                # self.api_common.add_meta_status(table_id=main_id, type_name="heatmap")

    def get_min_seq(self, table_path):
        """
        根据质控后的文件统计最小值
        :param table_path:
        :return:
        """
        origin_min_list = []
        origin_mean_list = []
        origin_max_list = []
        with open(table_path, 'r') as f:
            for line in f:
                line = line.strip().split("\t")
                if line[0] == "Sample_Name":
                    pass
                else:
                    if line[4] not in origin_min_list:
                        origin_min_list.append(int(float(line[4])))
                    if line[3] not in origin_mean_list:
                        origin_mean_list.append(int(float(line[3])))
                    if line[5] not in origin_max_list:
                        origin_max_list.append(int(float(line[5])))
            min_value = min(origin_min_list)
            if min_value < 90:
                min_value = 90
            else:
                min_value = min_value
            mean_value = min(origin_mean_list)
            if mean_value < 90:
                mean_value = 90
            else:
                mean_value = mean_value
            max_value = min(origin_max_list)
            if max_value < 90:
                max_value = 90
            else:
                max_value = max_value
        return(min_value, mean_value, max_value)

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        if self.pipline in ["pipeline1"]:
            link_dir(self.sample_rename.output_dir, self.output_dir + "/QC_Stat")
            if self.option("raw_sequence").is_set:
                link_file(self.sample_check.output_dir + "/raw_sequence.txt", os.path.join(self.output_dir, "QC_Stat", "raw_sequence.txt"))
            os.rename(self.output_dir + "/QC_Stat/info_path.xls", os.path.join(self.output_dir, "QC_Stat", "valid_sequence.txt"))
            if os.path.exists(os.path.join(self.output_dir, "Denoise_Stat")):
                shutil.rmtree(os.path.join(self.output_dir, "Denoise_Stat"))
            os.mkdir(os.path.join(self.output_dir, "Denoise_Stat"))
            if os.path.exists(os.path.join(self.output_dir, "ASV")):
                shutil.rmtree(os.path.join(self.output_dir, "ASV"))
            os.mkdir(os.path.join(self.output_dir, "ASV"))
            for file in os.listdir(self.denoise.output_dir):###链接降噪结果到对应的文件夹
                if re.search(r"{}".format(self.option("denoise_method")), file):
                    link_file(os.path.join(self.denoise.output_dir, file), os.path.join(self.output_dir, "Denoise_Stat", file))
                else:
                    link_file(os.path.join(self.denoise.output_dir, file), os.path.join(self.output_dir, "ASV", file))
            if os.path.exists(os.path.join(self.phylo.output_dir, "phylo.tre")):
                link_file(os.path.join(self.phylo.output_dir, "phylo.tre"), os.path.join(self.output_dir, "ASV", "ASV_phylo.tre"))
        else:
            if os.path.exists(os.path.join(self.phylo.output_dir, "phylo.tre")):
                if not os.path.exists(self.output_dir + "/ASV"):
                    os.mkdir(self.output_dir + "/ASV")
                link_file(os.path.join(self.phylo.output_dir, "phylo.tre"), os.path.join(self.output_dir, "ASV", "ASV_phylo.tre"))
        link_dir(self.anno.output_dir + "/ASVTaxon_summary", self.output_dir + "/ASVTaxon_summary")
        link_dir(self.anno.output_dir + "/Tax_assign", self.output_dir + "/Tax_assign")
        if os.path.exists(os.path.join(self.alpha.work_dir, "Estimators/output")):
            if len(os.listdir(os.path.join(self.alpha.work_dir, "Estimators/output"))) != 0:
                link_dir(os.path.join(self.alpha.work_dir, "Estimators/output"), self.output_dir + "/Alpha_diversity/Estimators")
        if os.path.exists(os.path.join(self.rarefaction.work_dir, "Rarefaction/output")):
            if len(os.listdir(os.path.join(self.rarefaction.work_dir, "Rarefaction/output"))) != 0:
                link_dir(os.path.join(self.rarefaction.work_dir, "Rarefaction/output"), self.output_dir + "/Alpha_diversity/Rarefaction")
        if os.path.exists(os.path.join(self.pan_core.work_dir, "PanCore/output")):
            if len(os.listdir(os.path.join(self.pan_core.work_dir, "PanCore/output"))) != 0:
                link_dir(os.path.join(self.pan_core.work_dir, "PanCore/output"), self.output_dir + "/Pan_Core")
        if os.path.exists(os.path.join(self.rank_abundance.work_dir, "RankAbundance/output")):
            if len(os.listdir(os.path.join(self.rank_abundance.work_dir, "RankAbundance/output"))) != 0:
                link_dir(os.path.join(self.rank_abundance.work_dir, "RankAbundance/output"), self.output_dir + "/Rank_abundance")
        if os.path.exists(os.path.join(self.beta.work_dir, "BetaDiversity/output")):
            if len(os.listdir(os.path.join(self.beta.work_dir, "BetaDiversity/output"))) != 0:
                link_dir(os.path.join(self.beta.work_dir, "BetaDiversity/output"), self.output_dir + "/Beta_diversity")
        if os.path.exists(os.path.join(self.composition.work_dir, "CompositionAnalysis/output")):
            if len(os.listdir(os.path.join(self.composition.work_dir, "CompositionAnalysis/output"))) != 0:
                link_dir(os.path.join(self.composition.work_dir, "CompositionAnalysis/output"), self.output_dir + "/CompositionAnalysis")
        self.logger.info("完成结果文件目录设置!")
        self.end()

    def send_files(self):
        """
        上传结果文件目录到s3存储
        :return:
        """
        repaths = [
            [".", "", "基础分析结果文件夹", 0, ""],
            ["QC_stat", "", "样本数据统计文件目录", 0, ""],
            ["QC_Stat/valid_sequence_info.txt", "txt", "各样本优化序列信息统计表", 0, ""],
            ["Denoise_Stat", "", "单个样本碱基质量统计目录", 0, ""],
            ["Denoise_Stat/DADA2_stats.qza", "", "DADA2降噪后各样本序列信息统计", 0, ""],
            ["Denoise_Stat/DADA2_stats.qzv", "txt", "DADA2降噪后各样本序列信息统计", 0, ""],
            ["Denoise_Stat/DADA2_sequence_info.txt", "txt", "降噪后各样本序列信息统计表", 0, ""],
            ["Denoise_Stat/Deblur_stats.qza", "", "Deblur降噪后各样本序列信息统计", 0, ""],
            ["Denoise_Stat/Deblur_stats.qzv", "txt", "Deblur降噪后各样本序列信息统计", 0, ""],
            ["Denoise_Stat/Deblur_sequence_info.txt", "txt", "降噪后各样本序列信息统计表", 0, ""],
            ["ASV", "", "ASV聚类结果文件目录", 0, ""],
            ["ASV/ASV_reps.fasta", "sequence.fasta", "ASV代表序列文件", 0, ""],
            ["ASV/ASV_reps.qza", "metaasv.qza", "ASV代表序列qza文件", 0, ""],
            ["ASV/ASV_reps.qzv", "metaasv.qzv", "ASV代表序列qzv文件", 0, ""],
            ["ASV/ASV_table.xls", "xls", "ASV代表序列丰度表", 0, ""],
            ["ASV/ASV_md5.xls", "xls", "ASV代表序列的md5值", 0, ""],
            ["ASV/ASV_table.biom", 'meta.otu.biom', "biom格式的ASV代表序列丰度表", 0, ""],
            ["ASV/ASV_table.qza", "meta.otu.otu_table", "ASV代表序列丰度表", 0, ""],
            ["ASV/ASV_table.qzv", "meta.otu.otu_table", "ASV代表序列丰度表", 0, ""],
            ["ASV/ASV_phylo.tre", "graph.newick_tree", "ASV代表序列进化树文件", 0, ""],
            ["Tax_assign", "", "ASV物种注释结果目录", 0, ""],
            ["Tax_assign/seqs_tax_assignments.txt", "taxon.seq_taxon", "ASV物种注释结果文件", 0, ""],
            ["ASVTaxon_summary", "", "ASV物种分类统计结果目录", 0, ""],
            ["ASVTaxon_summary/asv_taxon.biom", "meta.otu.biom", "Biom格式的单级物种分类学统计结果", 0, ""],
            ["ASVTaxon_summary/asv_taxon.xls", "meta.otu.otu_table", "ASV物种分类统计表", 0, ""],
            ["ASVTaxon_summary/asv_summary.xls", "meta.otu.otu_table", "各样本中ASV数目统计", 0, ""],
            ["ASVTaxon_summary/asv_summary_a", "meta.otu.tax_summary_dir", "样本中各分类学水平物种的绝对丰度统计表", 0, ""],
            ["ASVTaxon_summary/asv_summary_r", "meta.otu.tax_summary_dir", "样本中各分类学水平物种的相对丰度统计表", 0, ""],
            ["Alpha_diversity", "", "Alpha多样性分析结果文件", 0, ""],
            ["Alpha_diversity/Estimators", "", "Alpha多样性指数分析结果目录", 0, ""],
            ["Alpha_diversity/Estimators/estimators.xls", "xls", "单个样本多样性指数表", 0, ""],
            ["Alpha_diversity/Rarefaction", "", "稀释曲线分析结果目录", 0, ""],
            ["Beta_diversity", "", "Beta多样性分析结果文件", 0, ""],
            ["Beta_diversity/Distance", "", "距离矩阵计算结果目录", 0, ""],
            ["Beta_diversity/Hcluster", "", "样本层级聚类分析结果目录", 0, ""],
            ["Beta_diversity/Hcluster/hcluster.tre", "graph.newick_tree", "样本层级聚类树文件", 0, ""],
            ["Beta_diversity/Nmds", "", "NMDS分析结果目录", 0, ""],
            ["Beta_diversity/Nmds/nmds_sites.xls", "xls", "样本各维度坐标", 0, ""],
            ["Beta_diversity/Nmds/nmds_stress.xls", "xls", "样本特征拟合度值", 0, ""],
            ["Beta_diversity/Pca", "", "PCA分析结果目录", 0, ""],
            ["Beta_diversity/Pca/pca_importance.xls", "xls", "主成分解释度表", 0, ""],
            ["Beta_diversity/Pca/pca_rotation.xls", "xls", "PCA主成分贡献度表", 0, ""],
            ["Beta_diversity/Pca/pca_rotation_all.xls", "xls", "全部物种主成分贡献度表", 0, ""],
            ["Beta_diversity/Pca/pca_sites.xls", "xls", "样本各成分轴坐标", 0, ""],
            ["Beta_diversity/Pcoa", "", "PCoA分析结果目录", 0, ""],
            ["Beta_diversity/Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, ""],
            ["Beta_diversity/Pcoa/pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比", 0, ""],
            ["Beta_diversity/Pcoa/pcoa_sites.xls", "xls", "样本坐标表", 0, ""],
            ["Pan_Core", "", "Pan/core物种分析结果目录", 0, ""],
            ["Pan_Core/Pan.richness.xls", "xls", "Pan物种分析结果表", 0, ""],
            ["Pan_Core/Core.richness.xls", "xls", "Core物种分析结果表", 0, ""],
            ["Rank_abundance", "", "Rank_abundance曲线分析结果目录", 0, ""],
            ["Rank_abundance/Rank_abundance.xls", "xls", "Rank_abundance曲线表", 0, ""],
            ["CompositionAnalysis", "", "群落组成分析结果文件", 0, ""],
            ["CompositionAnalysis/CommunityBarPie", "", "群落Bar图和Pie图结果目录", 0, ""],
            ["CompositionAnalysis/CommunityBarPie/taxa.table.xls", "xls", "各分组/样本物种丰度结果表", 0, ""],
            ["CompositionAnalysis/CommunityBarPie/taxa.percents.table.xls", "xls", "各分组/样本物种分布比例", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap", "", "群落Heatmap图结果目录", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap/Heatmap.taxa.table.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap/Sample_hcluster.tre", "", "样本聚类树", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap/Species_hcluster.tre", "xls", "物种聚类树", 0, ""],
        ]
        regexps = [
            [r"Alpha_diversity/Rarefaction/.+/otu.*.r_#.xls" , "", "每个样本的不同指数稀释性曲线表", 0, ""],
            [r'Beta_diversity/Distance/%s.*\.xls$' % self.option('dis_method'), 'meta.beta_diversity.distance_matrix', '样本距离矩阵文件', 0, ""],
            ["ASVTaxon_summary/tax_summary_a/.+\.biom$", "meta.otu.biom", "Biom格式的单级物种分类学统计结果(absolute)", 0, ""],
            ["ASVTaxon_summary/tax_summary_a/.+\.xls$", "xls", "单级物种分类统计表(absolute)", 0, ""],
            ["ASVTaxon_summary/tax_summary_a/.+\.full\.xls$", "xls", "多级物种分类统计表(absolute)", 0, ""],
            ["ASVTaxon_summary/tax_summary_r/.+\.biom$", "meta.otu.biom", "Biom格式的单级物种分类学统计结果(relative)", 0, ""],
            ["ASVTaxon_summary/tax_summary_r/.+\.xls$", "xls", "单级物种分类学统计表(relative)", 0, ""],
            ["ASVTaxon_summary/tax_summary_r/.+\.full\.xls$", "xls", "多级物种分类学统计表(relative)", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap/Sample_%s.xls" % self.option('dis_method'), "", "样本距离矩阵文件", 0, ""],
            ["CompositionAnalysis/CommunityHeatmap/Species_%s.xls" % self.option('dis_method'), "xls", "物种距离矩阵文件", 0, ""],
        ]
        for i in self.option("rarefy_indices").split(","):
            dir_code_list = {
                "sobs": "",
                "ace": "",
                "chao": "",
                "shannon": "",
                "simpson": "",
                "shannoneven": "",
                "simpsoneven": "",
                "coverage":""}
            repaths.append(["./Alpha_diversity/Rarefaction/{}".format(i), "文件夹", "稀释曲线分析结果目录", 0, dir_code_list[i]])
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

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
        in_fastq = str(json_dict['options']['fastq_file']).split("||")[-1]
        true_file_name = in_fastq.split("||")[-1].split("{")[0]
        new_file_name = in_fastq.split("||")[-1].split("{")[1].rstrip("}")
        if new_file_name not in s3_origin:
            s3_origin[new_file_name] = true_file_name
        return s3_origin

    def get_dir_name(self):
        """
        从mapping_file文件中获取s3_name与改名后的文件名称的对应关系
        :return:
        """
        s3_origin = {}
        mapping_file = os.path.join(self.work_dir, "remote_input/fastq_file/mapping_file.txt")
        with open(mapping_file, 'r') as mm:
            js = mm.read()
            json_dict = json.loads(js)
        fastq_list = json_dict['fastq_file']
        for fastq_dict in fastq_list:
            new_file_name = fastq_dict['alias']
            true_file_name = fastq_dict['file_path']
            if new_file_name not in s3_origin:
                s3_origin[new_file_name] = true_file_name
        return s3_origin