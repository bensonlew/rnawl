# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

"""多样性基础分析"""

from biocluster.workflow import Workflow
import os,re
import json
import shutil
import gevent
import shutil
from bson.objectid import ObjectId
from mbio.packages.metaasv.search_polution_by_list import check_pollution_pip
from mbio.packages.metaasv.common_function import link_dir,check_file,link_file,get_group_from_table


class MetaAsvApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        此功能负责 只是导表，不再重运行工作流
        """
        self._sheet = wsheet_object
        super(MetaAsvApiWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fastq_type', 'type': 'string'},  # 文件类型
            {'name': 'fastq_file', 'type': 'infile', 'format': 'sequence.fastq,sequence.fastq_dir'},# 输入的fastq文件或fastq文件夹
            {"name": "raw_sequence", "type": "infile", "format": "sequence.raw_sequence_txt"},  ## 原始序列信息表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  ##group分组表
            {'name': 'asv_upload', 'type': 'infile', 'format': 'sequence.fasta'},  ##流程2参数
            {'name': 'asv_abundance', 'type': 'infile', 'format': 'meta.otu.otu_table'},  ##流程1参数
            {'name': 'denoise_method', 'type': 'string'},  # 降噪方法选择
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'database_type', 'type': 'string'},  # 数据库类型
            {'name': 'anno_method', 'type': 'string'},  # 注释方法选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'identity', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'coverage', 'type': 'float', 'default': 0.8},  # 相似性值，范围0-1.
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            {"name": "estimate_indices", "type": "string",
             "default": "ace,chao,shannon,simpson,coverage,sobs,shannoneven,simpsoneven"},  ## 多样性指数类型
            {"name": "rarefy_indices", "type": "string",
             "default": "ace,chao,shannon,simpson,coverage,sobs,shannoneven,simpsoneven"},  # 稀释曲线指数类型
            {"name": "beta_analysis", "type": "string", "default": "pca,hcluster,pcoa,nmds"},  ##排序回归分析
            {"name": "composition_type", "type": "string", "default": "barpie,heatmap"},  ##组成分析
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},  ## 距离算法
            {"name": "file_list", "type": "string", "default": ""},  ##前端传过来用于改名和合并文件
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.bin_name = ''
        self.remote_dir = self._sheet.output + '/'
        if self.option("asv_upload").is_set and self.option("asv_abundance").is_set:
            self.pipline = "pipeline2"
        else:
            self.pipline = "pipeline1"

    def check_options(self):
        """
        检查参数
        """
        pass

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        pass

    def end(self):
        """
        结束
        :return:
        """
        self.run_api()
        self.send_files()
        super(MetaAsvApiWorkflow, self).end()

    def run_api(self):
        """
        导表
        :return:
        """
        self.logger.info("开始运行导表！")
        if self.pipline in ["pipeline1"]:
            link_dir("/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metaasv/upload",self.output_dir)
            self.export_specimen()
            self.export_group()
            self.export_datastat()
            # self.export_sample_check()
            self.export_asv()
            self.export_phylo_tree()
            self.export_alpha_diversity()
            self.export_rarefaction()
            if os.path.exists(os.path.join(self.output_dir, "Pan_Core")):
                self.export_pan_core()
            if os.path.exists(os.path.join(self.output_dir, "Rank_abundance")):
                self.export_rank_abundance()
            self.export_beta()
            if os.path.exists(os.path.join(self.output_dir, "CompositionAnalysis")):
                self.export_composition()
        elif self.pipline in ["pipeline2"]:
            link_dir("/mnt/ilustre/users/sanger-dev/workspace/20200602/MetaAsv_pipeline2_multi_blast123/output",self.output_dir)
            self.export_specimen2()
            self.export_group()
            self.export_asv()
            self.export_phylo_tree()
            self.export_alpha_diversity()
            self.export_rarefaction()
            if os.path.exists(os.path.join(self.output_dir, "Pan_Core")):
                self.export_pan_core()
            if os.path.exists(os.path.join(self.output_dir, "Rank_abundance")):
                self.export_rank_abundance()
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
        result = api_sample_check.check_sample(self._sheet.id)
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
            main_id = api_sample_check.add_seq_sample(params, self._sheet.id, name="Sample_check_Origin", query_id=self.option("query_id"))
            info_path = self.seq_extract.work_dir + "/info.txt"
            api_sample_check.add_seq_sample_detail(info_path, main_id)

    def export_asv(self):
        """
        运行asv导表
        :return:
        """
        self.logger.info("开始运行asv导表")
        self.api_common = self.api.api("metaasv.common_api")
        otu_path = self.output_dir + "/ASVTaxon_summary/asv_taxon.xls"
        otu_absolute = self.output_dir + "/ASVTaxon_summary/tax_summary_a/asv_taxon_asv.full.xls"
        if  self.pipline in ["pipeline2"]:
            rep_path = self.option("asv_upload").prop["path"]
        else:
            rep_path = self.output_dir + "/ASV/ASV_reps.fasta"
        if not os.path.isfile(otu_path):
            self.logger.error("找不到报告文件:{}".format(otu_path))
            self.set_error("找不到报告文件")
        params = {
            "group_id": self.group_id,
            "size": "",
            "submit_location": 'asv',
            "filter_json": "",
            "task_type": "2",
        }
        self.asv_id = self.api_common.add_otu_table(otu_path, major=True, rep_path=rep_path, spname_spid=self.spname_spid, params=params)
        api_otu_level = self.api.api("metaasv.sub_sample")
        api_otu_level.add_sg_otu_detail_level(otu_absolute, self.asv_id, 9)
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
            "level_id": self.level_id,
            'submit_location': 'pan_core',
            "asv_id": str(self.asv_id),
            'task_type': "2",
            'group_id': self.group_id
        }
        main_id = api_pan_core.add_pan_core(params=params, level_id=self.level_id, asv_id=str(self.asv_id),spname_spid=self.spname_spid,group_id=self.group_id)
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
                'level_id': self.level_id,
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
                'level_id': self.level_id,
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
            link_file(os.path.join(self.output_dir, "Beta_diversity/Distance/%s_asv_taxon_asv.xls"%self.option("dis_method")), self.output_dir + "/Beta_diversity/Hcluster/%s_asv_taxon_asv.xls"%self.option("dis_method"))
            if self.option("group").is_set:
                group_file = self.option("group").prop["path"]
            else:
                group_file = os.path.join(self.work_dir, "group_file.xls")
            main_id = beta_diversity.add_beta_multi_analysis_result(dir_path,"hcluster", main=True, otu_id=str(self.asv_id), level=self.level_id, params=params,spname_spid=self.spname_spid, group_id=self.group_id, group_file=group_file)
            # self.api_common.add_meta_status(table_id=main_id, type_name='hcluster')  # 主表写入没有加name，所以此处table_name固定

        for ana in self.option('beta_analysis').split(','):
            if ana in ['pca', 'pcoa', 'nmds']:
                api_betam = self.api.api("metaasv.beta_diversity")
                params = {
                    'level_id': self.level_id,
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
                'level_id': self.level_id,
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
                insert_asv_percents_table = os.path.join(self.output_dir, "CompositionAnalysis/CommunityHeatmap/taxa.mongo.xls")
                main_id = api_heatmap.add_heatmap(params,from_otu_table=self.asv_id,spname_spid=self.spname_spid, group_id=self.group_id)
                if os.path.exists(insert_asv_table):
                    api_heatmap.add_heatmap_detail(insert_asv_table, main_id, "absolute",specimen_sorts=sample_list,species_sorts=species_list)
                if os.path.exists(insert_asv_percents_table):
                    api_heatmap.add_heatmap_detail(insert_asv_percents_table, main_id, "relative",specimen_sorts=sample_list,species_sorts=species_list)
                # self.api_common.add_meta_status(table_id=main_id, type_name="heatmap")

    def send_files(self):
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

    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        task_info = self.api.api('task_info.metaasv_task_info')
        task_info.add_task_info()
        task_info.update_sg_task(self._sheet.id, self.option("database_type"), self.option("database"))
        self.logger.info("<<<<<<<<<<<<<<<<<<<")
        self.logger.info(self.work_dir)
        gevent.spawn_later(5, self.end)
        super(MetaAsvApiWorkflow, self).run()