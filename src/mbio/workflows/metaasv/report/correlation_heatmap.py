# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
from mbio.packages.metaasv.common_function import link_dir, calculate_abundance
import re


class CorrelationHeatmapWorkflow(Workflow):
    """
    metaasv 相关性heatmap图分析
    环境因子与物种丰度之间相关性
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(CorrelationHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的ASVfile
            {"name": "asv_id", "type": "string"},#asv表ID
            {"name": "update_info", "type": "string"},
            {"name": "env_file", "type": "infile", 'format': "meta.otu.group_table"},  # 输入的env file
            {"name": "env_id", "type": "string"},##env主表ID
            {"name": "main_id", "type": "string"},##主表ID
            {"name": "env_labs", "type": "string"},##筛选的环境因子
            {"name": "level", "type": "int"},##分类学水平
            {"name": "group_id", "type": "string"},
            # {"name": "submit_location", "type": "string", "default": "cor_heatmap"},
            # {"name": "task_type", "type": "int", "default": 2},
            {"name": "method", "type": "string", "default": "pearsonr"},##相关性系数
            {"name": "env_cluster", "type": "string", "default": ""},#环境因子聚类方式
            {"name": "species_cluster", "type": "string", "default": ""},##物种聚类方式
            {"name": "group_detail", "type": "string"},##
            {"name": "top_species", "type": "int", "default": 0},##top物种数
            {"name": "env_distance", "type": "string", "default": "bray_curtis"},##环境因子计算距离方式
            {"name": "species_distance", "type": "string", "default": "bray_curtis"},##物种计算距离方式
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.correlation = self.add_module('metaasv.pearsons_correlation')
        self.params = {}
        self.name_to_name = {}
        self.env_name = {}

    def run_correlation(self):
        """
        计算相关性
        :return:
        """
        options = {
            'otutable': self.option('otu_file'),
            'envtable': self.option('env_file'),
            "method": self.option('method'),
            "env_cluster": self.option("env_cluster"),
            "species_cluster": self.option("species_cluster"),
            "top_species": self.option('top_species'),
            "species_distance": self.option("species_distance"),
            "env_distance": self.option("env_distance"),
            }
        self.correlation.set_options(options)
        self.correlation.on("end", self.set_db)
        self.correlation.run()

    def run(self):
        self.run_correlation()
        super(CorrelationHeatmapWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        link_dir(self.correlation.output_dir, self.output_dir)
        ###测试用
        if os.path.exists(os.path.join(self.work_dir, "abundance.absolute.xls")):
            os.remove(os.path.join(self.work_dir, "abundance.absolute.xls"))
        os.link(self.option("otu_file").prop['path'], os.path.join(self.work_dir, "abundance.absolute.xls"))
        self.logger.info("连接结果文件成功！")
        calculate_abundance(os.path.join(self.work_dir, "abundance.absolute.xls"), abundance_method="absolute")
        new_species_tree = ""
        env_tree = ""
        new_env_tree = ""
        env_list = []
        species_list = []
        api_correlation = self.api.api("metaasv.correlation_heatmap")
        corr_path = glob.glob(self.output_dir+"/*correlation*")
        pvalue_path = glob.glob(self.output_dir+"/*pvalue*")

        env_tree_path = self.output_dir + "/env_tree.tre"
        species_tree_path = self.output_dir + "/species_tree.tre"

        if os.path.exists(env_tree_path):
            with open(env_tree_path, "r") as f:
                env_tree = f.readline().strip()
                env_list_all = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', env_tree)
                env_list = [i[1] for i in env_list_all]
        if os.path.exists(species_tree_path):
            with open(species_tree_path, "r") as f:
                species_tree = f.readline().strip()
                species_list_all = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                species_list = [i[1] for i in species_list_all]

        main_id = ObjectId(self.option("main_id"))
        api_correlation.add_correlation_detail(corr_path[0], "correlation", main_id)
        api_correlation.add_correlation_detail(pvalue_path[0], "pvalue", main_id, env_list=env_list, species_list=species_list)
        if os.path.exists(env_tree_path):
            api_correlation.insert_tree_table(env_tree_path,main_id, "specimen")
        if os.path.exists(species_tree_path):
            api_correlation.insert_tree_table(species_tree_path,main_id, "species")
        if os.path.exists(self.correlation.output_dir + "/pearsons_correlation.xls"):
            os.rename(self.correlation.output_dir + "/pearsons_correlation.xls",self.correlation.output_dir + "/" + self.option('method') + "_correlation.xls")
        if os.path.exists(self.correlation.output_dir + "/pearsons_pvalue.xls"):
            os.rename(self.correlation.output_dir + "/pearsons_pvalue.xls",self.correlation.output_dir + "/" + self.option('method') + "_pvalue.xls")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.correlation.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性Heatmap分析结果目录", 0, ""],   # add 2 lines by hongdongxuan 20170324
            ["./pearsons_correlation.xls", "xls", "相关性系数表", 0, ""],
            ["./pearsons_pvalue.xls", "xls", "相关性P值", 0, ""],
            ["./species_tree.tre", "tre", "物种层级聚类树", 0, ""],
            ["./env_tree.tre", "tre", "环境因子层级聚类树", 0, ""]
        ])
        # print self.get_upload_files()
        super(CorrelationHeatmapWorkflow, self).end()
