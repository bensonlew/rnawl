# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from mbio.packages.metaasv.common_function import link_dir,link_file


class NetworkCorfdWorkflow(Workflow):
    """
    metaasv 双因素相关性网络
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NetworkCorfdWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "sequence.profile_table"},  # 物种注释表或功能注释表
            {"name": "level", "type": "int", "default": 9},
            {"name": "asv_id", "type": 'string'},
            {"name": "env_id","type":"string"},
            {"name": "env_file", "type": "infile", "format": "sequence.profile_table"},  # 环境因子表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "env_labs", "type": "string"},  # 环境因子
            {"name": "fac1_top", "type": "int", "default": 0},  # 总丰度top
            {"name": "pvalue", "type": "float", "default": 1},##此参数页面不展示
            {"name": "coefficient", "type": "string", "default": "spearman"},  # 相关性系数spearman,pearson,kendall
            {"name": "coefficient_value", "type": "float", "default": 0.0},   # 相关性阈值，此参数页面不展示
            {"name": "update_info", "type": "string"},
            {"name": "color_level", "type": "string"}, # 颜色显示水平
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_tax_samples = self.add_tool("metaasv.sort_samples")
        self.get_correlation = self.add_tool('meta.association_model.correlation_fd')  # 计算相关性
        self.network_cor_tool = self.add_tool("meta.association_model.network_cor")  # 网络相关系数计算

    def run(self):
        self.logger.info("start NetworkCorfd!")
        self.sort_tax_samples.on("end", self.run_cal_correlation)
        self.run_tax_sort_samples()
        super(NetworkCorfdWorkflow, self).run()

    def run_tax_sort_samples(self):
        """
        功能：排序、选top物种、分组合并
        :return:
        """
        abund_table = self.option("otu_table").path
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group_table').path,
            "top" : self.option('fac1_top'),
            "fill_zero": "false",
            "abundance_method": "all"
        })
        self.sort_tax_samples.run()

    def run_cal_correlation(self):
        """
        计算相关性
        :return:
        """
        tax_abund_table =  self.sort_tax_samples.option("out_otu_table")
        self.logger.info("start cal_correlation")
        options = ({
            'coefficient': self.option('coefficient'),
            "coefficient_value": self.option('coefficient_value'),
            'p_value': self.option('pvalue')
        })
        self.profile_table1 = tax_abund_table  #self.get_abu_table.option("outprofile")
        self.profile_table2 = self.option("env_file")
        options["trans_t2"] = True
        options["table1"] = self.profile_table1.path
        options["table2"] = self.profile_table2
        self.get_correlation.set_options(options)
        self.get_correlation.on("end", self.run_network)
        self.get_correlation.run()

    def run_network(self):
        """
        计算网络相关性系数
        :return:
        """
        self.logger.info("开始计算网络相关系数!")
        correlation_file = self.get_correlation.option("correlation_file")
        self.logger.info("correlation_file_path:{}".format(correlation_file.prop['path']))
        options = {
            "correlation_file": correlation_file,
        }
        self.network_cor_tool.set_options(options)
        self.network_cor_tool.on('end', self.set_db)
        self.network_cor_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        self.logger.info("start link")
        link_dir(self.network_cor_tool.output_dir, self.output_dir)
        self.logger.info("start setoutput")
        api_run_network = self.api.api("metaasv.network_corfd")
        cor_file = self.get_correlation.option("correlation_file").prop["path"]
        name = self.option("coefficient")
        link_cor = self.output_dir + "/" + name + "_corr_edge.txt"
        link_file(cor_file, link_cor)

        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        self.logger.info(main_id)
        attributes_file = self.output_dir + "/corr_network_attributes.txt"
        profile1 = self.profile_table1.prop["path"]
        profile2 = self.profile_table2.prop["path"]
        api_run_network.add_network_corfd(attributes_file, main=False, main_table_id=main_id)
        api_run_network.add_network_corfd_link(main_id, self.output_dir)
        api_run_network.add_network_corfd_degree(main_id, self.output_dir)
        api_run_network.add_network_corfd_node(main_id, self.output_dir, profile1, profile2,int(self.option('color_level')) ) # ,level1=level1, level2=level2
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性网络分析结果目录",0,""],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表",0,""],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件",0,""],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表",0,""],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表",0,""],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表",0,""],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表",0,""],
        ])
        regexps = ([
               [r".*_corr_edge\.txt", "txt", "物种或功能相似性网络边文件",0,""]
            ])
        result_dir.add_regexp_rules(regexps)
        super(NetworkCorfdWorkflow, self).end()
