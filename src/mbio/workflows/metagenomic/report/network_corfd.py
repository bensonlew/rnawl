# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2017.12.27

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class NetworkCorfdWorkflow(Workflow):
    """
    宏基因组NetworkCorfd注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NetworkCorfdWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "infile", "format": "sequence.profile_table"},  # 基因list
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},  # gene丰度表
            {"name": "fac1_anno", "type": "infile", "format": "sequence.profile_table"},  # 物种注释表或功能注释表
            {"name": "fac2_anno", "type": "infile", "format": "sequence.profile_table"},  # 物种注释表或功能注释表
            {"name": "env_file", "type": "infile", "format": "sequence.profile_table"},  # 环境因子表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            #{"name": "samples", "type": "string"},
            {"name": "fac1_database", "type": "string"},  # 物种选择水平
            {"name": "fac2_database", "type": "string"},  # 物种选择水平
            {"name": "fac1_level", "type": "string"},  # 物种选择水平
            {"name": "fac2_level", "type": "string"},  # 功能选择水平
            {"name": "env_id", "type": "string"},  # 环境因子
            {"name": "env_labs", "type": "string"},  # 环境因子
            {"name": "fac1_top", "type": "string", "default": "10"},  # 总丰度top
            {"name": "fac2_top", "type": "string", "default": "10"},  # 总丰度top
            {"name": "pvalue", "type": "float", "default": 0.05},
            {"name": "coefficient", "type": "string", "default": "spearman"},  # 相关性系数spearman,pearson,kendall
            {"name": "coefficient_value", "type": "float", "default": 0.5},   # 相关性阈值
            {"name": "update_info", "type": "string"},
            {"name": "fac1_color_level", "type": "string"},             # 颜色显示水平
            {"name": "fac2_color_level", "type": "string"},             # 颜色显示水平
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_table_id", "type": "string"},
            {"name": "fac1_second_level", "type": "string"},
            {"name": "fac1_lowestlevel", "type": "string"},
            {"name": "fac2_second_level", "type": "string"},
            {"name": "fac2_lowestlevel", "type": "string"},
            {"name": "clean_stat", "type": "infile", "format": "meta.profile"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "method", "type": "string", "default": ""}, # geneset_table是什么丰度计算方法
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.get_taxon_table = self.add_tool('meta.association_model.creat_level_table')  # 创建新level丰度表并筛选top
        #self.get_fun_table = self.add_tool('meta.association_model.creat_level_table')  # 创建新level丰度表并筛选top
        self.get_correlation = self.add_tool('meta.association_model.correlation_fd')  # 计算相关性
        self.network_cor_tool = self.add_tool("meta.association_model.network_cor")  # 网络相关系数计算

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start NetworkCorfd!")
        if self.option("fac1_anno").is_set and self.option("fac2_anno").is_set:
            self.get_abu_table = self.add_tool('meta.association_model.creat_level_table')
            self.get_abu_table_fact = self.add_tool('meta.association_model.creat_level_table')
            self.on_rely([self.get_abu_table, self.get_abu_table_fact], self.run_cal_correlation)
            self.run_get_abu_table()
            self.run_get_abu_table_fact()
        if self.option("fac1_anno").is_set and self.option("env_file").is_set:
            self.get_abu_table = self.add_tool('meta.association_model.creat_level_table')
            self.get_abu_table.on("end", self.run_cal_correlation)
            self.run_get_abu_table()
        super(NetworkCorfdWorkflow, self).run()

    def run_get_abu_table(self):
        options = {
            'anno_file': self.option('fac1_anno'),
            'gene_profile': self.option('gene_profile'),
            'gene_list': self.option('gene_list'),
            'level': self.option('fac1_level'),
            'top': self.option('fac1_top'),  ### 为与all共用改成str类型
            "group_table": self.option('group_table'),
            "Total": 0,
            "database": self.option("fac1_database")
        }
        if self.option("fac1_second_level"):
            options["levelname"] = self.option('fac1_second_level')
        if self.option("fac1_lowestlevel"):
            options["lowestlevel"] = self.option('fac1_lowestlevel')
        if self.option("fac1_color_level"):
            options["colorlevel"] = self.option('fac1_color_level')
        self.get_abu_table.set_options(options)
        self.get_abu_table.run()

    def run_get_abu_table_fact(self):
        options = {
            'anno_file': self.option('fac2_anno'),
            'gene_profile': self.option('gene_profile'),
            'gene_list': self.option('gene_list'),
            'level': self.option('fac2_level'),
            'top': self.option('fac2_top'),  ### 为与all共用改成str类型
            "group_table": self.option('group_table'),
            "Total": 0,
            "database": self.option("fac2_database")
        }
        if self.option("fac2_second_level"):
            options["levelname"] = self.option('fac2_second_level')
        if self.option("fac2_lowestlevel"):
            options["lowestlevel"] = self.option('fac2_lowestlevel')
        if self.option("fac2_color_level"):
            options["colorlevel"] = self.option('fac2_color_level')
        self.get_abu_table_fact.set_options(options)
        self.get_abu_table_fact.run()

    def run_cal_correlation(self):
        self.logger.info("start cal_correlation")
        options = ({
            #'table1': profile_table1,
            #'table2': profile_table2,
            'coefficient': self.option('coefficient'),
            "coefficient_value": self.option('coefficient_value'),
            'p_value': self.option('pvalue')
        })
        if self.option("fac1_anno").is_set and self.option("fac2_anno").is_set:
            self.profile_table1 = self.get_abu_table.option("outprofile")
            self.profile_table2 = self.get_abu_table_fact.option("outprofile")
            self.level1 = self.option("fac1_level")
            self.level2 = self.option("fac2_level")
        elif self.option("fac1_anno").is_set and self.option("env_file").is_set:
            self.profile_table1 = self.get_abu_table.option("outprofile")
            self.profile_table2 = self.option("env_file")
            options["trans_t2"] = True
            self.level1 = self.option("fac1_level")
            self.level2 = "env"
        options["table1"] = self.profile_table1
        options["table2"] = self.profile_table2
        self.get_correlation.set_options(options)
        self.get_correlation.on("end", self.run_network)
        self.get_correlation.run()

    def run_network(self):
        self.logger.info("开始计算网络相关系数!")
        correlation_file = self.get_correlation.option("correlation_file")
        self.logger.info(correlation_file)
        self.logger.info(correlation_file.prop["path"])
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
        self.logger.info("start setoutput")
        api_run_network = self.api.api("metagenomic.network_corfd")
        cor_file = self.get_correlation.option("correlation_file").prop["path"]
        method = self.option("coefficient")
        #name = self.option("level") + "." + self.option("coefficient")
        name =  self.option("coefficient")
        link_cor = self.output_dir + "/" + name +  "_corr_edge.txt"
        if os.path.exists(link_cor):
            os.remove(link_cor)
        os.link(cor_file, link_cor)
        all_files = os.listdir(self.network_cor_tool.output_dir)
        self.logger.info("start link")
        for i in all_files:
            old = os.path.join(self.network_cor_tool.output_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        self.logger.info(main_id)
        attributes_file = self.output_dir + "/corr_network_attributes.txt"
        profile1 = self.profile_table1.prop["path"]
        profile2 = self.profile_table2.prop["path"]
        level1 = self.level1
        level2 = self.level2
        if self.option("fac2_second_level"):
            fac2_second_level = self.option("fac2_second_level")
        api_run_network.add_network_corfd(attributes_file, main=False, main_table_id=main_id)
        api_run_network.add_network_corfd_link(main_id, self.output_dir)
        api_run_network.add_network_corfd_degree(main_id, self.output_dir)
        api_run_network.add_network_corfd_node(main_id, self.output_dir, profile1, profile2,level1=level1, level2=level2,fac2_second_level=fac2_second_level)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "network_corfd")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "network_corfd")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    @envname_restore
    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性网络分析结果目录", 0, ""],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表", 0, ""],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件", 0, ""],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表", 0, ""],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表", 0, ""],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表", 0, ""],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表", 0, ""],
            ["CorrNetwork.pdf", "pdf", "双因素相关性网络图"]
        ])
        regexps = ([
               [r".*_corr_edge\.txt", "txt", "物种或功能相似性网络边文件", 0, ""]
            ])
        result_dir.add_regexp_rules(regexps)
        super(NetworkCorfdWorkflow, self).end()
