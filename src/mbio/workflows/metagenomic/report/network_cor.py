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
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class NetworkCorWorkflow(Workflow):
    """
    宏基因组NetworkCor注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NetworkCorWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "infile", "format": "sequence.profile_table"},  # 基因list
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},  # gene丰度表
            {"name": "anno_file", "type": "infile", "format": "sequence.profile_table"},  # 注释表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            #{"name": "samples", "type": "string"},
            {"name": "level", "type": "string"},  # 选择水平
            {"name": "top", "type": "int", "default": 10},  # 总丰度top
            {"name": "pvalue", "type": "float", "default": 0.05},
            {"name": "coefficient", "type": "string", "default": "spearman"},  # 相关性系数spearman,pearson,kendall
            {"name": "coefficient_value", "type": "float", "default": 0.5},   # 相关性阈值
            {"name": "update_info", "type": "string"},
            {"name": "anno_type", "type": "string"},               # 选择数据库nr、kegg、cog、vfdb、card、ardb、cazy
            {"name": "color_level", "type": "string"},             # 颜色显示水平
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_table_id", "type": "string"},
            {"name": "second_level", "type": "string"},
            {"name": "lowestlevel", "type": "string"},
            {"name": "clean_stat", "type": "infile", "format": "meta.profile"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "method", "type": "string", "default": ""}, # geneset_table是什么丰度计算方法
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.get_abund_table = self.add_tool('meta.association_model.creat_level_table')  # 创建新level丰度表并筛选top
        self.get_correlation = self.add_tool('meta.association_model.correlation')  # 计算相关性
        self.network_cor_tool = self.add_tool("meta.association_model.network_cor")  # 网络相关系数计算

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start NetworkCor!")
        self.run_get_abund_table()
        super(NetworkCorWorkflow, self).run()

    def run_get_abund_table(self):
        self.logger.info(self.option('level'))
        if self.option("anno_type") == "go":
            self.get_abund_table = self.add_module("annotation.anno_go_stat")
            options = {
                'gene_anno': self.option('anno_file').prop["path"],
                'go1234level_out': self.option('go1234level_out'),
                'reads_profile_table': self.option('gene_profile').prop["path"],
                'level': self.option('level'),
            }
        else:
            options = {
                'anno_file': self.option('anno_file'),
                'gene_profile': self.option('gene_profile'),
                'gene_list': self.option('gene_list'),
                'level': self.option('level'),
                'top': str(self.option('top')),  ### 为与all共用改成str类型
                "group_table": self.option('group_table'),
                "Total": 0,
                "database": self.option("anno_type"),

            }
            if self.option("anno_type") == "nr":
                options["anno_type"] = "nr"
            if self.option("second_level"):
                options["levelname"] = self.option('second_level')
            if self.option("lowestlevel"):
                options["lowestlevel"] = self.option('lowestlevel')
            if self.option("color_level"):
                options["colorlevel"] = self.option('color_level')
            if self.option("method") == "ppm" and self.option("clean_stat").is_set:
                options["method"] = self.option("method")
                options["clean_stat"] = self.option("clean_stat")
        self.get_abund_table.set_options(options)
        self.get_abund_table.on("end", self.run_cal_correlation)
        self.get_abund_table.run()

    def run_cal_correlation(self):
        self.logger.info("start cal_correlation")
        if self.option("anno_type") == "go":
            profile_table = self.get_abund_table.option("out_table").prop['path']
        else:
            profile_table = self.get_abund_table.option("outprofile").prop['path']
        self.get_correlation.set_options({
            'profile_table': profile_table,
            'coefficient': self.option('coefficient'),
            "coefficient_value": self.option('coefficient_value'),
            'p_value': self.option('pvalue')
        })
        self.get_correlation.on("end", self.run_networkcor)
        self.get_correlation.run()

    def run_networkcor(self):
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
        api_networkcor = self.api.api("metagenomic.network_cor")
        cor_file = self.get_correlation.option("correlation_file").prop["path"]
        method = self.option("coefficient")
        name = self.option("level") + "." + self.option("coefficient")
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12802201")
        self.logger.info(main_id)
        attributes_file = self.output_dir + "/corr_network_attributes.txt"
        if self.option("anno_type") == "go":
            profile = self.get_abund_table.option("out_table").prop['path']
        else:
            profile = self.get_abund_table.option("outprofile").prop['path']
        anno_type = self.option("anno_type")
        if self.option("second_level"):
            level = self.option("second_level")
            self.option("color_level",None)
        else:
            level = self.option("level")
        api_networkcor.add_network_cor(anno_type, attributes_file, main=False, main_table_id=main_id)
        api_networkcor.add_network_cor_link(main_id, self.output_dir, anno_type)
        api_networkcor.add_network_cor_degree(main_id, self.output_dir, anno_type=anno_type)
        if self.option("color_level"):
            color_level = self.option("color_level")
            api_networkcor.add_network_cor_node(main_id, profile, level, self.output_dir, anno_type=anno_type,
                                                color_level=color_level)
        else:
            api_networkcor.add_network_cor_node(main_id, profile, level, self.output_dir, anno_type=anno_type)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "network_cor")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "network_cor")
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

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性网络分析结果目录", 0, "120223"],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表", 0, "120224"],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件", 0, "120225"],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表", 0, "120226"],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表", 0, "120227"],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表", 0, "120228"],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表", 0, "120229"],
            ["CorrNetwork.pdf", "pdf", "单因素相关性网络图"]
        ])
        regexps = ([
               [r".*_corr_edge\.txt", "txt", "物种或功能相似性网络边文件", 0, "120230"]
            ])
        result_dir.add_regexp_rules(regexps)
        super(NetworkCorWorkflow, self).end()
