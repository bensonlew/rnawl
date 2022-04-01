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
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class NetworkWorkflow(CommTableWorkflow):
    """
    宏基因组分布网络注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "infile", "format": "sequence.profile_table"},  # 基因list
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},  # gene丰度表
            {"name": "anno_table", "type": "infile", "format": "sequence.profile_table"},  # 注释表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "level", "type": "string"},  # 选择水平
            {"name": "anno_type", "type": "string"},  # 选择数据库nr、kegg、cog、vfdb、card、ardb、cazy
            {"name": "group_method", "type": "int"},  # 分组求和方法：无：0，求和：1，均值：2，中位数：3
            {"name": "top", "type": "int", "default": 10},  # 总丰度top
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_table_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.get_abund_table = self.add_tool('meta.association_model.creat_level_table')  # 丰度计算和筛选
        self.network_tool = self.add_tool("meta.association_model.network")  # 网络系数计算

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start Network!")
        if self.option("anno_type") == "go":
            if self.option("method") == "ppm" and self.option("clean_stat").is_set:
                self.run_cal_ppm(rely=self.run_go_abu(is_run=True), is_run=True)
            else:
                self.run_go_abu(rely=self.run_network, is_run=True)
        else:
            if self.option("method") == "ppm" and self.option("clean_stat").is_set:
                self.run_cal_ppm(rely=self.run_get_abund_table, is_run=True)
            else:
                self.run_get_abund_table()
        super(NetworkWorkflow, self).run()

    def run_get_abund_table(self):
        self.logger.info("start get abund table>>>>>>>>>>>")
        self.get_abund_table.set_options({
            'anno_file': self.option('anno_table'),
            'gene_profile': self.option('gene_profile'),
            'gene_list': self.option('gene_list'),
            'level': self.option('level'),
            'top': str(self.option('top')),
            'group_method': self.option('group_method'),
            "group_table": self.option('group_table'),
            "Total": 0
        })
        self.get_abund_table.on("end", self.run_network)
        self.get_abund_table.run()

    def run_network(self):
        if self.option("anno_type") == "go":
            self.table = self.abundance.option("out_table").prop['path']
        else:
            self.table = self.get_abund_table.option("outprofile").prop["path"]
        options = {
            "profile_table": self.table,
        }
        self.network_tool.set_options(options)
        self.network_tool.on('end', self.set_db)
        self.network_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_network = self.api.api("metagenomic.network")
        all_files = os.listdir(self.network_tool.output_dir)
        #abu_file = self.get_abund_table.option("outprofile").prop["path"]
        for i in all_files:
            if "real_dc" in i:
                n = i.replace("real_dc", "network")
            elif "real" in i:
                n = i.replace("real", "network")
            else:
                n = i
            old = os.path.join(self.network_tool.output_dir, i)
            link = os.path.join(self.output_dir, n)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        anno_type = self.option("anno_type")
        group_file = self.option("group_table").prop["path"]
        level = self.option("level")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12802101")
        self.logger.info(main_id)
        attributes_file = self.output_dir + "/network_attributes.txt"
        api_network.add_network(anno_type, attributes_file, main=False, main_table_id=main_id)
        api_network.add_network_node(main_id, group_file, self.output_dir, anno_type, level, self.table)
        api_network.add_network_link(main_id, self.output_dir, anno_type)
        api_network.add_network_degree(main_id, self.output_dir, anno_type)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "network")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "network")
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
            [".", "", "分布网络分析结果目录", 0, "120214"],
            ["network_attributes.txt", "txt", "网络单值属性表", 0, "120215"],
            ["network_centrality.txt", "txt", "网络中心系数表", 0, "120216"],
            ["network_degree.txt", "txt", "网络度统计总表", 0, "120217"],
            ["network_tax_fun_degree.txt", "txt", "网络物种或功能节点度分布表", 0, "120218"],
            ["network_sample_degree.txt", "txt", "网络sample节点度分布表", 0, "120219"],
            ["network_nodes_degree.txt", "txt", "网络所有节点度分布表", 0, "120220"],
            ["network_edge_table.txt", "txt", "网络边的属性表", 0, "120221"],
            ["network_nodes_table.txt", "txt", "网络节点属性表", 0, "120222"],
            ["network.pdf", "pdf", "物种/功能分布网络图"]
        ])
        super(NetworkWorkflow, self).end()
