# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# last_modifiy = modified 2018.10.31


from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class AnnoQsWorkflow(Workflow):
    """
    宏基因组QS注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoQsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_qs_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string"},
            {"name": "identity", "type": "float", "default": 0.0},
            {"name": "align_length", "type": "int", "default": 0},
            {"name": "level_select", "type": "string", "default": "all"},
            {"name": "abu_num", "type": "string", "default": "all"},
            {"name": "abu_proportion", "type": "string", "default": "all"},
            {"name": "gene_type", "type": "string"},
            {"name": "group", "type": "int", "default": 1},  # 1 所有样本，2 剔除样本
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "align_database", "type": "string", "default": "all"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filter = self.add_module("annotation.mg_fun_select")
        self.profile = "T"
        self.abu_act = "F"
        self.fun_select = "F"
        if self.option("group") == 1 and self.option("gene_type") == "Origin":
            self.profile = "F"
        if self.option("abu_num") != "all" or self.option("abu_proportion") != "all":
            self.abu_act = "T"
        if self.option("identity") != 0 or self.option("align_length") != 0 or self.option("level_select") != "all":
            self.fun_select = "T"

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_QS workflow")
        self.run_filter()
        super(AnnoQsWorkflow, self).run()

    def run_filter(self):
        self.logger.info("start run anno_select module!")
        options = {
            "gene_anno": self.option("gene_qs_anno"),
            "database": "qs",
            "geneset_table": self.option("geneset_table"),
            "lowest_level_profile": self.option("lowest_level_profile"),
            "gene_profile": self.option("gene_profile"),
            "samples": self.option("samples"),
            "identity": self.option("identity"),
            "align_length": self.option("align_length"),
            "level_select": self.option("level_select"),
            "abu_num": self.option("abu_num"),
            "abu_proportion": self.option("abu_proportion"),
            "abu_filter": self.abu_act ,
            "fun_filter": self.fun_select,
            "gene_filter": self.profile,
            "group_table":self.option("group_table")
        }
        self.filter.set_options(options)
        self.filter.on('end', self.set_db)
        self.filter.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno_qs = self.api.api("metagenomic.qs_anno")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        self.logger.info(old_anno_file)
        anno_file_path = self.output_dir + "/gene_qs_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        qs_profile_dir = self.filter.output_dir
        all_files = os.listdir(qs_profile_dir)
        for i in all_files:
            old = os.path.join(qs_profile_dir, i)
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12804701")
        self.logger.info(main_id)
        self.remote_dir = self._sheet.output
        web_path = self.remote_dir + "gene_qs_anno.xls"
        api_anno_qs.add_anno_qs(geneset_id,specimen,web_path,main=False, main_table_id=main_id)
        api_anno_qs.add_anno_qs_class(main_id,qs_profile_dir)
        # api_anno_qs.add_qs_graph(main_id, qs_profile_dir + '/anno_qs_graph.xls')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "QS功能注释结果目录", 0, "120058"],
            ["qs_class_profile.xls", "xls", "各样品QS的class水平丰度表",0,"120266"],
            ["gene_qs_anno.xls", "xls", "各样品gene对应注释表",0,"120267"],
            ["qs_lowest_profile.xls", "xls", "各样品QS的最低水平丰度表",0,"120268"],
            # ["anno_qs_graph.xls", "xls", "各样品QS画图的数据表",0,"120269"]
        ])
        super(AnnoQsWorkflow, self).end()
