# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifiy = modified 2018.10.23

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class AnnoTcdbWorkflow(Workflow):
    """
    宏基因组tcdb注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoTcdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_tcdb_anno", "type": "infile", "format": "sequence.profile_table"},
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
            # {"name": "submit_location", "type": "string"},
            {"name": "main_table_id", "type": "string"},
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
            self.abu_act = "T"
        if self.option("identity") != 0 or self.option("align_length") != 0 or self.option("level_select") != "all":
            self.fun_select = "T"

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_tcdb workflow")
        self.run_filter()
        super(AnnoTcdbWorkflow, self).run()

    def run_filter(self):
        self.logger.info("start run anno_select module!")
        options = {
            "gene_anno": self.option("gene_tcdb_anno"),
            "database": "tcdb",
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
            "gene_filter": self.profile
        }
        self.filter.set_options(options)
        self.filter.on('end', self.set_db)
        self.filter.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno_tcdb = self.api.api("metagenomic.anno_tcdb")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_tcdb_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        tcdb_profile_dir = self.filter.output_dir
        all_files = os.listdir(tcdb_profile_dir)
        for i in all_files:
            old = os.path.join(tcdb_profile_dir, i)
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
        sanger_type, sanger_path = self._sheet.output.split(':')
        sanger_prefix = Config().get_netdata_config(sanger_type)
        self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        web_path = self.remote_dir + "/gene_tcdb_anno.xls"
        api_anno_tcdb.add_anno_tcdb(geneset_id, specimen, web_path, main=False,
                                    main_table_id=main_id)
        api_anno_tcdb.add_anno_tcdb_arg(main_id, tcdb_profile_dir, update_main=False)
        api_anno_tcdb.add_anno_tcdb_type(main_id, tcdb_profile_dir)
        api_anno_tcdb.add_anno_tcdb_class(main_id, tcdb_profile_dir)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TCDB功能注释结果目录", 0, "120084"],
            ["gene_tcdb_anno.xls", "xls", "每条基因的TCDB功能注释表", 0, "120087"],
            ["tcdb_ARG_profile.xls", "xls", "各样品ARDB ARG丰度表", 0, "120088"],
            ["tcdb_type_profile.xls", "xls", "各样品ARDB Type丰度表", 0, "120089"],
            ["gene_tcdb_class_stat.xls", "xls", "ARDB Class基因信息统计表", 0, "120090"],
            ["tcdb_class_profile.xls", "xls", "各样品ARDB Class丰度表", 0, "120086"]
        ])
        super(AnnoTcdbWorkflow, self).end()

