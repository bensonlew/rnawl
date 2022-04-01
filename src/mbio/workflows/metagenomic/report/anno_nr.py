# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2017.10.11

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class AnnoNrWorkflow(Workflow):
    """
    宏基因组nr注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoNrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_nr_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "string", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string"},
            {"name": "identity", "type": "float", "default": 0},
            {"name": "align_length", "type": "float", "default": 0},
            {"name": "level_select", "type": "string", "default": "all"},
            {"name": "abu_num", "type": "string", "default": "all"},
            {"name": "abu_proportion", "type": "string", "default": "all"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "gene_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "group", "type": "int", "default": 1},
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
        self.logger.info("start run anno_nr workflow")
        self.run_filter()
        super(AnnoNrWorkflow, self).run()

    def run_filter(self):
        self.logger.info("start run anno_select module!")
        options = {
            "gene_anno": self.option("gene_nr_anno"),
            "database": "nr",
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
        self.logger.info("start nr database")
        api_anno_nr = self.api.api("metagenomic.mg_anno_nr")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        self.logger.info(self.option("gene_nr_anno"))
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_nr_anno.xls"
        self.logger.info(anno_file_path)
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        nr_profile_dir = self.filter.output_dir
        all_files = os.listdir(nr_profile_dir)
        for i in all_files:
            old = os.path.join(nr_profile_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.logger.info("正在写入mongo数据库")
        self.logger.info(nr_profile_dir)
        main_id = self.option("main_table_id")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir, "gene_nr_anno.xls")
        api_anno_nr.add_anno_nr(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_nr.add_anno_nr_detail(main_id, nr_profile_dir, update_main=False)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "NR功能注释结果目录", 0, "120035"],
            ["gene_nr_anno.xls", "xls", "每条基因的物种注释表", 0, "120036"],
            ["tax_d.xls", "xls", "域注释丰度表", 0, "120038"],
            ["tax_k.xls", "xls", "界注释丰度表", 0, "120039"],
            ["tax_p.xls", "xls", "门注释丰度表", 0, "120040"],
            ["tax_c.xls", "xls", "纲注释丰度表", 0, "120041"],
            ["tax_o.xls", "xls", "目丰注释度表", 0, "120042"],
            ["tax_f.xls", "xls", "科注释丰度表", 0, "120043"],
            ["tax_g.xls", "xls", "属注释丰度表", 0, "120044"],
            ["tax_s.xls", "xls", "种注释丰度表", 0, "120045"],
        ])
        super(AnnoNrWorkflow, self).end()
