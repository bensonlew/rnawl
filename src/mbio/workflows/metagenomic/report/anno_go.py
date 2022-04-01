# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# last_modifiy = modified 2018.11.14

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

class AnnoGoWorkflow(Workflow):
    """
    宏基因组go注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoGoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_go_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string"},
            {"name": "level_select", "type": "string", "default": "all"},
            {"name": "abu_num", "type": "string", "default": "all"},
            {"name": "abu_proportion", "type": "string", "default": "all"},
            {"name": "gene_type", "type": "string"},
            {"name": "group", "type": "int", "default": 1},  # 1 所有样本，2 剔除样本
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0}
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
        if self.option("level_select") != "all":
            self.fun_select = "T"

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_GO workflow")
        self.run_filter()
        super(AnnoGoWorkflow, self).run()

    def run_filter(self):
        self.logger.info("start run anno_select module!")
        options = {
            "gene_anno": self.option("gene_go_anno"),
            "database": "go",
            "geneset_table": self.option("geneset_table"),
            "lowest_level_profile": self.option("lowest_level_profile"),
            "go1234level_out":self.option("go1234level_out"),
            "gene_profile": self.option("gene_profile"),
            "samples": self.option("samples"),
            "level_select": self.option("level_select"),
            "abu_num": self.option("abu_num"),
            "abu_proportion": self.option("abu_proportion"),
            "abu_filter": self.abu_act,
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
        api_anno_go = self.api.api("metagenomic.anno_go")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        self.logger.info(old_anno_file)
        anno_file_path = self.output_dir + "/all.go.annotation.xls"
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12804401")
        self.logger.info(main_id)
        self.remote_dir = self._sheet.output
        web_path = self.remote_dir + "all.go.annotation.xls"
        api_anno_go.add_anno_go(geneset_id,specimen,web_path, main=False, main_table_id=main_id)
        api_anno_go.add_go_func(main_id,qs_profile_dir + '/all.go1.function.xls',59)
        api_anno_go.add_go_func(main_id, qs_profile_dir + '/all.go12.function.xls',60)
        api_anno_go.add_go_func(main_id, qs_profile_dir + '/all.go123.function.xls',61)
        api_anno_go.add_go_func(main_id, qs_profile_dir + '/all.go1234.function.xls',62)
        if self.option("save_pdf"):
        # self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        # if self.pdf_status:
            name = get_name(self.option("main_table_id"), "anno_go")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "anno_go")
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
        if self.option("save_pdf"):
        # if self.pdf_status:
            # if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "QS功能注释结果目录",0,"120241"],
            ["all.go1.function.xls", "xls", "level1水平各样品丰度表",0,"120242"],
            ["all.go12.function.xls", "xls", "level2水平各样品丰度表",0,"120243"],
            ["all.go123.function.xls", "xls", "level3水平各样品丰度表",0,"120244"],
            ["all.go1234.function.xls", "xls", "level4水平各样品丰度表",0,"120245"],
            ["all.go.annotation.xls", "xls", "各样品go注释gene对应注释表",0,"120246"],
            ["go_level_bar.pdf", "pdf", "GO注释分类统计柱形图"],
        ])
        super(AnnoGoWorkflow, self).end()