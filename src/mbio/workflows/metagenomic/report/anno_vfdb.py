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
from comfun import ComfunWorkflow

class AnnoVfdbWorkflow(ComfunWorkflow):
    """
    宏基因组vfdb注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoVfdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_vfdb workflow")
        self.run_filter(self.set_db)
        super(AnnoVfdbWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno_vfdb = self.api.api("metagenomic.mg_anno_vfdb")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_vfdb_anno.xls"
        self.logger.info(anno_file_path)
        self.logger.info(old_anno_file)
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        vfdb_profile_dir = self.filter.output_dir
        all_files = os.listdir(vfdb_profile_dir)
        for i in all_files:
            old = os.path.join(vfdb_profile_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.logger.info("正在写入mongo数据库")
        self.logger.info(vfdb_profile_dir)
        self.logger.info(self.option("align_database"))
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12800901")
        self.logger.info(main_id)
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir,"gene_vfdb_anno.xls")
        api_anno_vfdb.add_anno_vfdb(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        if self.option("align_database") == "core":
            api_anno_vfdb.add_anno_vfdb_vfs(main_id, vfdb_profile_dir, "core", update_main=False)
            api_anno_vfdb.add_anno_vfdb_pie(main_id, vfdb_profile_dir)
        elif self.option("align_database") == "predict":
            api_anno_vfdb.add_anno_vfdb_vfs(main_id, vfdb_profile_dir, "predict", update_main=False)
        elif self.option("align_database") == "all":
            if os.path.exists(vfdb_profile_dir + "/vfdb_core_VF_profile.xls"):
                api_anno_vfdb.add_anno_vfdb_vfs(main_id, vfdb_profile_dir, "core", update_main=False)
            if os.path.exists(vfdb_profile_dir + "/vfdb_predict_VF_profile.xls"):
                api_anno_vfdb.add_anno_vfdb_vfs(main_id, vfdb_profile_dir, "predict", update_main=False)
            api_anno_vfdb.add_anno_vfdb_pie(main_id, vfdb_profile_dir)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "VFDB毒力因子注释结果目录", 0, "120071"],
            ["vfdb_level_pie.xls", "xls", "VFDB两级分类的丰度统计表", 0, "120080"],
            ["gene_vfdb_anno.xls", "xls", "每条基因的VFDB功能注释表", 0, "120143"],
        ])
        regexps = ([
            [r"vfdb_.*_Gi_profile\.xls", "xls", "各样品基因丰度表", 0, "120144"],
            [r"vfdb_.*_VF_profile\.xls", "xls", "各样品毒力因子丰度表", 0, "120145"],
            ])
        result_dir.add_regexp_rules(regexps)
        super(AnnoVfdbWorkflow, self).end()
