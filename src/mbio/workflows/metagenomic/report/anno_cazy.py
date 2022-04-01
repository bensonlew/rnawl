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


class AnnoCazyWorkflow(ComfunWorkflow):
    """
    宏基因组cazy注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoCazyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_cazy workflow")
        self.run_filter(self.set_db)
        super(AnnoCazyWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno_cazy = self.api.api("metagenomic.mg_anno_cazy")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_cazy_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        cazy_profile_dir = self.filter.output_dir
        all_files = os.listdir(cazy_profile_dir)
        for i in all_files:
            old = os.path.join(cazy_profile_dir, i)
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12804901")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir, "gene_cazy_anno.xls")
        api_anno_cazy.add_anno_cazy(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_cazy.add_anno_cazy_family(main_id, cazy_profile_dir, update_main=False)
        api_anno_cazy.add_anno_cazy_class(main_id, cazy_profile_dir)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "CAZy碳水化合物活性酶注释结果目录", 0, "120064"],
            ["gene_cazy_anno.xls", "xls", "每条基因的CAZy功能注释表", 0, "120066"],
            ["cazy_family_profile.xls", "xls", "各样品CAZy Family丰度表", 0, "120067"],
            ["cazy_class_profile.xls", "xls", "各样品CAZy Class丰度表", 0, "120065"],
            ["gene_cazy_class_stat.xls", "xls", "CAZy Class基因信息统计表", 0, "120069"],
            ["gene_cazy_family_stat.xls", "xls", "CAZy Family基因信息统计表", 0, "120070"]
        ])
        super(AnnoCazyWorkflow, self).end()
