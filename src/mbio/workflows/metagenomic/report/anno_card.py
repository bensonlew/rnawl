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

class AnnoCardWorkflow(ComfunWorkflow):
    """
    宏基因组card注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoCardWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_card workflow")
        self.run_filter(self.set_db)
        super(AnnoCardWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno_card = self.api.api("metagenomic.mg_anno_card")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_card_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        card_profile_dir = self.filter.output_dir
        all_files = os.listdir(card_profile_dir)
        for i in all_files:
            old = os.path.join(card_profile_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        self.logger.info(main_id)
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12800301")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir, "gene_card_anno.xls")
        self.logger.info(web_path)
        api_anno_card.add_anno_card(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_card.add_anno_card_aro(main_id, card_profile_dir, update_main=False)
        #api_anno_card.add_anno_card_class(main_id, card_profile_dir)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "CARD抗性基因功能注释结果目录", 0, "120091"],
            ["card_Drug_class_profile.xls", "xls", "各样品Drug class丰度表", 0, "120093"],
            ["card_ARO_gene_number.xls", "xls", "CARD每个ARO比对基因信息表", 0, "120094"],
            ["card_ARO_profile.xls", "xls", "各样品CARD ARO丰度表", 0, "120096"],
            ["gene_card_anno.xls", "xls", "每条基因的CARD功能注释表", 0, "120095"],
            ["card_Antibiotic_class_profile.xls", "xls", "各样品Antibiotic class丰度表", 0, ""],
            ["card_Resistance_Mechanism_profile.xls", "xls", "各样品抗性机制丰度表", 0, ""],
            ["card_ARO_profile_all.xls", "xls", "各样品ARO丰度总表", 0, ""],
        ])
        super(AnnoCardWorkflow, self).end()
