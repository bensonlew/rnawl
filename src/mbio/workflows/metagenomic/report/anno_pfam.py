# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' 20180930

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import unittest
from comfun import ComfunWorkflow


class AnnoPfamWorkflow(ComfunWorkflow):
    """
    宏基因组pfam个性化注释交互分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoPfamWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_pfam workflow")
        self.run_filter(self.set_db)
        super(AnnoPfamWorkflow, self).run()

    def set_db(self):
        """
        将pfam注释结果导入mongo数据库
        :return:
        """
        api_anno_pfam = self.api.api("metagenomic.mg_anno_pfam")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option('final_gene_anno').prop['path']
        anno_file_path = self.output_dir + "/gene_pfam_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        pfam_profile_dir = self.filter.output_dir
        all_files = os.listdir(pfam_profile_dir)
        for i in all_files:
            old = os.path.join(pfam_profile_dir, i)
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
                self.set_error('main_id必须为Object对象或其对应的字符串！', code="12804101")
        #sanger_type, sanger_path = self._sheet.output.split(':')
        #sanger_prefix = Config().get_netdata_config(sanger_type)
        #self.remove_dir =os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir + 'gene_pfam_anno.xls')
        self.logger.info(web_path)
        api_anno_pfam.add_anno_pfam(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_pfam.add_anno_pfam_stat(main_id, pfam_profile_dir)
        api_anno_pfam.add_anno_pfam_detail(main_id, pfam_profile_dir, "pfam",update_main=False)
        api_anno_pfam.add_anno_pfam_detail(main_id, pfam_profile_dir, "type",update_main=False)
        api_anno_pfam.add_anno_pfam_detail(main_id, pfam_profile_dir, "clan",update_main=False)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Pfam结构域注释结果目录",0,"120252"],
            ["gene_pfam_anno.xls", "xls", "每条基因的Pfam功能注释表",0,"120253"],
            ["pfam_type_profile.xls", "xls", "各样品的Pfam Type丰度表",0,"120254"],
            ["pfam_acc_profile.xls", "xls", "各样品的Pfam Pfam丰度表/最低层级表",0,"120255"],
            ["pfam_clan_profile.xls", "xls", "各样品的Pfam CLAN丰度表",0,"120256"],
            ["gene_pfam_anno_stat.xls", "xls", "Pfam 注释基因统计",0,"120257"],
        ])
        super(AnnoPfamWorkflow, self).end()
