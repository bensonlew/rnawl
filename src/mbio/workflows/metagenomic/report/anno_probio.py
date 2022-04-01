# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified

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

class AnnoProbioWorkflow(ComfunWorkflow):
    """
    宏基因组ardb注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoProbioWorkflow, self).__init__(wsheet_object)
        options = [
            #{"name": "name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_probio workflow")
        self.run_filter(self.set_db)
        super(AnnoProbioWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno = self.api.api("metagenomic.probiotics")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        profile_dir = self.filter.output_dir
        all_files = os.listdir(profile_dir)
        for i in all_files:
            old = os.path.join(profile_dir, i)
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12804001")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir, "gene_probio_anno.xls")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        new_anno_file = self.output_dir + "/gene_probio_anno.xls"
        anno_detail_path = self.output_dir + "/probio_anno.xls"
        abun_file_path = self.output_dir + "/probio_abun.xls"
        if os.path.exists(new_anno_file):
            os.remove(new_anno_file)
        os.link(old_anno_file, new_anno_file)
        api_anno.add_anno_probio(geneset_id, specimen, web_path, main=False,
                                    main_table_id=main_id)
        api_anno.add_anno_probio_detail(main_id, anno_detail_path, update_main=False)
        api_anno.add_anno_probio_abun(main_id, abun_file_path)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "益生菌功能注释结果目录", 0, ""],
            ["gene_probio_anno.xls", "xls", "每条基因的益生菌注释表", 0, ""],
            ["probio_anno.xls", "xls", "益生菌注释详情表", 0, ""],
            ["probio_abun.xls", "xls", "益生菌丰度表", 0, ""],
        ])
        super(AnnoProbioWorkflow, self).end()

