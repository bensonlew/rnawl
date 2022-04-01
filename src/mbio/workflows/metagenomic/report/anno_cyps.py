# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'


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


class AnnoCypsWorkflow(ComfunWorkflow):
    """
    宏基因注p450注释交互分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoCypsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_cyps workflow")
        self.run_filter(self.set_db)
        super(AnnoCypsWorkflow, self).run()

    def set_db(self):
        """
        将注释的结果导入mongo数据库中
        :return:
        """
        api_anno_cyps = self.api.api("metagenomic.mg_anno_cyps")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        anno_file_path = self.output_dir + "/gene_cyps_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        cyps_profile_dir = self.filter.output_dir
        all_files = os.listdir(cyps_profile_dir)
        for i in all_files:
            old = os.path.join(cyps_profile_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.logger.info("开始写入mongo数据库")
        main_id = self.option("main_table_id")
        self.logger.info(main_id)
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error("main_id必须为ObjectId对象或其他对应的字符串！", code="12803701")
        #sanger_type, sanger_path = self._sheet.output.split(':')
        #sanger_prefix = Config().get_netdata_config(sanger_type)
        #self.remote_dir = os.path.join(sanger_prefix[sanger_type + "_path"], sanger_path)
        self.remote_dir = self._sheet.output
        web_path = os.path.join(self.remote_dir + "gene_cyps_anno.xls")
        self.logger.info(web_path)
        api_anno_cyps.add_anno_cyps(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_cyps.add_anno_cyps_detail(main_id, cyps_profile_dir)
        api_anno_cyps.add_anno_cyps_abu(main_id, cyps_profile_dir, "homo", update_main=False)
        api_anno_cyps.add_anno_cyps_abu(main_id, cyps_profile_dir, "super", update_main=False)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "P450蛋白功能注释结果目录",0,"120235"],
            ["gene_cyps_anno.xls", "xls", "每条蛋白的P450功能注释表",0,"120236"],
            ["cyps_homo_profile.xls", "xls", "各样品P450 Homologous_family丰度表",0,"120237"],
            ["cyps_super_profile.xls", "xls", "各样品P450 Superfamily丰度表",0,"120238"],
            ["cyps_sid_profile.xls", "xls", "样品P450注释的最低层级表",0,"120239"],
            ["cyps_anno_stat.xls", "xls", "P450注释信息统计表",0,"120240"],
        ])
        super(AnnoCypsWorkflow, self).end()
