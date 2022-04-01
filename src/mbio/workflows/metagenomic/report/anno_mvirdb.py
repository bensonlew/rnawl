# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifiy = modified 20181115


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
from mainapp.models.mongo.metagenomic import Metagenomic
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class AnnoMvirdbWorkflow(ComfunWorkflow):
    """
    宏基因组TCDB注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoMvirdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.metagenomic = Metagenomic()
        self.metagenomic._config = Config()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno workflow")
        self.run_filter(self.set_db)
        super(AnnoMvirdbWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_anno = self.api.api("metagenomic.mvirdb")
        #geneset_id = self.option("geneset_id")
        #specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        self.logger.info(old_anno_file)
        anno_file_path = self.output_dir + "/function_select_anno.xls"
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
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
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12805001")
        self.logger.info(main_id)
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.remote_dir = self._sheet.output.split(":")[-1].lstrip("/")
        self.remote_dir = self._sheet.output
        #web_path = os.path.join(self.remote_dir, "gene_cog_anno.xls")
        api_anno.add_mvirdb_detail(anno_file_path ,main_id)
        for i in all_files:
            i_basename = os.path.basename(i)
            if '_abund_out.xls' in  i_basename:
                level = i_basename.split('_',1)[0].lower()
                api_anno.add_mvirdb_abund(self.output_dir + '/' + i,main_id,level)
        anno_file = os.path.join(self._sheet.output , 'function_select_anno.xls')
        self.metagenomic.common_update_one('anno_mvir', main_id, {'anno_file': anno_file})
        if self.option("save_pdf"):
        # self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        # if self.pdf_status:
            name = get_name(self.option("main_table_id"), "anno_mvir")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "anno_mvir")
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
            [".", "", "MvirDB功能注释结果目录", 0, "120058"],
            ["Factor_abund_out.xls", "xls", "各样品Factor丰度表", 0, "120060"],
            ["Type_abund_out.xls", "xls", "各样品Type丰度表",0,"120247"],
            ["Factor_gene_stat.xls", "xls", "各样品Factor基因列表",0,"120248"],
            ["Type_gene_stat.xls", "xls", "各样品Type基因列表",0,"120249"],
            ["function_select_anno.xls", "xls", "每条基因的MvirDB功能注释表",0,"120250"],
            ["Mvirdb_type_bar.pdf", "pdf", "Type丰度图"]
        ])

        super(AnnoMvirdbWorkflow, self).end()
