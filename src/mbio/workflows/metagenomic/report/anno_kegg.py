# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.07.16

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


class AnnoKeggWorkflow(ComfunWorkflow):
    """
    宏基因组kegg注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno_kegg workflow")
        self.run_filter(self.set_db)
        super(AnnoKeggWorkflow, self).run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        self.logger.info("start api")
        api_anno_kegg = self.api.api("metagenomic.mg_anno_kegg")
        geneset_id = self.option("geneset_id")
        specimen = self.option("samples")
        old_anno_file = self.filter.option("final_gene_anno").prop["path"]
        self.logger.info(old_anno_file)
        anno_file_path = self.output_dir + "/gene_kegg_anno.xls"
        self.logger.info("start link gene anno")
        if os.path.exists(anno_file_path):
            os.remove(anno_file_path)
        os.link(old_anno_file, anno_file_path)
        self.logger.info("gene anno ok")
        # KEGG profile_dir 与其他数据库不一致，直接使用计算后的目录,因为有pathway_img目录
        # kegg_profile_dir = self.filter.work_dir + "/RecalAnnoAbu/output"
        kegg_profile_dir = self.filter.recal_abu_tool.option("result_dir").prop["path"]
        self.logger.info(kegg_profile_dir)
        all_files = os.listdir(kegg_profile_dir)
        need_link = ["kegg_enzyme_profile.xls", "kegg_gene_profile.xls", "kegg_KO_profile.xls",
                     "kegg_level1_profile.xls", "kegg_level2_profile.xls", "kegg_level3_profile.xls",
                     "kegg_module_profile.xls", "kegg_pathway_eachmap.xls", "kegg_pathway_profile.xls", "pathway_img"]  # 增加pathway_img
        self.logger.info(all_files)
        for i in all_files:
            if i in need_link:
                old = os.path.join(kegg_profile_dir, i)
                link = os.path.join(self.output_dir, i)
                if os.path.isfile(link):
                    os.remove(link)
                elif os.path.isdir(link):
                    shutil.rmtree(link)
                if os.path.isfile(old):
                    os.link(old, link)
                elif os.path.isdir(old):
                    os.mkdir(link)
                    all_pic = os.listdir(old)
                    for pics in all_pic:
                        old_pic = os.path.join(old, pics)
                        link_pic = os.path.join(link, pics)
                        os.link(old_pic, link_pic)
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        self.remote_dir = self._sheet.output  # .split(":")[-1].lstrip("/")
        web_path = os.path.join(self.remote_dir ,"gene_kegg_anno.xls")
        api_anno_kegg.add_anno_kegg(geneset_id, specimen, web_path, main=False, main_table_id=main_id)
        api_anno_kegg.add_anno_kegg_gene(main_id, kegg_profile_dir, update_main=False)
        api_anno_kegg.add_anno_kegg_orthology(main_id, kegg_profile_dir)
        api_anno_kegg.add_anno_kegg_module(main_id, kegg_profile_dir)
        api_anno_kegg.add_anno_kegg_enzyme(main_id, kegg_profile_dir)
        #self.logger.info(self.option("img_link_dir"))  # 不用这个路径了
        img_link_dir = os.path.join(self.remote_dir, "pathway_img")
        # api_anno_kegg.add_anno_kegg_pathway(main_id, kegg_profile_dir, self.option("img_link_dir"))
        api_anno_kegg.add_anno_kegg_pathway(main_id, kegg_profile_dir, img_link_dir)

        wf_dir = os.path.abspath(os.curdir)
        try:
            os.chdir(self.output_dir)
            os.system('tar -zcf pathway_img.tar.gz pathway_img')
            os.system('rm -r pathway_img')
            os.chdir(wf_dir)
        except:
            os.chdir(wf_dir)
            self.set_error("pathway_img 压缩失败")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "KEGG功能注释结果目录", 0, "120046"],
            ["gene_kegg_anno.xls", "xls", "各样本KEGG基因注释表", 0, "120047"],
            ["kegg_gene_profile.xls", "xls", "各样品KEGG基因丰度表", 0, "120051"],
            ["kegg_KO_profile.xls", "xls", "各样品KO丰度表", 0, "120052"],
            ["kegg_enzyme_profile.xls", "xls", "各样品KEGG酶丰度表", 0, "120050"],
            ["kegg_module_profile.xls", "xls", "各样品KEGG Module丰度表", 0, "120056"],
            ["kegg_pathway_profile.xls", "xls", "各样品Pathway丰度表", 0, "120057"],
            ["kegg_pathway_eachmap.xls", "xls", "各样本Pathway map表", 0, "120049"],
            ["kegg_level1_profile.xls", "xls", "各样品Pathway level1丰度表", 0, "120053"],
            ["kegg_level2_profile.xls", "xls", "各样品Pathway level2丰度表", 0, "120054"],
            ["kegg_level3_profile.xls", "xls", "各样品Pathway level3丰度表", 0, "120055"],
            ["pathway_img.tar.gz", "", "", 1, ""]
        ])
        super(AnnoKeggWorkflow, self).end()
