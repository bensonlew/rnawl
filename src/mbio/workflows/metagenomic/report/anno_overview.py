# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2017.11.2


from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from biocluster.file import exists
from mbio.api.database.metagenomic.mg_anno_overview import MgAnnoOverview

class AnnoOverviewWorkflow(Workflow):
    """
    宏基因组注释总览交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoOverviewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "select", "type": "string"},
            # {"name": "select_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_length_table", "type": "infile", "format": "annotation.mg_anno_dir"},
            #{"name": "task_id", "type": "string"},
            {"name": "type", "type": "int", "default": 1,"choose":[1,2,3,4]},
            # 1为总览部分重建基因集，2为kegg部分重建基因集，3，基因集管理页面创建，4，基因集管理页面合并
            {"name": "task_id", "type": "string"},
            {"name": "geneset_list", "type": "string"},  #合并的geneset_list
            {"name": "database_list", "type": "string", "default": "cog,genus,phylum,kegg"},  #type为1,3,4时使用
            {"name": "samples", "type": "string"},  #type为2时使用
            {"name": "update_info", "type": "string"},
            {"name": "gene_kegg_anno", "type": "infile", "format": "sequence.profile_table"},  # type为2时使用
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_relative_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "select_type", "type": "string"},
            {"name":"sum_anno_table","type":"infile", "format": "annotation.mg_anno_dir"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gene_select = self.add_tool("meta.annotation.overview_select")
        self.reset = self.add_tool("meta.annotation.reset_geneset")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        gene_dir = self.option("sum_anno_table").prop["path"]
        if self.option("type") ==2:
            self.overview_file = self.option("gene_kegg_anno")
        else:
            if exists(gene_dir + "/gene_overview_anno.xls"):
                self.overview_file = os.path.join(gene_dir,"gene_overview_anno.xls")
            elif exists(gene_dir + "/anno_overview.xls"):
                self.overview_file = os.path.join(gene_dir, "anno_overview.xls")
            else: ###加逻辑判断，弥补找不到overview表的情况
                overview = MgAnnoOverview(self)
                task_id = self.option("task_id")
                self.overview_file = overview.from_mongo_to_file(task_id)
        if self.option("type") != 4:
            self.run_select_gene()
        else:
            self.reset_geneset()
        super(AnnoOverviewWorkflow, self).run()

    def run_select_gene(self):
        # gene_length_table = self.option("gene_length_table").prop["path"]
        # gene_dir = "/".join(gene_length_table.split("/")[0:len(gene_length_table.split("/")) - 2])
        gene_dir = self.option("sum_anno_table").prop["path"]
        type = self.option("type")
        if type == 1:
            options = {
                "select": self.option("select"),
                "type": self.option("type"),
                "database_list":self.option("database_list"),
                "overview_file": self.overview_file
            }
        elif type == 2:
            options = {
                "select": self.option("select"),
                "type": self.option("type"),
                "gene_kegg_anno": self.option("gene_kegg_anno")
            }
        elif type == 3:
            options = {
                "select": self.option("select"),
                "type": self.option("type"),
                "select_type": self.option("select_type"),
                "task_id": self.option("task_id"),
                "overview_file": self.overview_file
            }
        elif type ==4:
            options = {
                "type": self.option("type"),
                "task_id": self.option("task_id"),
                "overview_file": self.overview_file
            }
        self.logger.info(gene_dir)
        self.logger.info(os.path.join(gene_dir, "anno_overview.xls"))
        self.gene_select.set_options(options)
        self.gene_select.on('end', self.reset_geneset)
        self.gene_select.run()

    def reset_geneset(self):
        gene_length_table = os.path.join(self.option("gene_length_table").prop["path"], "gene_profile/gene.uniGeneset.fa.length.txt")
        gene_profile = self.option("gene_profile").prop["path"]
        gene_relative_profile = self.option("gene_relative_profile").prop["path"]
        options = {
            "overviewfile": self.overview_file,
            #"select_genes": select_genes,
            "gene_length_table": gene_length_table,
            "gene_profile": gene_profile,
            "gene_relative_profile": gene_relative_profile,
            "geneset_list": self.option("geneset_list"),
            "mytype": self.option("type"),
        }
        if self.option("type") != 4:
            options["select_genes"] = self.gene_select.option("gene_list").prop["path"]
        self.reset.set_options(options)
        self.reset.on('end', self.set_db)
        self.reset.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_geneset = self.api.api("metagenomic.geneset")
        #specimen = self.option("samples")
        #params = self.option("params")
        #name = self.option("geneset")
        self.logger.info(self.reset.output_dir)
        length_path = self.reset.output_dir + "/length_distribute"
        top100_profile = self.reset.output_dir + "/gene_profile/top100_reads_number.xls"
        top100_pro_relative = self.reset.output_dir + "/gene_profile/top100_reads_number_relative.xls"
        if not os.path.exists(self.output_dir + "/gene_profile"):
            os.mkdir(self.output_dir + "/gene_profile")
        if not os.path.exists(self.output_dir + "/uniGeneset"):
            os.mkdir(self.output_dir + "/uniGeneset")
        gene_profile_all = os.listdir(self.reset.output_dir + "/gene_profile")
        for i in gene_profile_all:
            old = os.path.join(self.reset.output_dir + "/gene_profile", i)
            link = os.path.join(self.output_dir + "/gene_profile", i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        if os.path.exists(self.output_dir + "/uniGeneset/geneCatalog_stat.xls"):
            os.remove(self.output_dir + "/uniGeneset/geneCatalog_stat.xls")
        os.link(self.reset.output_dir + "/uniGeneset/geneCatalog_stat.xls",
                self.output_dir + "/uniGeneset/geneCatalog_stat.xls")
        type = self.option("type")
        anno_file = "anno_overview.xls"
        anno_file_old = os.path.join(self.reset.output_dir, anno_file)
        self.logger.info(anno_file)
        anno_file_new = os.path.join(self.output_dir, anno_file)
        if os.path.exists(anno_file_new):
            os.remove(anno_file_new)
        os.link(anno_file_old, anno_file_new)
        #overview_profile_dir = self.recal_anno_abu.output_dir
        self.logger.info("正在写入mongo数据库")
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # web_path = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # web_path = self._sheet.output.split(":")[-1].lstrip("/")
        web_path = self._sheet.output
        main_id = self.option("main_table_id")
        task_id = self.option("task_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12800801")
        try:
            api_geneset.add_geneset(self.output_dir, web_path, type=2, main_id=main_id, task_id=task_id)
            api_geneset.add_geneset_bar(main_id, length_path)
            api_geneset.add_geneset_readsn(main_id, top100_profile)
            api_geneset.add_geneset_readsr(main_id, top100_pro_relative)
        except Exception as e:
            self.set_error("重建基因集导表失败——%s", variables=(e), code="12800802")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "创建基因集目录", 0, "120146"],
            ["/uniGeneset/geneCatalog_stat.xls", "xls", "去冗余前后基因数目和长度统计表", 0, "120022"],
            ["/gene_profile/reads_number_relative.xls", "xls", "基因在各个样品中的相对丰度表", 0, "120028"],
            ["/gene_profile/reads_number.xls", "xls", "基因在各个样品中的丰度表", 0, "120029"],
            ["/gene_profile/top100_reads_number.xls", "xls", "丰度前100的基因丰度表", 0, "120031"],
            ["/gene_profile/top100_reads_number_relative.xls", "xls", "丰度前100的基因相对丰度表", 0, "120032"],
            ["/gene_profile/reads_profile.tar.gz", "gz", "基因reads数相对丰度与基因reads数丰度的压缩文件", 0, "120148"],
            ["/gene_profile/gene.uniGeneset.fa.length.txt", "txt", "基因长度表", 0, "120034"],
        ])
        regexps = ([
            [r"anno_.*\.xls", "xls", "筛选的基因注释结果",0,"120251"]
        ])
        result_dir.add_regexp_rules(regexps)
        super(AnnoOverviewWorkflow, self).end()
