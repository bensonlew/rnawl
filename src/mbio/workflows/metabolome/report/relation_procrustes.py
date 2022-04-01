# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modifiy = 2021.11.24

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os, glob
import json
import types


class RelationProcrustesWorkflow(Workflow):
    """
    Procrustes description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationProcrustesWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},
            #{"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "log10", "type": "bool", "default": False},
            {"name": "metab_set_table", "type": "string"},
            {"name": "group_detail", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "trans_exp_main_id", "type": "string"},  # 转录表达量表id
            {"name": "trans_geneset_main_id", "type": "string"},  # 转录基因集表id
            {"name": "method", "type": "string", "default": "pca"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.procrustes = self.add_module("metabolome.procrustes")
        self.metab_tool = self.add_tool("metabolome.select_table")
        self.trans_tool = self.add_tool("metabolome.relation.trans_select_table")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run RelationProcrustes workflow")
        self.select_tools = [self.metab_tool, self.trans_tool]
        self.on_rely(self.select_tools, self.run_procrustes)
        self.procrustes.on('end', self.set_db)
        self.select_metab_table()
        self.select_trans_table()
        super(RelationProcrustesWorkflow, self).run()

    def select_metab_table(self):
        exp_profile = self.option("metab_table").prop["path"]
        options = {
            "origin_table": exp_profile,
            "st": "F",
            "group": self.option("group_detail"),
            "select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": 0,
            "log10": self.option("log10")
        }
        if self.option("metab_set_table") != "all":
            options["select_genes"] = self.option("metab_set_table")
        else:
            pass
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_trans_table(self):
        options = {
            "group": self.option("group_detail"),
            "trans_exp_main_id": self.option("trans_exp_main_id"),
            "task_id" : "_".join(self._sheet.id.split("_")[0:2])
        }
        if self.option("trans_geneset_main_id") != "all":
            options["trans_geneset_main_id"] = self.option("trans_geneset_main_id")
        else:
            pass
        self.trans_tool.set_options(options)
        self.trans_tool.run()

    def run_procrustes(self):
        self.logger.info("start run procrustes !")
        metab_table = self.metab_tool.option("select_table")
        trans_table = self.trans_tool.option("select_table")
        options = dict(
            metab_table = metab_table,
            asso_table = trans_table,
            method = self.option('method')
        )
        self.procrustes.set_options(options)
        self.procrustes.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        self.dump_tool = self.api.api("metabolome.relation_procrustes")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('ERROR: wrong format of main_id, ObjectId or string expected!')
        #if main_id is None:
        #    main_id = self.dump_tool.add_procrustes()  # 测试用导入主表
        summary = os.path.join(self.procrustes.output_dir, "procrustes_results.txt")
        self.dump_tool.add_procrustes_detail(summary, main_id)
        reg1 = os.path.join(self.procrustes.output_dir,"*transformed_reference.txt")
        reg2 = os.path.join(self.procrustes.output_dir,"*transformed_q*.txt")
        oldfile1 = glob.glob(reg1)
        oldfile2 = glob.glob(reg2)
        transform_ref = self.link_file(oldfile1[0], "metab.transformed_reference.txt",0)
        transform_query = self.link_file(oldfile2[0], "gene.transformed_query.txt",0)
        self.dump_tool.add_procrustes_graph(transform_ref, transform_query, main_id)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetprocrustes",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def link_file(self, oldfile, newfile, mode):
        if mode==1:
            oldfile = os.path.join(self.procrustes.output_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)
        return newfile

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析普氏分析结果文件夹", 0, "150079"],
            ["metab.transformed_reference.txt", "txt", "代谢物拓扑转换后结果", 0, "150080"],
            ["gene.transformed_query.txt", "txt", "基因拓扑转换后结果", 0, "150081"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        super(RelationProcrustesWorkflow, self).end()
