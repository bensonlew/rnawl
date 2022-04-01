# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class ExpCorrWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "samples", "type": "string", "default": ""},
            {"name": "dist", "type": "string", "default": "euclidean"},
            {"name": "coefficient", "type": "string", "default": "pearson"},
            {"name": "cluster", "type": "string", "default": "complete"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "transform", "type": "string","default":"UV"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.profile = self.add_tool("metabolome.select_table")
        self.corr_cluster = self.add_tool("metabolome.compare.corr_tree")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Exp_cor workflow")
        if self.option("samples") == "All" and self.option("transform") in ["none","None"]:
            self.run_corr_tre()
        else:
            self.select_profile()
        super(ExpCorrWorkflow, self).run()

    def select_profile(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        options = {
            "origin_table": exp_profile,
            "samples": self.option("samples"),
            "st": "F"
        }
        if self.option("transform") in ["UV", "Par","Ctr"]:
            options["scale"] = True
            options["scale_method"] = self.option("transform")
        self.profile.set_options(options)
        self.profile.on('end', self.run_corr_tree)
        self.profile.run()

    def run_corr_tree(self):
        self.logger.info("start run corr_cluster !")
        if self.option("samples") == "All" and self.option("transform") in ["none","None"]:
            profile = self.option("metab_table").prop["path"]
        else:
            profile = self.profile.option("select_table").prop["path"]
        options = {
            'exp': profile,
            'scd': self.option("dist"),
            'corr_method': self.option("coefficient"),
        }
        if self.option("cluster") != "none":
            options['scm'] = self.option("cluster")
            options['sct'] = "hierarchy"
        self.corr_cluster.set_options(options)
        self.corr_cluster.on('end', self.set_db)
        self.corr_cluster.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.exp_corr")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700601")
        corr_file = self.link_file("corr.xls", "sample_corr.xls")
        p_file = self.link_file("pvalue.xls", "corr_pvalue.xls")
        oldtreefile = os.path.join(self.corr_cluster.output_dir, "corr.cluster_tree.xls")
        if os.path.exists(oldtreefile):
            t_file = self.link_file("corr.cluster_tree.xls", "sample_corr_tree.xls")
            list_file = None
        else:
            t_file = None
            list_file = corr_file
        api_name.add_exp_corr(main_id=main_id, tree_file=t_file,list_file=list_file)
        api_name.add_exp_corr_detail(main_id, corr_file, p_file)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "expcorr",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.corr_cluster.output_dir, oldfile)
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
        result_dir.add_relpath_rules([
            [".", "", "样本相关性结果文件夹", 0, "150009"],
            ["sample_corr_tree.xls", "xls", "样本相关性树文件", 0, "150010"],
            ["sample_corr.xls", "xls", "样本相关性表", 0, "150011"],
            ["corr_pvalue.xls", "xls", "样本相关性P值表", 0, "150012"],
        ])
        super(ExpCorrWorkflow, self).end()
