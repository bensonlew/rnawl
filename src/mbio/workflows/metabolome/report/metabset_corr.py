# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os, glob
import re
import shutil
import json
import types


class MetabsetCorrWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_id", "type": "string"},
            #{"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "samples", "type": "string", "default": ""},
            {"name": "dist", "type": "string", "default": "euclidean"},
            {"name": "coefficient", "type": "string", "default": "pearson"},
            {"name": "cluster", "type": "string", "default": "complete"},
            {"name": "cluster_type", "type": "string", "default": "hierarchy"},
            {"name": "n_cluster", "type": "int", "default": ""},
            {"name": "top", "type": "int", "default": ""},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
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
        if self.option("samples") == "All":
            self.run_corr_tre()
        else:
            self.select_profile()
        super(MetabsetCorrWorkflow, self).run()

    def select_profile(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "top": self.option("top"),
            "st": "F",
            "samples": self.option("samples"),
            "select_columns": "metab_id"
        }
        self.profile.set_options(options)
        self.profile.on('end', self.run_corr_tree)
        self.profile.run()

    def run_corr_tree(self):
        self.logger.info("start run corr_cluster !")
        profile = self.profile.option("select_table").prop["path"]
        metab_abu = self.option("metab_table").prop["path"]
        #metab_des = metab_abu.replace("metab_abund.txt", "metab_desc.txt")
        metab_des = self.option("metab_desc").prop["path"]
        if not os.path.exists(metab_des):
            self.set_error('metab_trans路径不存在-%s!', variables=(metab_des), code="14701201")
        options = {
            'exp': profile,
            'sct': self.option("cluster_type"),
            'corr_method': self.option("coefficient"),
            'file_tran': True,
            'metab_trans': metab_des
        }
        if self.option("cluster_type") == "hierarchy":
            options['scm'] = self.option("cluster")
            options['scd'] = self.option("dist")
        if self.option("cluster_type") == "kmeans":
            options['n_cluster'] = self.option("n_cluster")
            options['scd'] = self.option("dist")
        self.corr_cluster.set_options(options)
        self.corr_cluster.on('end', self.set_db)
        self.corr_cluster.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.metabset_corr")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        t_file = None
        matab_tree_id = None
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14701202")
        corr_file = self.link_file("corr.xls", "metab_corr.xls")
        p_file = self.link_file("pvalue.xls", "corr_pvalue.xls")
        if os.path.exists(self.corr_cluster.output_dir + "/corr.cluster_tree.xls"):
            t_file = self.link_file("corr.cluster_tree.xls", "metab_corr_tree.xls")
            matab_tree_id = self.link_file("metab_id.corr_tree.xls", "metab_id.corr_tree.xls")
        list_file = corr_file
        api_name.add_metabset_corr(main_id=main_id, tree_file=matab_tree_id, list_file=list_file)
        api_name.add_metabset_corr_detail(main_id, corr_file, p_file)
        subclusters = glob.glob(self.corr_cluster.output_dir + '/*subcluster*')
        for each in subclusters:
            if os.path.exists(each):
                each = each.split("/")[-1]
                self.link_file(each, each)
        if os.path.exists(self.corr_cluster.output_dir + "/corr.kmeans_cluster.xls"):
            self.link_file("corr.kmeans_cluster.xls", "corr.kmeans_cluster.xls")
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetcorr",
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
        relpath_rules = [
            [".", "", "代谢物相关性结果文件夹", 0, "150031"],
            ["metab_corr_tree.xls", "xls", "代谢物相关性树文件", 0, "150032"],
            ["metab_corr.xls", "xls", "代谢物相关性表", 0, "150033"],
            ["corr_pvalue.xls", "xls", "代谢物相关性P值表", 0, "150034"],
            ["corr.kmeans_cluster.xls", "", "代谢物相关性kmeans分类文件", 0, "150067"],
            ["metab_id.corr_tree.xls", "", "代谢物id相关性树文件", 0, "150076"],
        ]
        regexps = [
            [r"corr.subcluster_.*\.xls", "xls", "代谢物相关性kmeans各子类结果表", 0, "150068"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(MetabsetCorrWorkflow, self).end()
