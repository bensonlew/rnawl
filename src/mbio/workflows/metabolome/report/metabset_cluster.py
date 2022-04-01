# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0710

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class MetabsetClusterWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_id", "type": "string"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "sam_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "group_method", "type": "string", "default": "average"},
            {"name": "top_meta", "type": "int", "default": 50},
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "n_cluster", "type": "int", "default": 10},
            {"name": "sam_dist", "type": "string", "default": "euclidean"},
            {"name": "sam_cluster", "type": "string", "default": "complete"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "scale", "type": "string", "default": "scale"},  ## 是否标准化聚类数据，scale, unscale
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.profile = self.add_tool("metabolome.select_table")
        self.metab_cluster = self.add_tool("metabolome.metabset.metab_cluster")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run metabset_cluster workflow")
        self.profile.on('end', self.run_cluster)
        self.metab_cluster.on('end', self.set_db)
        self.select_profile()
        super(MetabsetClusterWorkflow, self).run()

    def select_profile(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        if self.option("group_method") == "average":
            group_method = 2
        elif self.option("group_method") == "median":
            group_method = 3
        elif self.option("group_method") == "sum":
            group_method = 1
        elif self.option("group_method") == "no":
            group_method = 0
        else:
            self.set_error("group_method方法不对-%s", variables=(self.option("group_method")), code="14700701")
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "top": self.option("top_meta"),
            "st": "F",
            "group_method": group_method,
            "group":self.option("group_table"),
            "select_columns": "metab_id",
            "merge": "nomerge"
        }
        if self.option("scale") == "scale":
            options["scale"] = True
        else:
            options["scale"] = False
        self.profile.set_options(options)
        self.profile.run()

    def run_cluster(self):
        self.logger.info("start run metabset_cluster !")
        profile = self.profile.option("select_table")
        metab_abu = self.option("metab_table").prop["path"]
        #metab_des = metab_abu.replace("metab_abund.txt","metab_desc.txt")
        metab_des = self.option("metab_desc").prop["path"]
        if not os.path.exists(metab_des):
            self.set_error('metab_trans路径不存在-%s!', variables=(metab_des), code="14700702")
        options = {
            'exp': profile,
            'sct': self.option("sam_cluster_method"),
            'mct': self.option("metab_cluster_method"),
            'metab_trans': metab_des
        }
        if  self.option("metab_cluster_method") == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['mcd'] = self.option("metab_dist")
            options['n_cluster'] = self.option("n_cluster")
        if  self.option("sam_cluster_method") == "hierarchy":
            options['scm'] = self.option("sam_cluster")
            options['scd'] = self.option("sam_dist")
        if  self.option("metab_cluster_method") == "kmeans":
            options['n_cluster'] = self.option("n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("scale") == "scale":
            options["before_scale"] = self.profile.option("select_origin_abu")
        self.logger.info(options)
        self.metab_cluster.set_options(options)
        self.metab_cluster.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.metabset_cluster")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700703")
        self.output_dir = self.metab_cluster.output_dir
        sam_tree = None
        matab_tree = None
        list_file = None
        matab_tree_id = None
        if os.path.exists(self.output_dir + "/sample.cluster_tree.xls"):
            sam_tree = os.path.join(self.output_dir,"sample.cluster_tree.xls")
        if os.path.exists(self.output_dir + "/metab.cluster_tree.xls"):
            matab_tree = os.path.join(self.output_dir,"metab.cluster_tree.xls")
            matab_tree_id = os.path.join(self.output_dir,"metab_id.cluster_tree.xls")
        if not os.path.exists(self.output_dir + "/sample.cluster_tree.xls"):
            list_file = os.path.join(self.output_dir,"cluster_exp.xls")
        if self.option("scale") == "scale":
            expression_file = os.path.join(self.output_dir,"cluster_scale_exp.xls")
        else:
            expression_file = os.path.join(self.output_dir,"cluster_exp.xls")
        api_name.add_metabset_cluster(sam_tree=sam_tree, matab_tree=matab_tree_id, main_id=main_id,list_file=list_file)
        api_name.add_metabset_cluster_detail(main_id, expression_file)
        self.output_dir = self.metab_cluster.output_dir
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetcluster",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules =[
            [".", "", "代谢物聚类结果文件夹", 0, "150027"],
            ["sample.cluster_tree.xls", "xls", "样本聚类树文件", 0, "150029"],
            ["metab.cluster_tree.xls", "xls", "代谢物聚类树文件", 0, "150028"],
            ["cluster_exp.xls", "xls", "代谢物表达量", 0, "150030"],
            ["metab.kmeans_cluster.xls", "xls", "代谢物kmeans分类文件", 0, "150065"],
            ["metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
            ["cluster_scale_exp.xls", "xls", "标准化后代谢物表达量", 0, "150077"],
        ]
        regexps = [
            [r"metab.subcluster_.*\.xls", "xls", "代谢物kmeans各子类结果表", 0, "150066"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(MetabsetClusterWorkflow, self).end()
