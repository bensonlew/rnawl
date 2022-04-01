# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

"""距离矩阵层级聚类"""

import datetime
from biocluster.workflow import Workflow
from bson import ObjectId
import re
import os
import json
import shutil
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json

class HclusterWorkflow(CommTableWorkflow):
    """
    报告中调用距离矩阵计算样本层级聚类数使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HclusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "profile_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "distance_method", "type": "string", "default": "bray_curtis"},
            {"name": "hcluster_method", "type": "string", "default": "average"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "group_id", "type": "string", "default": ""},
            # {"name": "abu_id", "type": "string","default": ""},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        self.hcluster = self.add_tool("meta.beta_diversity.hcluster")
        #self.abundance = self.add_tool("meta.create_abund_table")
        self.sam = self.add_tool("meta.otu.sort_samples_mg")

    def run(self):
        # self.IMPORT_REPORT_DATA = True
        # self.IMPORT_REPORT_DATA_AFTER_END = False
        # self.abundance.on('end', self.sort_sample)
        # self.sam.on('end', self.run_dist)
        # self.dist.on('end', self.run_hcluster)
        # self.hcluster.on('end', self.set_db)
        if self.option("group_table").is_set:
            self.sam.on('end', self.run_dist)
            if not self.option("profile_table").is_set:
                #self.abundance.on('end', self.sort_sample)
                self.run_abundance(self.sort_sample)
        else:
            if not self.option("profile_table").is_set:
                #self.abundance.on('end', self.run_dist)
                self.run_abundance(self.run_dist)
        # super(HclusterWorkflow, self).run()
        # self.abundance.on('end', self.sort_sample)
        # self.sam.on('end', self.run_dist)
        self.dist.on('end', self.run_hcluster)
        self.hcluster.on('end', self.set_db)
        if self.option("profile_table").is_set:
            if self.option("group_table").is_set:
                self.sort_sample()
            else:
                self.run_dist()
        else:
            #self.run_abundance()
            self.abundance.run()
        super(HclusterWorkflow, self).run()

    def run_dist(self):
        if self.option("group_table").is_set:
            otutable = self.sam.option("out_otu_table")
        else:
            if self.option("profile_table").is_set:
                otutable = self.option("profile_table")
            else:
                otutable = self.abundance.option('out_table')
        options = {
            'method': self.option('distance_method'),
            'otutable': otutable
        }
        self.dist.set_options(options)
        self.dist.run()

    def run_hcluster(self):
        options = {
            'linkage': self.option('hcluster_method'),
            'dis_matrix': self.dist.option('dis_matrix')
        }
        self.hcluster.set_options(options)
        self.hcluster.run()

    def sort_sample(self):
        if self.option("profile_table").is_set:
            otutable = self.option("profile_table")
        else:
            otutable = self.abundance.option('out_table')
        options = {
            'in_otu_table': otutable,
            'group_table': self.option("group_table")
        }
        self.sam.set_options(options)
        self.sam.run()

    def set_db(self):
        params_json = json.loads(self.option('params'))
        matrix_path = self.dist.output_dir + '/' + os.listdir(self.dist.output_dir)[0]
        final_matrix_path = os.path.join(self.output_dir, os.listdir(self.dist.output_dir)[0])
        shutil.copy2(matrix_path, final_matrix_path)
        if not os.path.isfile(matrix_path):
            self.set_error("找不到报告文件:%s", variables=(matrix_path), code="12801601")
        self.api_dist = self.api.api('metagenomic.distance_metagenomic')
        dist_id = self.api_dist.add_dist_table(matrix_path, main=True, name=None,params=params_json, task_id=self.option('task_id'))
        newick_fath = self.hcluster.output_dir + "/hcluster.tre"
        final_newick_path = os.path.join(self.output_dir, "hcluster.tre")
        shutil.copy2(newick_fath, final_newick_path)
        if not os.path.isfile(newick_fath):
            self.set_error("找不到报告文件:%s", variables=(newick_fath), code="12801602")
        self.api_tree = self.api.api('metagenomic.hcluster_tree')
        self.api_tree.add_hcluster_tree(newick_fath, main=False, tree_id=self.option('main_id'),
                                        update_dist_id=dist_id)
        # self.api_abu = self.api.api('metagenomic.abund_table_path')
        # if self.option("abu_id"):
        #     self.api_abu.add_abund_table_path(self.abundance.output_dir + '/new_abund_table.xls',self.option("abu_id"),main=False)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "hcluster_tree")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = json.loads(self.option('params'))["submit_location"]
            # self.logger.info(submit_loc)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir_hucluster = self.add_upload_dir(self.output_dir)
        result_dir_hucluster.add_relpath_rules([
            [".", "", "样本层次聚类分析结果目录", 0, "120104"],
            ["./hcluster.tre", "tre", "样本层次聚类树结果表", 0, "120105"],  # modified by hongdongxuan 20170321
            ["./hcluster.pdf", "pdf", "层级聚类图"],
            ["./distance_Heatmap.pdf", "pdf", "距离Heatmap图"],
        ])
        result_dir_hucluster.add_regexp_rules([
            [r'%s.*\.xls' % self.option('distance_method'), 'xls', '样本距离矩阵文件', 0, "120107"]
        ])
        super(HclusterWorkflow, self).end()
