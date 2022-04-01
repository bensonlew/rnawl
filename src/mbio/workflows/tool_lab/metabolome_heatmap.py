# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last_modifiy = 2021.04.16

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from mbio.packages.tool_lab.common_function import rename_name,rename_name_back


class MetabolomeHeatmapWorkflow(Workflow):
    """
    小工具 代谢聚类热图
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabolomeHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "sam_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "group_method", "type": "string", "default": "none"},
            {"name": "top_meta", "type": "int", "default": 50},
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "n_cluster", "type": "int", "default": 10},
            {"name": "sam_dist", "type": "string", "default": "euclidean"},
            {"name": "sam_cluster", "type": "string", "default": "complete"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "scale", "type": "string", "default": "scale"},  ## 是否标准化聚类数据，scale, unscale

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.profile = self.add_tool("tool_lab.select_table")
        self.metab_cluster = self.add_tool("tool_lab.metab_cluster")
        if self.option("metab_cluster_method") == "none":
            self.metab_cluster_method = ""
        else:
            self.metab_cluster_method = self.option("metab_cluster_method")
        if self.option("sam_cluster_method") == "none":
            self.sam_cluster_method = ""
        else:
            self.sam_cluster_method = self.option("sam_cluster_method")
        self.table_file = self.work_dir + "/table_new.txt"
        with open(self.option("table").prop["path"],"r") as f, open(self.table_file,"w") as t:
            data = f.readlines()
            if data[0].split("\t")[1] == "metab_id":
                for i in data:
                    t.write(re.sub(r",|'|/|:","_",i.strip().split("\t")[0].replace("(","").replace(")","").replace("'","").replace('"',"")) + "\t" + "\t".join(i.strip().split("\t")[2:])+ "\n")
            else:
                for i in data:
                    t.write(re.sub(r",|'|/|:","_",i.strip().split("\t")[0].replace("(","").replace(")","").replace("'","").replace('"',"")) + "\t" + "\t".join(i.strip().split("\t")[1:]) + "\n")


    def run(self):
        self.logger.info("start run metabset_cluster workflow")
        self.profile.on('end', self.run_cluster)
        self.metab_cluster.on('end', self.set_db)
        self.select_profile()
        super(MetabolomeHeatmapWorkflow, self).run()


    def select_profile(self):
        self.logger.info("start profile!")
        self.data_table = self.work_dir + "/table_new2.txt"
        self.name_dict = rename_name(self.table_file, self.data_table)
        if self.option("group_method") != "no":
            options = {
                "origin_table": self.data_table,
                "group_method": self.option("group_method"),
                "group":self.option("group_table"),
            }
        else:
            options = {
                "origin_table": self.data_table,
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
        with open(self.work_dir + "/desc.txt","w") as t, open(self.profile.option("select_table").prop["path"],"r") as f:
            data = f.readlines()
            t.write("Metabolite\tmetab_id\n")
            num = 1
            for i in data[1:]:
                t.write(i.split("\t")[0] + "\t" + i.split("\t")[0] + "\n")
        options = {
            'exp': profile,
            'sct': self.sam_cluster_method,
            'mct': self.metab_cluster_method,
            'metab_trans': self.work_dir + "/desc.txt"
        }
        if  self.metab_cluster_method == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['mcd'] = self.option("metab_dist")
            options['n_cluster'] = self.option("n_cluster")
        if  self.sam_cluster_method == "hierarchy":
            options['scm'] = self.option("sam_cluster")
            options['scd'] = self.option("sam_dist")
        if  self.metab_cluster_method == "kmeans":
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
        api_name = self.api.api("tool_lab.metabset_cluster")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        self.output_dir = self.metab_cluster.output_dir
        sam_tree = None
        list_file = None
        matab_tree_id = None
        sample_class = None
        metab_class = None
        for file in os.listdir(self.output_dir):
            rename_name_back(self.output_dir + "/" + file,self.name_dict)
        if os.path.exists(self.output_dir + "/sample.cluster_tree.xls"):
            sam_tree = os.path.join(self.output_dir,"sample.cluster_tree.xls")
        if os.path.exists(self.output_dir + "/metab.cluster_tree.xls"):
            matab_tree_id = os.path.join(self.output_dir,"metab_id.cluster_tree.xls")
        if os.path.exists(self.output_dir + "/cluster_exp.xls"):
            list_file = os.path.join(self.output_dir,"cluster_exp.xls")
        if self.option("group_method") == "none":
            sample_class = self.option("group_table").prop["path"]
        if os.path.exists(self.output_dir + "/metab." + self.metab_cluster_method + "_cluster.xls"):
            metab_class = os.path.join(self.output_dir + "/metab." + self.metab_cluster_method + "_cluster.xls")
        expression_file = os.path.join(self.output_dir,"cluster_exp.xls")
        main_table_id = api_name.add_metabolome_heatmap_detail(sam_tree=sam_tree, matab_tree=matab_tree_id, main_id=main_id,list_file=list_file,metab_class=metab_class,sample_class=sample_class)
        #api_name.add_metabolome_heatmap_main(main_table_id, expression_file)
        self.end()

    def end(self):
        super(MetabolomeHeatmapWorkflow, self).end()
