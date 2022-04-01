# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import types
from mainapp.models.mongo.metabolome import Metabolome
from bson.objectid import ObjectId

class RelationKeggClusterWorkflow(Workflow):
    """
    kegg富集热图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationKeggClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "metab_set_table", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "trans_geneset_main_id", "type": "string"},
            {"name": "trans_kegg_main_id", "type": "string"},
            {"name": "correct", "type": "string"},  # 多重检验方法
            {"name": "species", "type": "string", "default": "all"},  # 代谢富集背景为物种
            {"name": "select", "type": "string", "default": "pvalue"},  # 挑选那列做聚类分析 pvalue, qvaule,in_pop,in_study
            {"name": "pathway_method", "type": "string", "default": "hierarchy"},  # hierarchy, kmeans or none
            {"name": "pathway_dist", "type": "string", "default": "euclidean"},
            {"name": "pathway_ctype", "type": "string", "default": "complete"},
            {"name": "pathway_n_cluster", "type": "int","default": 10},
            {"name": "set_method", "type": "string", "default": "hierarchy"},  # hierarchy, kmeans or none
            {"name": "set_dist", "type": "string", "default": "euclidean"},
            {"name": "set_ctype", "type": "string", "default": "complete"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.enrich_cluster_module = self.add_module("metabolome.relation_enrich_cluster")
        self.task_id = "_".join(self._sheet.id.split("_")[0:2])
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = self.metablome.find_version_from_task(self.task_id)
        self.output_dir = self.enrich_cluster_module.output_dir

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Relation_kegg_cluster workflow")
        self.enrich_cluster_module.on('end', self.set_db)
        self.run_enrich_cluster()
        super(RelationKeggClusterWorkflow, self).run()

    def run_enrich_cluster(self):
        opts = {
            'anno_overview': self.option("anno_overview"),
            'ko_overview': self.option("ko_overview"),
            'metab_set_table': self.option("metab_set_table"),
            "trans_geneset_main_id": self.option("trans_geneset_main_id"),
            "trans_kegg_main_id": self.option("trans_kegg_main_id"),
            'correct': self.option('correct'),
            'bg': "species",
            'species' : self.option("species"),
            "version": self.version_value,
            "select": self.option("select"),
            "pathway_method": self.option("pathway_method"),
            "set_method": self.option("set_method"),      
            "task_id": self.task_id
        }
        if self.option("pathway_method") == "hierarchy":
            opts["pathway_dist"] =  self.option("pathway_dist")
            opts["pathway_ctype"] = self.option("pathway_ctype")
            opts["pathway_n_cluster"] = self.option("pathway_n_cluster")
        elif self.option("pathway_method") == "kmeans":
            opts["pathway_dist"] = self.option("pathway_dist")
            opts["pathway_n_cluster"] = self.option("pathway_n_cluster")
        else:
            pass
        if self.option("set_method") == "hierarchy":
            opts["set_dist"] =  self.option("set_dist")
            opts["set_ctype"] =  self.option("set_ctype")
        elif self.option("pathway_method") == "kmeans":
            opts["set_dist"] =  self.option("set_dist")
        else:
            pass
        self.enrich_cluster_module.set_options(opts)
        self.enrich_cluster_module.run()
        
    def set_db(self):
        api_name = self.api.api("metabolome.relation_kegg_cluster")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        # 导入主表数据
        pathway_tree = None
        set_tree = None
        if os.path.exists(self.output_dir + "/metabset.cluster_tree.xls"):
            set_tree = os.path.join(self.output_dir,"metabset.cluster_tree.xls")
        if os.path.exists(self.output_dir + "/pathway.cluster_tree.xls"):
            pathway_tree = os.path.join(self.output_dir,"pathway.cluster_tree.xls")
        self.logger.info("main_id为{}".format(main_id))
        main_table_id = api_name.add_relation_kegg_heatmap(main_id=main_id, pathway_tree=pathway_tree, set_tree=set_tree)
        
        # 导入详情表数据
        pathway_cluster_table = os.path.join(self.output_dir, "pathway.cluster.txt")
        exp_table = os.path.join(self.output_dir, "exp_table.xls")
        pvalue_table = os.path.join(self.output_dir, "pvalue_table.xls")
        id2db = os.path.join(self.output_dir, "id2db.xls")
        id2term = os.path.join(self.output_dir, "id2term.xls")
        select = self.option("select")
        chose_species = self.option("species").split(";")
        if len(chose_species) == 1:
            species = None
        else:
            species = chose_species[1]
        api_name.add_relation_kegg_heatmap_detail(main_table_id, pathway_table=pathway_cluster_table, exp_table=exp_table, pvalue_table=pvalue_table, id2db=id2db, id2term=id2term, select=select, species=species)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "relationkeggcluster",
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
            [".", "", "代谢物聚类结果文件夹", 0],
            ["col.cluster_tree.xls", "xls", "代谢集和基因集聚类树文件", 0],
            ["pathway.cluster_tree.xls", "xls", "通路聚类树文件", 0],
            ["pathway.kmeans_cluster.xls", "xls", "通路kmeans分类文件", 0],
            ["exp_table.xls", "xls", "矩阵值", 0]
        ]
        regexps = [
            [r"pathway.subcluster_.*\.xls", "xls", "通路各子类结果表", 0],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(RelationKeggClusterWorkflow, self).end()








