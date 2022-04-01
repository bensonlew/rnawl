# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modifiy = modified 2021.11.26

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import os, glob
from mainapp.models.mongo.metabolome import Metabolome
import types
from bson.objectid import ObjectId
from biocluster.config import Config


class RelationCorrWorkflow(Workflow):
    """
    相关性分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},  # 代谢集id
            {"name": "metab_anno", "type": "infile", "format": "sequence.profile_table"}, # 代谢物注释总表
            {"name": 'metab_hmdb_anno', "type": "infile", "format": "sequence.profile_table"},
            {"name": "trans_exp_main_id", "type": "string"},  # 转录表达量表id
            {"name": "trans_geneset_main_id", "type": "string"},  # 转录基因集表id
            {"name": "group_detail", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            #参数选项
            {"name": "log10", "type": "bool", "default": False},
            {"name": "group_method", "type": "int", "default": 0},  # 0、2、3
            {"name": "metab_top", "type": "int", "default": 50},
            {"name": "trans_top", "type": "int", "default": 50},
            {"name": "coefficient", "type": "string", "default": "pearson"},  # 相关性算法
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "metab_n_cluster", "type": "int", "default": 5},
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},  # 层级聚类算法
            {"name": "trans_dist", "type": "string", "default": "euclidean"},
            {"name": "trans_cluster", "type": "string", "default": "complete"},
            {"name": "trans_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "trans_n_cluster", "type": "int", "default": 5},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.metab_tool = self.add_tool("metabolome.select_table")
        self.trans_tool = self.add_tool("metabolome.relation.trans_select_table")
        self.asso_corr_tool = self.add_tool("metabolome.metabset.asso_corr")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Asso_cor workflow")
        self.select_tools = [self.metab_tool, self.trans_tool]
        self.on_rely(self.select_tools, self.run_asso_corr)
        self.asso_corr_tool.on('end', self.set_db)
        self.select_metab()
        self.select_trans()
        super(RelationCorrWorkflow, self).run()

    def select_metab(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "st": "F",
            "group": self.option("group_detail"),
            "select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": self.option("group_method"),
            'top': self.option("metab_top"),
            "log10": self.option("log10")
        }
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_trans(self):
        options = {
            "group": self.option("group_detail"),
            "group_method": self.option("group_method"),
            "trans_exp_main_id": self.option("trans_exp_main_id"),
            "trans_geneset_main_id": self.option("trans_geneset_main_id"),
            "task_id" : "_".join(self._sheet.id.split("_")[0:2]),
            'top': self.option("metab_top")
        }
        self.trans_tool.set_options(options)
        self.trans_tool.run()

    def run_asso_corr(self):  # 相关性计算
        self.logger.info("start run asso_corr !")
        metab_table = self.metab_tool.option("select_table")
        trans_table = self.trans_tool.option("select_table")
        metab_des = self.option("metab_desc").prop["path"]
        options = {
            'metab_table': metab_table,
            'asso_table': trans_table,
            'mct': self.option("metab_cluster_method"),
            'sct': self.option("trans_cluster_method"),
            'coefficient': self.option("coefficient"),
            'metab_trans': metab_des
        }
        # 层级聚类
        if self.option("metab_cluster_method") == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['metab_n_cluster'] = self.option("metab_n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("trans_cluster_method") == "hierarchy":
            options['scm'] = self.option("trans_cluster")
            options['asso_n_cluster'] = self.option("trans_n_cluster")
            options['scd'] = self.option("trans_dist")
        # kmeans
        if self.option("metab_cluster_method") == "kmeans":
            options['metab_n_cluster'] = self.option("metab_n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("trans_cluster_method") == "kmeans":
            options['asso_n_cluster'] = self.option("trans_n_cluster")
            options['scd'] = self.option("trans_dist")
        self.asso_corr_tool.set_options(options)
        self.asso_corr_tool.run()
        
    def set_db(self):
        api_name = self.api.api("metabolome.relation_corr")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        assoc_tree = None
        metab_tree_id = None
        corr_file = None
        metab_hmdb_anno_file = None
        gene_id2name = None
        if os.path.exists(self.asso_corr_tool.output_dir + "/asso.cluster_tree.xls"):
            assoc_tree = self.link_file("asso.cluster_tree.xls", "gene.cluster_tree.xls")
        if os.path.exists(self.asso_corr_tool.output_dir + "/metab.cluster_tree.xls"):
            metab_tree = self.link_file("metab.cluster_tree.xls", "metab.cluster_tree.xls")
            metab_tree_id = self.link_file("metab_id.cluster_tree.xls", "metab_id.cluster_tree.xls")
        if not os.path.exists(self.asso_corr_tool.output_dir + "/asso.cluster_tree.xls") or not os.path.exists(
                        self.asso_corr_tool.output_dir + "/metab.cluster_tree.xls"):
            listfile = corr_file
        # 将gene_id和gene_name的对应关系导入到主表当中
        if os.path.exists(os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")):
            gene_id2name = os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")
        trans_select_table = os.path.join(self.trans_tool.output_dir, "trans_select_table.xls")
        api_name.add_relation_corr(main_id=main_id, trans_tree=assoc_tree, metab_tree=metab_tree_id, gene_id2name=gene_id2name, trans_select_table=trans_select_table)
        # 导入绘图详情表
        corr_file = self.link_file("corr.xls", "corr.xls")
        p_file = self.link_file("pvalue.xls", "corr_pvalue.xls")
        metab_anno_file = self.option("metab_anno").prop['path']
        self.logger.info(self.option("metab_hmdb_anno"))
        if self.option("metab_hmdb_anno").is_set:
            metab_hmdb_anno_file = self.option("metab_hmdb_anno").prop['path']
        api_name.add_relation_corr_detail(main_id, corr_file, p_file, metab_anno=metab_anno_file, metab_hmdb_anno=metab_hmdb_anno_file) # 详情表中还要导入代谢物的kegg和hmdb的相关信息
        subclusters = glob.glob(self.asso_corr_tool.output_dir + '/*subcluster*')
        for each in subclusters:
            if os.path.exists(each):
                each = each.split("/")[-1]
                self.link_file(each, each)
        if os.path.exists(self.asso_corr_tool.output_dir + "/metab.kmeans_cluster.xls"):
            self.link_file("metab.kmeans_cluster.xls", "metab.kmeans_cluster.xls")
        if os.path.exists(self.asso_corr_tool.output_dir + "/trans.kmeans_cluster.xls"):
            self.link_file("trans.kmeans_cluster.xls", "gene.kmeans_cluster.xls")
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "relationcorr",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.asso_corr_tool.output_dir, oldfile)
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
            [".", "", "关联分析相关性结果文件夹", 0, "150069"],
            ["gene.cluster_tree.xls", "xls", "基因聚类树文件", 0, "150070"],
            ["metab.cluster_tree.xls", "xls", "代谢物聚类树文件", 0, "150028"],
            ["corr.xls", "xls", "基因与代谢物相关性系数表", 0, "150071"],
            ["corr_pvalue.xls", "xls", "基因与代谢物相关性系数P值表", 0, "150072"],
            ["metab.kmeans_cluster.xls", "xls", "代谢物相关性kmeans分类文件", 0, "150067"],
            ["gene.kmeans_cluster.xls", "xls", "基因相关性kmeans分类文件", 0, "150073"],
            ["metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
        ]
        regexps = [
            [r"metab\.subcluster_.*\.xls", "xls", "代谢物相关性kmeans各子类结果表", 0, "150068"],
            [r"gene\.subcluster_.*\.xls", "xls", "基因相关性kmeans各子类结果表", 0, "150074"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(RelationCorrWorkflow, self).end()