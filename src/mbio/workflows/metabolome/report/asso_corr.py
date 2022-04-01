# -*- coding: utf-8 -*-
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
import pandas as pd



class AssoCorrWorkflow(Workflow):
    """
    相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AssoCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_id", "type": "string"},
            {"name": "asso_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "asso_col_row", "type": "string"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "coefficient", "type": "string", "default": "pearson"},
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "asso_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "metab_top", "type": "int", "default": 50},
            {"name": "asso_top", "type": "int", "default": 50},
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "metab_n_cluster", "type": "int", "default": ""},
            {"name": "asso_n_cluster", "type": "int", "default": ""},
            {"name": "asso_dist", "type": "string", "default": "euclidean"},
            {"name": "asso_cluster", "type": "string", "default": "complete"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
            #{"name": "new_old_names_map", "type": "string"}  ###将新名称改成原始名称后进行分析 20201030
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.metab_tool = self.add_tool("metabolome.select_table")
        self.asso_tool = self.add_tool("metabolome.select_table")
        self.asso_corr_tool = self.add_tool("metabolome.metabset.asso_corr")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Asso_cor workflow")
        self.select_tools = [self.metab_tool, self.asso_tool]
        self.on_rely(self.select_tools, self.run_asso_corr)
        self.asso_corr_tool.on('end', self.set_db)
        self.select_meta()
        self.select_asso()
        super(AssoCorrWorkflow, self).run()

    def select_meta(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "top": self.option("metab_top"),
            "st": "F",
            "group": self.option("group_table"),
            "select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": 0
        }
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_asso(self):
        self.logger.info("start profile!")
        exp_profile = self.option("asso_table").prop["path"]
        '''
        ##将新名称改成原始名称后进行分析 20201030
        data = pd.read_table(exp_profile, sep='\t', header=0)
        new_old_names_map = eval(self.option('new_old_names_map'))
        if self.option("asso_col_row") == "row":
            col1 = data.columns[0]
            for i in data.index:
                p_name = data.loc[i,col1]
                if p_name in new_old_names_map:
                    data.loc[i,col1] = new_old_names_map[p_name]
        else:
            data.rename(columns=new_old_names_map,inplace=True)

        new_exp = os.path.join(self.work_dir, 'change_names.asso.xls')
        data.to_csv(new_exp, sep='\t',index=False)

        #group_file = self.option("group").prop["path"]
        '''

        options = {
            "origin_table": exp_profile,
            "top": self.option("asso_top"),
            "st": "F",
            "group": self.option("group_table"),
            "merge": "nomerge",
            "group_method": 0
        }
        #if self.option("asso_col_row") == "row":
        options["trans"] = True
        self.asso_tool.set_options(options)
        self.asso_tool.run()

    def run_asso_corr(self):
        self.logger.info("start run asso_corr !")
        metab_table = self.metab_tool.option("select_table")
        asso_table = self.asso_tool.option("select_table")
        metab_abu = self.option("metab_table").prop["path"]
        #metab_des = metab_abu.replace("metab_abund.txt", "metab_desc.txt")
        metab_des = self.option("metab_desc").prop["path"]
        options = {
            'metab_table': metab_table,
            'asso_table': asso_table,
            'mct': self.option("metab_cluster_method"),
            'sct': self.option("asso_cluster_method"),
            'coefficient': self.option("coefficient"),
            'metab_trans': metab_des
        }
        if self.option("metab_cluster_method") == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("asso_cluster_method") == "hierarchy":
            options['scm'] = self.option("asso_cluster")
            options['scd'] = self.option("asso_dist")
        if self.option("metab_cluster_method") == "kmeans":
            options['metab_n_cluster'] = self.option("metab_n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("asso_cluster_method") == "kmeans":
            options['asso_n_cluster'] = self.option("asso_n_cluster")
            options['scd'] = self.option("asso_dist")
        self.asso_corr_tool.set_options(options)
        self.asso_corr_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("metabolome.asso_corr")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        assoc_tree = None
        metab_tree = None
        listfile = None
        metab_tree_id = None
        corr_file = self.link_file("corr.xls", "corr.xls")
        p_file = self.link_file("pvalue.xls", "corr_pvalue.xls")
        if os.path.exists(self.asso_corr_tool.output_dir + "/asso.cluster_tree.xls"):
            assoc_tree = self.link_file("asso.cluster_tree.xls", "asso.cluster_tree.xls")
        if os.path.exists(self.asso_corr_tool.output_dir + "/metab.cluster_tree.xls"):
            metab_tree = self.link_file("metab.cluster_tree.xls", "metab.cluster_tree.xls")
            metab_tree_id = self.link_file("metab_id.cluster_tree.xls", "metab_id.cluster_tree.xls")
        if not os.path.exists(self.asso_corr_tool.output_dir + "/asso.cluster_tree.xls") or not os.path.exists(
                        self.asso_corr_tool.output_dir + "/metab.cluster_tree.xls"):
            listfile = corr_file
        api_name.add_association_corr(main_id=main_id, assoc_tree=assoc_tree, metab_tree=metab_tree_id, listfile=listfile)
        api_name.add_association_corr_detail(main_id, corr_file, p_file)
        subclusters = glob.glob(self.asso_corr_tool.output_dir + '/*subcluster*')
        for each in subclusters:
            if os.path.exists(each):
                each = each.split("/")[-1]
                self.link_file(each, each)
        if os.path.exists(self.asso_corr_tool.output_dir + "/metab.kmeans_cluster.xls"):
            self.link_file("metab.kmeans_cluster.xls", "metab.kmeans_cluster.xls")
        if os.path.exists(self.asso_corr_tool.output_dir + "/asso.kmeans_cluster.xls"):
            self.link_file("asso.kmeans_cluster.xls", "asso.kmeans_cluster.xls")

        #add v3 ,主表增加asso_sample
        asso_sample = self.get_asso_sample()
        api_name.db['association_corr'].update({"_id":main_id}, {"$set": {"asso_sample": asso_sample}})
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "associationcorr",
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

    #获取2张表共有的样本 v3
    def get_asso_sample(self):
        asso_table = pd.read_table(self.option("asso_table").path, sep='\t', index_col=0)
        if self.option("asso_col_row") == 'col':
            asso_table_sample = asso_table.columns.tolist()
        else:
            asso_table_sample = asso_table.index.tolist()
        group_detail = json.loads(self.option("group_detail"))
        select_sample = []
        for g in group_detail:
            select_sample.extend(group_detail[g])
        asso_sample = set(asso_table_sample) & set(select_sample)
        return ','.join(asso_sample)


    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析相关性结果文件夹", 0, "150069"],
            ["asso.cluster_tree.xls", "xls", "关联数据树文件", 0, "150070"],
            ["metab.cluster_tree.xls", "xls", "代谢物聚类树文件", 0, "150028"],
            ["corr.xls", "xls", "关联数据与代谢物相关性系数表", 0, "150071"],
            ["corr_pvalue.xls", "xls", "关联数据与代谢物相关性系数P值表", 0, "150072"],
            ["metab.kmeans_cluster.xls", "xls", "代谢物相关性kmeans分类文件", 0, "150067"],
            ["asso.kmeans_cluster.xls", "xls", "代谢物kmeans分类文件", 0, "150073"],
            ["metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
        ]
        regexps = [
            [r"metab\.subcluster_.*\.xls", "xls", "代谢物相关性kmeans各子类结果表", 0, "150068"],
            [r"asso\.subcluster_.*\.xls", "xls", "关联数据相关性kmeans各子类结果表", 0, "150074"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(AssoCorrWorkflow, self).end()