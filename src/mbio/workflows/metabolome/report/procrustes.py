# -*- coding: utf-8 -*-
# __author__ = 'linmeng.liu'
# last_modifiy = modified 2018.0816

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



class ProcrustesWorkflow(Workflow):
    """
    Procrustes description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProcrustesWorkflow, self).__init__(wsheet_object)
        '''
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        '''
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_set_id", "type": "string"},
            {"name": "asso_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "asso_col_row", "type": "string"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "method", "type": "string", "default": "pca"},
            {"name": "metab_dist", "type": "string", "default": "bray_curtis"},
            {"name": "asso_dist", "type": "string", "default": "bray_curtis"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
            #{"name": "new_old_names_map", "type": "string"}  ###将新名称改成原始名称后进行分析 20201030
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.procrustes = self.add_module("metabolome.procrustes")
        self.metab_tool = self.add_tool("metabolome.select_table")
        self.asso_tool = self.add_tool("metabolome.select_table")
        self.dump_tool = self.api.api("metabolome.procrustes")

    def run(self):
        # self.start_listener(); self.fire("start") # if you have no tools, you should use this line
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Procrustes workflow")
        self.select_tools = [self.metab_tool, self.asso_tool]
        self.on_rely(self.select_tools, self.run_procrustes)
        self.procrustes.on('end', self.set_db)
        self.select_metab()
        self.select_asso()
        super(ProcrustesWorkflow, self).run()


    def select_metab(self):
        self.logger.info("start select metab file!")
        exp_profile = self.option("metab_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "st": "F",
            "group": self.option("group_table"),
            "select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": 0
        }
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_asso(self):
        self.logger.info("start select asso profile!")
        exp_profile = self.option("asso_table").prop["path"]

        '''
        #group_file = self.option("group").prop["path"]
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
        '''

        options = {
            "origin_table": exp_profile,
            "st": "F",
            "group": self.option("group_table"),
            "merge": "nomerge",
            "group_method": 0
        }
        #if self.option("asso_col_row") == "row":
        options["trans"] = True
        self.asso_tool.set_options(options)
        self.asso_tool.run()


    def run_procrustes(self):
        self.logger.info("start run procrustes !")
        metab_table = self.metab_tool.option("select_table")
        asso_table = self.asso_tool.option("select_table")
        options = dict(
            metab_table = metab_table,
            metab_dist = self.option('metab_dist'),
            asso_table = asso_table,
            asso_dist = self.option('asso_dist'),
            method = self.option('method')
        )
        self.procrustes.set_options(options)
        self.procrustes.run()

    def set_db(self):
        """
        dump data to db
        add result info
        """
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('ERROR: wrong format of main_id, ObjectId or string expected!')
        reg1 = os.path.join(self.procrustes.output_dir,"*transformed_reference.txt")
        reg2 = os.path.join(self.procrustes.output_dir,"*transformed_q*.txt")
        oldfile1 = glob.glob(reg1)
        oldfile2 = glob.glob(reg2)
        self.logger.info("file_path1: %s" % oldfile1)
        self.logger.info("file_path2: %s" % oldfile2)
        summary = os.path.join(self.procrustes.output_dir, "procrustes_results.txt")
        #summary = self.link_file("procrustes_results.txt", "procrustes_summary.txt",1)
        transform_ref = self.link_file(oldfile1[0], "metab.transformed_reference.txt",0)
        transform_query = self.link_file(oldfile2[0], "asso.transformed_query.txt",0)
        self.dump_tool.add_procrustes_detail(summary, main_id)
        self.dump_tool.add_procrustes_graph(transform_ref, transform_query, main_id)
        #add v3 ,主表增加asso_sample
        asso_sample = self.get_asso_sample()
        self.dump_tool.db['metabset_procrustes'].update({"_id":main_id}, {"$set": {"asso_sample": asso_sample}})
        # self.option("save_pdf", 0)
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
            [".", "", "关联分析普氏分析结果文件夹", 0, "150079"],
            ["metab.transformed_reference.txt", "txt", "代谢物拓扑转换后结果", 0, "150080"],
            ["asso.transformed_query.txt", "txt", "关联数据拓扑转换后结果", 0, "150081"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        super(ProcrustesWorkflow, self).end()
