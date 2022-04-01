# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import re
import pandas as pd
import os, glob


class RelationO2plsWorkflow(Workflow):
    """
    O2plsda分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationO2plsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},  # 代谢物desc表
            {"name": "group_detail", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "trans_exp_main_id", "type": "string"},  # 转录表达量表id
            {"name": "scale", "type": "string", 'default': "UV"},
            {"name": "log10", "type": "bool", "default": False},
            {'name': 'update_info', 'type': 'string'},
            {"name": "main_table_id", "type": "string"},
            {"name": "sort_name", "type": "string", 'default': "yes"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.metab_tool = self.add_tool("metabolome.select_table")
        self.trans_tool = self.add_tool("metabolome.relation.trans_select_table")
        self.o2plsda = self.add_tool("tool_lab.o2plsda")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.select_tools = [self.metab_tool, self.trans_tool]
        self.on_rely(self.select_tools, self.run_o2plsda)
        self.o2plsda.on('end', self.set_db)
        self.select_metab_table()
        self.select_trans_table()
        super(RelationO2plsWorkflow, self).run()

    def sort_abund_table(self, exp_table, desc_table):
        df_exp = pd.read_table(exp_table, "\t", index_col=0)
        df_desc = pd.read_table(desc_table, "\t", index_col=0)
        total_index = df_desc.index.tolist()
        for i in df_desc.index:
            if re.match(r'^metab_\d+$', df_desc.loc[i, "Metabolite"]):
                total_index.remove(i)
        df_sort = df_exp.ix[total_index, :]
        df_sort = df_sort.reset_index()
        has_name_exp_profile = self.work_dir + "/hasname_exp_table.txt"
        df_sort.to_csv(has_name_exp_profile, "\t", index=False)
        return has_name_exp_profile

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
        if self.option("sort_name") == "yes":
            has_name_exp_profile = self.sort_abund_table(exp_profile, self.option('metab_desc').prop["path"])
            options["origin_table"] = has_name_exp_profile
        elif self.option("sort_name") == "no":
            options["origin_table"] = exp_profile
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_trans_table(self):  #筛选转录表
        options = {
            "group": self.option("group_detail"),
            "trans_exp_main_id": self.option("trans_exp_main_id"),
            "task_id" : "_".join(self._sheet.id.split("_")[0:2])
        }
        self.trans_tool.set_options(options)
        self.trans_tool.run()

    def run_o2plsda(self):
        metab_file = self.metab_tool.output_dir + "/select_table.xls"
        trans_file = self.trans_tool.output_dir + "/trans_select_table.xls"
        self.rename_name(metab_file, self.work_dir + "/x_data.name", self.work_dir + "/x_data.xls", "metab")
        self.rename_name(trans_file, self.work_dir + "/y_data.name", self.work_dir + "/y_data.xls", "gene")
        self.o2plsda.set_options({
            'x_list': self.work_dir + "/x_data.name",
            'y_list': self.work_dir + "/y_data.name",
            'x_data': self.work_dir + "/x_data.xls",
            'y_data': self.work_dir + "/y_data.xls",
            'group': self.option("group_detail").path,
            'oxoy': "Metabolome;Transcript",
            'scale': self.option('scale'),
            "x_method": "1",
            'y_method': "1"
        })
        self.o2plsda.run()

    def rename_name(self, file ,file2, file3, type):
        with open(file, "r") as f, open(file2, "w+") as g, open(file3, "w+") as s:
            lines = f.readlines()
            num = 0
            s.write(lines[0])
            for line in lines[1:]:
                num += 1
                lin = line.strip().split("\t")
                g.write("{}\t{}\n".format(type + str(num), lin[0]))
                lin[0] = type + str(num)
                s.write("\t".join(lin) + "\n")

    def set_db(self):
        self.api_o2plsda = self.api.api('metabolome.relation_o2pls')
        path1 = self.o2plsda.output_dir + "/O2PLS_Loadings.xls"
        path2 = self.o2plsda.output_dir + "/O2PLS_Scores.xls"
        path3 = self.o2plsda.output_dir + "/fit_summary.xls"
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('ERROR: wrong format of main_id, ObjectId or string expected!')
        metab_desc = self.option('metab_desc').prop["path"]
        gene_id2name = None
        if os.path.exists(os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")):
            gene_id2name = os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")
        self.api_o2plsda.add_o2plsda_detail(main_id, path1, path2, path3, "Metabolome;Transcript", metab_desc=metab_desc, trans_desc=gene_id2name)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.o2plsda.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "O2PLS结果输出目录"],
            ["./fit_summary.xls", "xls", "O2PLS统计情表"],
            ["./O2PLS_Loadings.xls", "xls", "O2PLS分析Loadings信息表"],
            ["./O2PLS_Scores.xls", "xls", "O2PLS分析Scores统计表"],
        ])
        super(RelationO2plsWorkflow, self).end()