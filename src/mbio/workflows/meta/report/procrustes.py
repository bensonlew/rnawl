# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os, glob
import re
import shutil
import json
import types
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


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
            {"name": "otu_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "otu_id", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "asso_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "method", "type": "string", "default": "pca"},
            {"name": "otu_dist", "type": "string", "default": "bray_curtis"},
            {"name": "asso_dist", "type": "string", "default": "bray_curtis"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.procrustes = self.add_module("meta.procrustes")
        self.otu_tool = self.add_tool("sequence.select_table")
        self.asso_tool = self.add_tool("sequence.select_table")
        self.dump_tool = self.api.api("procrustes")

    def run(self):
        # self.start_listener(); self.fire("start") # if you have no tools, you should use this line
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Procrustes workflow")
        self.select_tools = [self.otu_tool, self.asso_tool]
        self.on_rely(self.select_tools, self.run_procrustes)
        self.procrustes.on('end', self.set_db)
        self.select_otu()
        self.select_asso()
        super(ProcrustesWorkflow, self).run()


    def select_otu(self):
        self.logger.info("start select otu file!")
        exp_profile = self.option("otu_table").prop["path"]
        options = {
            "origin_table": exp_profile,
            #"select_genes": self.option("metab_set_table"),
            "st": "F",
            "group": self.option("group_table").path,
            #"select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": 0
        }
        self.otu_tool.set_options(options)
        self.otu_tool.run()

    def select_asso(self):
        self.logger.info("start select asso profile!")
        exp_profile = self.option("asso_table").prop["path"]
        #group_file = self.option("group").prop["path"]
        options = {
            "origin_table": exp_profile,
            "st": "F",
            "group": self.option("group_table").path,
            "merge": "nomerge",
            "group_method": 0
        }
        #if self.option("asso_col_row") == "row":
        options["trans"] = True
        self.asso_tool.set_options(options)
        self.asso_tool.run()


    def run_procrustes(self):
        self.logger.info("start run procrustes !")
        otu_table = self.otu_tool.option("select_table")
        #asso_table = self.asso_tool.option("select_table")
        asso_table = self.work_dir + "/" + "asso.xls"
        with open(self.asso_tool.option("select_table").prop["path"]) as f,open(asso_table,"w") as t:
            data = f.readlines()
            head = data[0]
            if head.startswith("\t"):
                t.write("asso_table" + head)
            else:
                t.write(head)
            for i in data[1:]:
                t.write(i)
        options = dict(
            otu_table = otu_table,
            otu_dist = self.option('otu_dist'),
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
                self.set_error('ERROR: wrong format of main_id, ObjectId or string expected!', code="12705001")
        reg1 = os.path.join(self.procrustes.output_dir,"*transformed_reference.txt")
        reg2 = os.path.join(self.procrustes.output_dir,"*transformed_q*.txt")
        oldfile1 = glob.glob(reg1)
        oldfile2 = glob.glob(reg2)
        self.logger.info("file_path1: %s" % oldfile1)
        self.logger.info("file_path2: %s" % oldfile2)
        #summary = os.path.join(self.procrustes.output_dir, "procrustes_results.txt")
        summary = self.link_file("procrustes_results.txt", "Procrustes.txt",1)
        transform_ref = self.link_file(oldfile1[0], "metab.transformed_reference.txt",0)
        transform_query = self.link_file(oldfile2[0], "asso.transformed_query.txt",0)
        self.dump_tool.add_procrustes_detail(summary, main_id)
        self.dump_tool.add_procrustes_graph(transform_ref, transform_query, main_id)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "sg_procrustes")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "procrustes",
                "interaction": 1,
                "main_table": "sg_procrustes",
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

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        if os.path.exists(self.output_dir+"/运行参数.txt"):
            os.rename(self.output_dir+"/运行参数.txt",self.output_dir+"/运行参数1.txt")
            with open(self.output_dir+"/运行参数1.txt") as f,open(self.output_dir+"/运行参数.txt","w") as t:
                data = f.readlines()
                for i in data:
                    if "bray_curtis" in i:
                        pass
                    else:
                        t.write(i)
            os.remove(self.output_dir+"/运行参数1.txt")
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析普氏分析结果文件夹",0,"110262"],
            ["metab.transformed_reference.txt", "txt", "代谢物拓扑转换后结果",0,"110263"],
            ["asso.transformed_query.txt", "txt", "关联数据拓扑转换后结果",0,"110264"],
            ["Procrustes.txt","txt","Procrustes分析结果表",0,""],
            ["Procrustes分析结果图.pdf","pdf", "Procrustes分析结果图", 0,""]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        super(ProcrustesWorkflow, self).end()
