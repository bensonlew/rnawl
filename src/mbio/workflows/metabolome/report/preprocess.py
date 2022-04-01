# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import os


class PreprocessWorkflow(Workflow):
    """
    代谢集预处理
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PreprocessWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ana_method", "type": "string", "default": "GC"},
            {"name": "pos_rout", "type": "infile", "format": "metabolome.metab_table_dir"},
            {"name": "neg_rout", "type": "infile", "format": "metabolome.metab_table_dir"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "fillna", "type": "string", "default": "min"},
            {"name": "rsd", "type": "string", "default": "30%"},
            {"name": "norm", "type": "string", "default": "sum"},
            {"name": "sample_name", "type": "string"},
            {"name": "inner_ref", "type": "string"},
            {"name": "scale", "type": "string", "default": "log2"},
            {"name": "log", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "rm_nan", "type":"float", "default":50},  #20190605
            {"name": "fill_type", "type":"string","default":"all"}, #v3 202003
            {"name" : "raw_exp_id", "type":"string","default":""},
            {"name": "run_cv_raw", "type": "string", "default": "False"}, # 决定是否计算原始的cv 和统计
            {"name": "task_version", "type" : "string", "default":"3.0"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.preprocess = self.add_tool("metabolome.preprocess")
        # self.output_dir = self.preprocess.output_dir

    def run(self):
        options = {
            "ana_method": self.option('ana_method'),
            "pos_rout": self.option("pos_rout"),
            "fillna": self.option("fillna"),
            "rsd": self.option("rsd"),
            "norm": self.option("norm"),
            "scale": self.option("scale"),
            "interactive": "True",
            "group_table": self.option("group_table"),
            "rm_nan" : self.option("rm_nan"),  #20190605
            "fill_type" : self.option("fill_type") #v3 202003
        }
        if self.option("ana_method") == "LC":
            options["neg_rout"] = self.option("neg_rout")
        if self.option("norm") == "sample":
            options["sample_name"] = self.option("sample_name")
        elif self.option("norm") == "inner":
            options["inner_ref"] = self.option("inner_ref")
        if self.option("scale") == "defined":
            options["log"] = self.option("log")
        if self.option("run_cv_raw")=="True" and self.option("task_version") !="1.0":
            options['raw_cv'] = 'T'
        else:
            options['raw_cv'] = 'F'
        options['task_version'] = self.option("task_version")
        self.preprocess.set_options(options)
        self.preprocess.on('end', self.set_db)
        self.preprocess.run()
        super(PreprocessWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        preprocess_api = self.api.api('metabolome.preprocess')
        table_path = self.preprocess.option("pos_out").path
        if self.option('ana_method') == "LC":
            table_path += "," + self.preprocess.option("neg_out").path
            table_path += ',' + self.preprocess.option("mix_out").path

        if self.option("task_version") !="1.0":
            preprocess_api.add_metab_table_detail(self.option('main_table_id'), table_path, raw_exp_id =self.option("raw_exp_id"))
            if self.option("run_cv_raw")=="True":
                #preprocess_api.add_metab_table_detail(self.option("raw_exp_id"), table_path, raw_exp_id =self.option("raw_exp_id"))
                preprocess_api.interaction_ori_cv(table_path,self.option("raw_exp_id"))
        else:
            preprocess_api.add_metab_table_detail(self.option('main_table_id'), table_path)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "exp",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.preprocess.output_dir))
        result_dir = self.add_upload_dir(self.preprocess.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "预处理结果文件夹", 0, "150001"],
            ["pos/", "", "预处理阳离子结果", 0, "150002"],
            ["pos/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["pos/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["neg/", "", "预处理阴离子结果", 0, "150005"],
            ["neg/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["neg/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["mix/", "", "预处理混合离子结果", 0, "150006"],
            ["mix/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["mix/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"]
        ])
        super(PreprocessWorkflow, self).end()
