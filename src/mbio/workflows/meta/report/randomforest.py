# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
# last modify by shaohua.yuan in 20171103; 2018.04.23 by zhujuan 新增十折交叉验证、AUC验证和预测样本分类功能

""""""

import datetime
from biocluster.workflow import Workflow
import re
import os
import json
import shutil
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class RandomforestWorkflow(Workflow):
    """
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RandomforestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "predict_sample", "type": "infile", "format": "meta.otu.otu_table"},  # 用于预测丰度表样品分类文件
            {"name": "ntree", "type": "int", "default": 500},
            {"name": "problem_type", "type": "int", "default": 1},  # change default 2 to 1 by shaohua.yuan
            {"name": "method", "type": "string", "default": "CV"},
            {"name": "randomforest_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "group_id2", "type": 'string'},  # 用于“选择预测样本”
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "group_detail2", "type": "string"},  # 用于“选择预测样本”
            {"name": "env_id", "type": "string"},  # 环境因子表id
            {"name": "env_labs", "type": "string"},  # 选择的环境因子
            {"name": "norm_method", "type": "string", "default": ""},  # 新增的数据标准化方式
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.randomforest = self.add_tool("meta.beta_diversity.randomforest")
        self.output_dir = self.randomforest.output_dir

    def change_otuname(self, tablepath, file_name):
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        with open(tablepath, "r") as f, open(newtable, "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                specimen = re.subn("^.*; ", "", lines[0])[0]
                g.write(specimen + "\t" + lines[1])
        return newtable

    def run_randomforest(self):
        newtable = self.change_otuname(self.option('otutable').prop['path'], "otutable")
        options = {
            'otutable': newtable,
            'method': self.option('method'),
            'level': self.option('level'),
            'grouptable': self.option('grouptable'),
            'ntree': self.option('ntree'),
            'problem_type': self.option('problem_type'),
            'norm_method': self.option('norm_method')
        }
        if self.option('predict_sample').is_set:
            newtable = self.change_otuname(self.option('predict_sample').prop['path'], "predict_sample")
            options["predict_sample"] = newtable
        self.randomforest.set_options(options)
        self.randomforest.on('end', self.set_db)
        self.output_dir = self.randomforest.output_dir
        self.randomforest.run()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
            if self.option('method') == "CV":
                if os.path.exists(self.output_dir+"/AUC或者十折交叉验证结果图.pdf"):
                    os.rename(self.output_dir+"/AUC或者十折交叉验证结果图.pdf",self.output_dir+"/十折交叉验证结果图.pdf")
            else:
                if os.path.exists(self.output_dir + "/AUC或者十折交叉验证结果图.pdf"):
                    os.rename(self.output_dir + "/AUC或者十折交叉验证结果图.pdf", self.output_dir + "/AUC交叉验证结果图.pdf")
        save_params(self.output_dir, self.id)
        repaths = [
            [".", "", "RandomForest分析结果目录", 0, "110169"],
            ["./randomForest_confusion_subRF_table.xls", "xls", "随机森林计算出的分类结果(根据挑选的重要性变量)", 0, "110175"],
            ["./randomForest_confusion_table.xls", "xls", "随机森林计算出的分类结果(所有变量)", 0, "110171"],
            ["./randomForest_imptance_table.xls", "xls", "所有物种（变量）的重要性衡量值", 0, "110174"],
            ["./randomForest_10-fold_CV.xls", "xls", "不同物种(变量)数下十折交叉验证的分类错误率表格", 0, "110177"],
            ["./randomForest_AUC.xls", "xls", "不同物种(变量)数下随机森林的AUC表", 0, "110172"],
            ["./randomForest_predict.xls", "xls", "随机森林预测样品结果表", 0, "110173"],
            ["./randomForest_subRF_pcoa_sites.xls", "xls", "根据挑选出的重要物种(变量)构建随机森林的Pcoa坐标", 0, "110175"],
            ["./物种重要性柱形图.pdf", "pdf", "重要性排行前50的物种或环境因子", 0, ""],
            ["./十折交叉验证结果图.pdf", "pdf", "Random Forest验证结果图", 0, ""],
            ["./AUC交叉验证结果图.pdf", "pdf", "Random Forest验证结果图", 0, ""],
        ]
        regexps = [
            [r'randomForest_top.*_vimp.xls$', 'xls', '重要物种（变量）丰度表格', 0, "110170"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(RandomforestWorkflow, self).end()

    def set_db(self):
        api_randomforest = self.api.randomforest
        datadim = self.output_dir + '/randomForest_subRF_pcoa_sites.xls'
        datavip = self.output_dir + '/randomForest_imptance_table.xls'
        if self.option('method') == "AUC":
            datamethod = self.output_dir + '/randomForest_AUC.xls'
            api_randomforest.add_randomforest_evaluate(file_path=datamethod, table_id=self.option("randomforest_id"))
        elif self.option('method') == "CV":
            datamethod = self.output_dir + '/randomForest_10-fold_CV.xls'
            api_randomforest.add_randomforest_evaluate(file_path=datamethod, table_id=self.option("randomforest_id"))
        if self.option('predict_sample').is_set:
            predict = self.output_dir + '/randomForest_predict.xls'
            api_randomforest.add_randomforest_evaluate(file_path=predict, table_id=self.option("randomforest_id"))
        if not os.path.isfile(datadim):
            self.logger.error("找不到报告文件:{}".format(datadim))
            self.set_error("找不到报告文件", code="12703601")
        if not os.path.isfile(datavip):
            self.logger.error("找不到报告文件:{}".format(datavip))
            self.set_error("找不到报告文件", code="12703601")
        api_randomforest.add_randomforest_dim(file_path=datadim, table_id=self.option("randomforest_id"))
        api_randomforest.add_randomforest_vip(file_path=datavip, table_id=self.option("randomforest_id"))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("randomforest_id"), "sg_randomforest")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("randomforest_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "randomforest",
                "interaction": 1,
                "main_table": "sg_randomforest",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_randomforest()
        super(RandomforestWorkflow, self).run()
