# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""Roc分析"""

from biocluster.workflow import Workflow
import os
from bson.objectid import ObjectId
import types
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class MgRocWorkflow(Workflow):
    """
    报告中调用Roc分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MgRocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "level_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "level_type", "type": "string"},
            {"name": "level_type_name", "type": "string"},
            {"name": "lowest_level", "type": "string", "default": ""},
            {"name": "confidence_interval", "type": "string"},
            {"name": "top_num", "type": "int"},
            {"name": "group_id", "type": "string"},
            {"name": "cal_method", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.abundance = self.add_tool("meta.create_abund_table")
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.roc = self.add_tool("statistical.roc")

    def run_abundance(self):
        options = {
            'anno_table': self.option('anno_table'),
            'geneset_table': self.option("geneset_table"),
            'level_type': self.option('level_type'),
            'level_type_name': self.option('level_type_name'),
            'gene_list': self.option('gene_list'),
            'lowest_level': self.option('lowest_level')
        }
        self.abundance.set_options(options)
        self.abundance.on("end", self.run_sort_samples)
        self.abundance.run( )

    def run_sort_samples(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group"),
            "sample_del_warn": "T"
        })
        self.sort_samples.on("end", self.run_roc)
        self.sort_samples.run()

    def run_roc(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        options={
            "abu_table": abund_table,
            "group_table": self.option("group"),
            "method": self.option("cal_method"),
            "confidence_interval":self.option("confidence_interval"),
            "top_n":self.option("top_num"),
        }
        self.logger.info(options)
        self.roc.set_options(options)
        self.roc.on('end', self.set_db)
        self.roc.run()


    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.roc.output_dir))
        result_dir = self.add_upload_dir(self.roc.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录", 0, "121001"],
            ["./roc_curve_smooth.xls", "xls", "ROC受试者工作特征曲线平滑处理后数据", 0, "121005"],
            ["./roc_curve.xls", "xls", "ROC受试者工作特征曲线数据", 0, "121002"],
            ["./roc_interval.xls", "xls", "ROC受试者工作特征曲线的置信区间", 0, "121006"],
            ["./roc_auc.xls", "xls", "ROC受试者工作特征曲线-AUC VALUE", 0, "121003"],
            ["./roc_auc_smooth.xls", "xls", "ROC受试者工作特征曲线-AUC VALUE平滑处理后数据", 0, "121004"],
            ["./best_loc.xls", "xls", "ROC受试者工作特征曲线最佳cut off点", 0, "121007"],
            ["roc.pdf", "pdf", "ROC曲线"]
        ])
        super(MgRocWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        api_roc = self.api.api("metagenomic.roc")
        main_id = self.option("main_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12803001")
        self.logger.info(main_id)
        if os.path.exists(self.roc.output_dir + '/roc_curve_smooth.xls'):
            api_roc.add_roc_curve(main_id,'smooth',self.roc.output_dir + '/roc_curve_smooth.xls')

        if os.path.exists(self.roc.output_dir + '/roc_curve.xls'):
            api_roc.add_roc_curve(main_id,'',self.roc.output_dir + '/roc_curve.xls')
        else:
            self.set_error('roc_curve.xls file is not exists！', code="12803003")
        if os.path.exists(self.roc.output_dir + '/roc_interval.xls'):
            api_roc.add_roc_interval(main_id,self.roc.output_dir + '/roc_interval.xls')

        if os.path.exists(self.roc.output_dir + '/roc_auc.xls'):
            api_roc.add_roc_auc(main_id,self.roc.output_dir + '/roc_auc.xls','')
        else:
            self.set_error('AUC.xls file is not exists！', code="12803005")
        if os.path.exists(self.roc.output_dir + '/roc_auc_smooth.xls'):
            api_roc.add_roc_auc(main_id,self.roc.output_dir + '/roc_auc_smooth.xls','smooth')

        if os.path.exists(self.roc.output_dir + '/best_loc.xls'):
            api_roc.add_roc_best_loc(main_id,self.roc.output_dir + '/best_loc.xls')
        else:
            self.set_error('best_loc.xls file is not exists！', code="12803007")
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "roc")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_id"), "roc")
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

    def run(self):
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            self.run_abundance()
        super(MgRocWorkflow, self).run()