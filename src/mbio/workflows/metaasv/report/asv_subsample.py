# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

"""otu样本序列数抽平"""
from biocluster.workflow import Workflow
import os
import json
import shutil
from mainapp.models.mongo.public.meta.meta import Meta
from mbio.packages.metaasv.common_function import link_dir,filter_asv_set


class AsvSubsampleWorkflow(Workflow):
    """
    Metaasv ASV分类学分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AsvSubsampleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的asv id
            {"name": "filter_json", "type": "string", "default": ""},  # 输入的json文件
            {"name": "size", "type": "string", "default": "min"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入Group表
            {"name": "level", "type": "string", "default": "9"},
            {"name": "update_info", "type": 'string'},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_id", "type": 'string', "default": ''}, ##传入主表
            {"name": "set_list", "type": 'string', "default": ""} ##传入ASV集情况
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filter_otu = self.add_tool("metaasv.filter_asv")
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.subsample = self.add_tool("meta.otu.sub_sample")
        self.tax_stat = self.add_tool("metaasv.asv_taxon_stat")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        # group_table_path = os.path.join(self.work_dir, "group_table.xls")
        # self.group_table_path = Meta().group_detail_to_table(self.option("group_detail"), group_table_path)
        # self.table2db = ""

    def filter_asv(self):
        """
        功能：过滤掉基因集
        """
        self.logger.info("开始过滤掉基因集")
        out_asv_path = os.path.join(self.work_dir, "input_asv.xls")
        if os.path.exists(out_asv_path):
            os.remove(out_asv_path)
        filter_asv_set(self.option("in_otu_table").prop["path"], self.option("set_list"), out_asv_path)
        self.logger.info("过滤基因集完成！")

    def run_sort_samples(self):
        """
        运行tool 得到排序后的结果
        """
        if self.option("set_list") in [""]:
            asv_path = self.option("in_otu_table").prop['path']
        else:
            asv_path = os.path.join(self.work_dir, "input_asv.xls")
        self.sort_samples.set_options({
            "in_otu_table": asv_path,
            "group_table": self.option("group_table")
        })
        if self.option("filter_json") not in ["", "[]"]:
            self.sort_samples.on("end", self.run_filter_otu)
        elif self.option("size") != "":
            self.sort_samples.on("end", self.run_subsample)
        else:
            self.sort_samples.on("end", self.result_format)
        self.sort_samples.run()

    def run_filter_otu(self):
        """
        运行tool filter_asv
        """
        self.logger.info("filter_json: {}".format(self.option("filter_json")))
        filter_json = json.loads(self.option("filter_json"))
        self.filter_otu.set_options({
            "in_otu_table": self.sort_samples.option("out_otu_table"),
            "filter_json": json.dumps(filter_json),
            "group_table": self.option("group_table")
        })
        if self.option("size") != "":
            self.filter_otu.on("end", self.run_subsample)
        else:
            self.filter_otu.on("end", self.result_format)
        self.filter_otu.run()

    def run_subsample(self):
        """
        运行tool sub_sample
        """
        if self.option("filter_json") not in ["", "[]"]:
            num_lines = sum(1 for line in open(self.filter_otu.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！")
            self.subsample.set_options({
                "in_otu_table": self.filter_otu.option("out_otu_table"),
                "size": self.option("size")
            })
        else:
            self.subsample.set_options({
                "in_otu_table": self.sort_samples.option("out_otu_table"),
                "size": self.option("size")
            })
        self.subsample.on("end", self.result_format)
        self.subsample.run()


    def result_format(self):
        """
        统计抽平结果 运行 otu_taxon_stat
        """
        if self.option("filter_json") not in ["", "[]"]:
            num_lines = sum(1 for line in open(self.filter_otu.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！")
        if self.option("size") != "":
            num_lines = sum(1 for line in open(self.subsample.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.logger.error("经过抽平之后的OTU表是空的，可能是因为进行物种筛选之后导致某些样本的序列数为0，然后按该样本的序列数进行了抽平！")
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！")
            final_file = self.subsample.option("out_otu_table")
        elif self.option("filter_json") not in ["", "[]"]:
            final_file = self.filter_otu.option("out_otu_table")
        else:
            final_file = self.sort_samples.option("out_otu_table")
        self.table2db = final_file.prop['path']
        self.tax_stat.set_options({"sub_otu_table": final_file})
        self.tax_stat.on("end", self.set_db)
        self.tax_stat.run()

    def set_db(self):
        """
        保存结果otu表到mongo数据库中
        """
        link_dir(self.tax_stat.output_dir, self.output_dir)
        api_otu = self.api.api("metaasv.sub_sample")
        self.logger.info("开始将信息导入asv_detail表和asv_specimen表中")
        otu_detail_level = self.output_dir + "/tax_summary_a/asv_taxon_asv.full.xls"
        api_otu.add_sg_otu_detail(self.table2db, self.option("input_otu_id"), self.option('main_id'))
        api_otu.add_sg_otu_detail_level(otu_detail_level, self.option('main_id'), self.option("level"), database=True)
        api_otu.add_sg_otu_seq_summary(self.table2db,self.option('main_id'))
        self.add_return_mongo_id("asv", self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ASV物种分类学分析结果文件", 0, ""],
            ["./tax_summary_a", "meta.otu.tax_summary_dir", "各分类学水平样本序列数统计表", 0, ""],
            ["./tax_summary_r", "meta.otu.tax_summary_dir", "各分类学水平样本序列数相对丰度百分比统计表", 0, ""],
            ["./ASV_taxon.biom", "meta.otu.biom", "biom格式的OTU物种分类统计表", 0, ""],
            ["./ASV_taxon.xls", "meta.otu.otu_table", "OTU物种分类统计表", 0, ""],
            ["./ASV_summary.xls", "meta.otu.otu_table", "基于 OTU 数量的统计", 0, ""]
        ])
        result_dir.add_regexp_rules([
            ["tax_summary_a/.+\.biom$", "meta.otu.biom", "ASV表的biom格式的文件(absolute)", 0, ""],
            ["tax_summary_a/.+\.xls$", "xls", "单级物种分类统计表(absolute)", 0, ""],
            ["tax_summary_a/.+\.full\.xls$", "xls", "多级物种分类统计表(absolute)", 0, ""],
            ["tax_summary_r/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件", 0, ""],
            ["tax_summary_r/.+\.xls$", "xls", "单级物种分类统计表", 0, ""],
            ["tax_summary_r/.+\.full\.xls$", "xls", "多级物种分类统计表", 0, ""]
            ])
        super(AsvSubsampleWorkflow, self).end()

    def run(self):
        if self.option("set_list") not in [""]:
            self.filter_asv()
        self.run_sort_samples()
        super(AsvSubsampleWorkflow, self).run()
