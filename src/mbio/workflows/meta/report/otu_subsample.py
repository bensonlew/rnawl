# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""otu样本序列数抽平"""
from biocluster.workflow import Workflow
import os
import shutil
from mainapp.models.mongo.public.meta.meta import Meta
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class OtuSubsampleWorkflow(Workflow):

    """
    报告中调用otu抽平时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OtuSubsampleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "filter_json", "type": "string", "default": ""},  # 输入的json文件
            {"name": "size", "type": "string", "default": "min"},
            {"name": "group_detail", "type": "string"},
            {"name": "level", "type": "string", "default": "9"},
            {"name": "output_otu_id", "type": "string"},  # 结果的otu id
            {"name": "update_info", "type": 'string'},
            {"name": "main_id", "type": 'string', "default": ''},
            {"name": "params", "type": 'string', "default": ''},
            {"name": "group_id", "type": 'string'},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # PDF保存
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filter_otu = self.add_tool("meta.otu.filter_otu")
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.subsample = self.add_tool("meta.otu.sub_sample")
        self.tax_stat = self.add_tool("meta.otu.otu_taxon_stat")
        group_table_path = os.path.join(self.work_dir, "group_table.xls")
        meta = Meta()
        meta._config = self.config  # 兼容不同mongo库版本
        self.group_table_path = meta.group_detail_to_table(self.option("group_detail"), group_table_path)
        self.table2db = ""

    def run_sort_samples(self):
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path
        })
        if self.option("filter_json") not in ["", "[]"]:
            self.sort_samples.on("end", self.run_filter_otu)
        elif self.option("size") != "":
            self.sort_samples.on("end", self.run_subsample)
        else:
            self.sort_samples.on("end", self.result_format)
        self.sort_samples.run()

    def run_filter_otu(self):
        self.filter_otu.set_options({
            "in_otu_table": self.sort_samples.option("out_otu_table"),
            "filter_json": self.option("filter_json")
        })
        if self.option("size") != "":
            self.filter_otu.on("end", self.run_subsample)
        else:
            self.filter_otu.on("end", self.result_format)
        self.filter_otu.run()

    def run_subsample(self):
        if self.option("filter_json") not in ["", "[]"]:
            num_lines = sum(1 for line in open(self.filter_otu.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！", code="12703001")
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


    def result_format(self):  # add by zhengyuan
        """
        统计抽平结果
        """
        if self.option("filter_json") not in ["", "[]"]:
            num_lines = sum(1 for line in open(self.filter_otu.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！", code="12703001")
        if self.option("size") != "":
            num_lines = sum(1 for line in open(self.subsample.option("out_otu_table").prop["path"]))
            if num_lines < 2:
                self.logger.error("经过抽平之后的OTU表是空的，可能是因为进行物种筛选之后导致某些样本的序列数为0，然后按该样本的序列数进行了抽平！")
                self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！", code="12703001")
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
        lst = os.listdir(self.tax_stat.output_dir)
        for l in lst:
            path = os.path.join(self.tax_stat.output_dir, l)
            final_path = os.path.join(self.output_dir, l)
            if os.path.isdir(path):
                shutil.copytree(path, final_path)
            else:
                shutil.copy(path, final_path)
        api_otu = self.api.sub_sample
        # output_otu_id = api_otu.add_sg_otu(self.sheet.params, self.option("size"), self.option("input_otu_id"))
        #if not os.path.isfile(final_file):
        #    raise Exception("找不到报告文件:{}".format(final_file))
        self.logger.info("开始将信息导入sg_otu_detail表和sg_otu_specimen表中")
        api_otu.add_sg_otu_detail(self.table2db, self.option("input_otu_id"), self.option('main_id'))
        api_otu.add_sg_otu_detail_level(self.table2db, self.option('main_id'), self.option("level"))
        api_otu.add_sg_otu_seq_summary(self.table2db,self.option('main_id'))   #guanqing.zou 20180505
        self.add_return_mongo_id("sg_otu", self.option('main_id'))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_otu")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "otu_statistic",
                "interaction": 1,
                "main_table": "sg_otu",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            if os.path.exists(self.output_dir + "/OTU_Rank_Abundance曲线图.pdf"):
                os.remove(self.output_dir + "/OTU_Rank_Abundance曲线图.pdf")
            if os.path.exists(self.figsave.output_dir + "/9_Rank_Abundance曲线图.pdf"):
                os.link(self.figsave.output_dir + "/9_Rank_Abundance曲线图.pdf",self.output_dir + "/OTU_Rank_Abundance曲线图.pdf")
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "OTU分类统计结果目录", 0, "110013"],
            ["./tax_summary_a", "meta.otu.tax_summary_dir", "各分类学水平样本序列数统计表", 0, "110019"],
            ["./tax_summary_r", "meta.otu.tax_summary_dir", "各分类学水平样本序列数相对丰度百分比统计表", 0, "110015"],
            ["./otu_taxon.biom", "meta.otu.biom", "biom格式的OTU物种分类统计表", 0, "110024"],
            ["./otu_taxon.xls", "meta.otu.otu_table", "OTU物种分类统计表", 0, "110025"],
            ["./otu_summary.xls", "meta.otu.otu_table", "基于 OTU 数量的统计", 0, "110023"],
            ["./OTU_Rank_Abundance曲线图.pdf", "pdf", "各样本的Rank_Abundance曲线图", 0, ""]
            #["./otu_taxon.subsample.xls", "xls", "抽平后的OTU表格"]  # add by hongdongxuan 20170324
        ])
        result_dir.add_regexp_rules([
            ["tax_summary_a/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件(absolute)", 0, "110021"],
            ["tax_summary_a/.+\.xls$", "xls", "单级物种分类统计表(absolute)", 0, "110022"],
            ["tax_summary_a/.+\.full\.xls$", "xls", "多级物种分类统计表(absolute)", 0, "110020"],
            ["tax_summary_r/.+\.biom$", "meta.otu.biom", "OTU表的biom格式的文件", 0, "110018"],  # add by zhouxuan (3 line) 20161129
            ["tax_summary_r/.+\.xls$", "xls", "单级物种分类统计表", 0, "110017"],
            ["tax_summary_r/.+\.full\.xls$", "xls", "多级物种分类统计表", 0, "110016"]
            ])
        # result_dir.add_regexp_rules([
        #     ['\.subsample\.', 'meta.otu.otu_table', "抽平后的otu表格"]   # modified by hongdongxuan 20170324
        # ])
        super(OtuSubsampleWorkflow, self).end()

    def run(self):
        self.run_sort_samples()
        super(OtuSubsampleWorkflow, self).run()
