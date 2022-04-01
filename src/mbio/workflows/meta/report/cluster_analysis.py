# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""群落组成分析模块"""
import os
import json
import gevent
import datetime
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mainapp.models.mongo.public.meta.meta import Meta
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class ClusterAnalysisWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ClusterAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "level", "type": "string", "default": "9"},  # 输入的OTU level
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "group_detail", "type": "string"}, # 输入的group_detail 示例如下
            # {"A":["578da2fba4e1af34596b04ce","578da2fba4e1af34596b04cf","578da2fba4e1af34596b04d0"],"B":["578da2fba4e1af34596b04d1","578da2fba4e1af34596b04d3","578da2fba4e1af34596b04d5"],"C":["578da2fba4e1af34596b04d2","578da2fba4e1af34596b04d4","578da2fba4e1af34596b04d6"]}
            #{"name": "method", "type": "string", "default": ""}  # 聚类方式， ""为不进行聚类
            {"name": "method", "type": "string", "default": ""},  #guanqing 20180411
            {"name": "combine_value", "type": "string", "default": ""} #guanqing 20180411
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        # self.matrix = self.add_tool("meta.beta_diversity.distance_calc") # 2016.12.1 zhouxuan
        # self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        group_table_path = os.path.join(self.work_dir, "group_table.xls")
        meta = Meta()
        meta._config = self.config  # 兼容不同mongo库版本
        self.group_table_path = meta.group_detail_to_table(self.option("group_detail"), group_table_path)

    # def check_options(self):  # 2016.12.1 zhouxuan
    #     if self.option('method') not in ['average', 'single', 'complete', ""]:
    #         raise OptionError('错误的层级聚类方式：%s' % self.option('method'))

    def run_sort_samples(self):
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path,
            "method": self.option("method"),   # guanqing 20180411
            "others": self.option("combine_value")  #guanqing 20180411
        })
        # if self.option("method") != "":  # 2016.12.1 zhouxuan
        #     self.sort_samples.on("end", self.run_matrix)
        # else:
        self.sort_samples.on("end", self.set_db)
        self.output_dir = self.sort_samples.output_dir  # modify by zhouxuan 2016.11.23
        self.sort_samples.run()

    # def run_matrix(self):  # 2016.12.1 zhouxuan
    #     trans_otu = os.path.join(self.work_dir, "otu.trans")
    #     self.sort_samples.option("out_otu_table").transposition(trans_otu)
    #     self.matrix.set_options({
    #         "method": "bray_curtis",
    #         "otutable": trans_otu
    #     })
    #     self.matrix.on('end', self.run_cluster)
    #     self.matrix.run()
    #
    # def run_cluster(self):
    #     options = {
    #         "dis_matrix": self.matrix.option('dis_matrix'),
    #         "linkage": self.option("method")
    #     }
    #     self.hcluster.set_options(options)
    #     self.hcluster.on('end', self.set_db)
    #     self.hcluster.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        # newick_id = ""
        # myParams = json.loads(self.sheet.params)  # 2016.12.1 zhouxuan
        # if self.option("method") != "":
        #     api_heat_cluster = self.api.heat_cluster
        #     name = "heat_cluster_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        #     newick_id = api_heat_cluster.create_newick_table(self.sheet.params, self.option("method"), myParams["otu_id"], name)
        #     self.hcluster.option("newicktree").get_info()
        #     api_heat_cluster.update_newick(self.hcluster.option("newicktree").prop['path'], newick_id)
        #     self.add_return_mongo_id("sg_newick_tree", newick_id, "", False)
        api_otu = self.api.cluster_analysis
        # new_otu_id = api_otu.add_sg_otu(self.sheet.params, self.option("input_otu_id"), None, newick_id)
        #api_otu.add_sg_otu_detail(self.sort_samples.option("out_otu_table").prop["path"], self.option("main_id"), self.option("input_otu_id")) 
        api_otu.add_sg_otu_detail(self.sort_samples.option("level_otu_table").prop["path"], self.option("main_id"), self.option("input_otu_id")) # guanqing.zou 20180411
        # self.add_return_mongo_id("sg_otu", new_otu_id)
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
                "submit_loc": "otu_group_analyse",
                "interaction": 1,
                "main_table": "sg_otu",
            })
            self.figsave.run()
        else:
            self.end()


    def end(self):   # modify by zhouxuan 2016.11.23
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种组成分析结果目录", 0, "110074"],
            ["taxa.table.xls", "xls", "各样本物种丰度结果表", 0, "110076"],  #modified by hongdongxuan 20170321
            ["taxa.precents.table.xls", "xls", "各样本物种相对丰度结果表", 0, "110075"],  #add by wangzhaoyue 2017.03.06
            ["群落柱形图.pdf", "pdf", "物种群落柱形图", 0, ""]
        ])
        super(ClusterAnalysisWorkflow, self).end()


    def run(self):
        self.run_sort_samples()
        super(ClusterAnalysisWorkflow, self).run()
