# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""otu 样本和并和抽平模块"""
import os
import shutil
from biocluster.workflow import Workflow
from mainapp.models.mongo.public.meta.meta import Meta
from mbio.packages.meta.save_params import save_params


class OtuStaticWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OtuStaticWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "level", "type": "string", "default": "9"},  # 输入的OTU level，这个分析中，不会输入， 所以值只会是default的9
            {"name": "group_detail", "type": "string"},  # 输入的group_detail 示例如下
            # {"A":["578da2fba4e1af34596b04ce","578da2fba4e1af34596b04cf","578da2fba4e1af34596b04d0"],"B":["578da2fba4e1af34596b04d1","578da2fba4e1af34596b04d3","578da2fba4e1af34596b04d5"],"C":["578da2fba4e1af34596b04d2","578da2fba4e1af34596b04d4","578da2fba4e1af34596b04d6"]}
            {"name": "size", "type": "string", "default": ""},  # 抽平的样本的大小, ""为不进行抽平
            {"name": "method", "type": "string", "default": ""}  # 样本合并的方式, ""为不进行样本合并
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.sub_sample = self.add_tool("meta.otu.sub_sample")
        group_table_path = os.path.join(self.work_dir, "group_table.xls")
        meta = Meta()
        meta._config = self.config  # 兼容不同mongo库版本
        self.group_table_path = meta.group_detail_to_table(self.option("group_detail"), group_table_path)

    def run_sort_samples(self):
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path,
            "method": self.option("method")
        })
        if self.option("size") != "":
            self.sort_samples.on("end", self.run_sub_sample)
        else:
            self.sort_samples.on("end", self.set_db)
        self.sort_samples.run()

    def run_sub_sample(self):
        self.sub_sample.set_options({
            "in_otu_table": self.sort_samples.option("out_otu_table"),
            "size": int(self.option("size"))
        })
        self.sub_sample.on("end", self.set_db)
        self.sub_sample.run()

    def set_db(self):
        out_otu = os.path.join(self.output_dir, "out_otu.xls")
        if self.option("size") == "":
            shutil.copy2(self.sort_samples.option("out_otu_table").prop["path"], out_otu)
        else:
            shutil.copy2(self.sub_sample.option("out_otu_table").prop["path"], out_otu)
        api_otu = self.api.sub_sample
        output_otu_id = api_otu.add_sg_otu(self.sheet.params, self.option("size"), self.option("input_otu_id"))
        if not os.path.isfile(out_otu):
            self.logger.error("找不到报告文件:{}".format(out_otu))
            self.set_error("找不到报告文件", code="12702901")
        self.logger.info("开始讲信息导入sg_otu_detail表和sg_otu_specimen表中")
        api_otu.add_sg_otu_detail(out_otu, self.option("input_otu_id"), output_otu_id)
        self.add_return_mongo_id("sg_otu", output_otu_id)
        self.end()

    def end(self):
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["out_otu.xls", 'meta.otu.otu_table', "经过otu统计后的表格"]
        ])
        super(OtuStaticWorkflow, self).end()

    def run(self):
        self.run_sort_samples()
        super(OtuStaticWorkflow, self).run()
