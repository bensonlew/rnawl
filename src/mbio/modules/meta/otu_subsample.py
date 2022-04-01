# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class OtuSubsampleModule(Module):

    def __init__(self, work_id):
        super(OtuSubsampleModule, self).__init__(work_id)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "size", "type": "string", "default": "min"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {'name': 'otu_taxon_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出的otu表文件
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_dir'}, # 输出的otu_taxon_dir(absolute)文件夹
        ]
        self.add_option(options)
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.subsample = self.add_tool("meta.otu.sub_sample")
        self.tax_stat = self.add_tool("meta.otu.otu_taxon_stat")

    def run_sort_samples(self):
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.option("group")
        })
        self.sort_samples.on("end", self.run_subsample)
        self.sort_samples.run()

    def run_subsample(self):
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
        num_lines = sum(1 for line in open(self.subsample.option("out_otu_table").prop["path"]))
        if num_lines < 2:
            self.logger.error("经过抽平之后的OTU表是空的，可能是因为进行物种筛选之后导致某些样本的序列数为0，然后按该样本的序列数进行了抽平！")
            self.set_error("经过OTU过滤之后的OTU表是空的，请重新填写筛选的条件！", code="12703001")
        final_file = self.subsample.option("out_otu_table")
        self.table2db = final_file.prop['path']
        self.tax_stat.set_options({"sub_otu_table": final_file})
        self.tax_stat.on("end", self.set_output)
        self.tax_stat.run()

    def set_output(self):
        lst = os.listdir(self.tax_stat.output_dir)
        for l in lst:
            path = os.path.join(self.tax_stat.output_dir, l)
            final_path = os.path.join(self.output_dir, l)
            if os.path.exists(final_path):
                if os.path.isdir(final_path):
                    shutil.rmtree(final_path)
                else:
                    os.remove(final_path)
            if os.path.isdir(path):
                shutil.copytree(path, final_path)
            else:
                shutil.copy(path, final_path)
        self.option("otu_taxon_table").set_path(self.output_dir + "/otu_summary.xls")
        self.option("otu_taxon_dir").set_path(self.output_dir + "/tax_summary_a")
        self.end()

    def run(self):
        super(OtuSubsampleModule, self).run()
        self.run_sort_samples()

    def end(self):
        super(OtuSubsampleModule, self).end()

