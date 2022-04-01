# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180503

import os
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class GeneVariantSitesWorkflow(Workflow):
    """
    交互分析：基因详情页-变异位点
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneVariantSitesWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pop_final_vcf", "type": "string"},  # pop_final_vcf
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "is_same", "type": "string", "default": "true"},  # 页面样本间拷贝数变异是否相同
            {"name": "group1", "type": "string"},  # group：sample1,sample2
            {"name": "group2", "type": "string"},  # group：sample1,sample2
            {"name": "variant_type", "type": "string", "default": "all"},  # 变异类型，snp、indel, all
            {"name": "location", "type": "string"},  # 基因组区域, chr1,1,10000
            {"name": "funtype", "type": "string", "default": "all"},  # 变异位点功效,逗号分隔，如High,Low,Moderate,Modifier
            {"name": "efftype", "type": "string", "default": "all"},  # 变异位点类型,逗号分隔，如5_prime_UTR_variant,3_prime_UTR_variant
            {"name": "maf1", "type": "string", "default": "0,1.2"},  # group1的基因型频率,用-分隔开前后值
            {"name": "maf2", "type": "string", "default": "0,1.2"},  # group2的基因型频率,用-分隔开前后值
            {"name": "dep1", "type": "string", "default": "-1,100000000"},  # 2,5
            {"name": "dep2", "type": "string", "default": "-1,100000000"},  # 2,6
            {"name": "miss1", "type": "string", "default": "0"},  # group1的缺失率,用-分隔开前后值
            {"name": "miss2", "type": "string", "default": "0"},  # group2的缺失率,用-分隔开前后值
            {"name": "len1", "type": "string", "default": "-1,100000000"},  # 变异长度，若为snp则固定为1，为indel时则是变异长度前面一个值
            {"name": "len2", "type": "string", "default": "-1,100000000"},  # 变异长度，若为snp则固定为1，为indel时则是变异长度后面一个值
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="14500601")

    def run_sample_compare(self):
        """
        样本比较分析
        """
        len1 = self.option("len1")
        options = {
            "vcf_file": self.option("pop_final_vcf"),
            "sample1": self.option("sample1"),
            "sample2": self.option("sample2"),
            "is_same": self.option("is_same"),
            "funtype": self.option("funtype"),
            "efftype": self.option("efftype"),
            "len1": self.option("len1"),
            "len2": self.option("len2"),
            "dep1": self.option("dep1"),
            "dep2": self.option("dep2"),
            "location": self.option("location")
        }
        self.sample_compare = self.add_module("wgs.sample_compare")
        self.sample_compare.set_options(options)
        self.sample_compare.on("end", self.set_output)
        self.sample_compare.run()
