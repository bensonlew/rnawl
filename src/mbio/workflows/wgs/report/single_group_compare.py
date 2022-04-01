# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180612

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SingleGroupCompareWorkflow(Workflow):
    """
    交互分析：cnv的差异比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SingleGroupCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pop_final_vcf", "type": "string"},  # pop_final_vcf
            {"name": "group", "type": "string"},  # group：sample1,sample2
            {"name": "analysis_type", "type": "string"},  # 分析类型，snp、indel
            {"name": "location", "type": "string"},  # 基因组区域, chr1,1,10000
            {"name": "funtype", "type": "string"},  # 变异位点功效,逗号分隔，如High,Low,Moderate,Modifier
            {"name": "efftype", "type": "string"},  # 变异位点类型,逗号分隔，如5_prime_UTR_variant,3_prime_UTR_variant
            {"name": "maf", "type": "string"},  # group的基因型频率,用-分隔开前后值
            {"name": "dep", "type": "string"},  # group的平均测序深度,用-分隔开前后值
            {"name": "miss", "type": "string"},  # group的缺失率,用-分隔开前后值
            {"name": "len", "type": "string"},  # 变异长度，若为snp则固定为1，为indel时则是变异长度前面一个值
            {"name": "step", "type": "int", "default": 100000},  # 染色体分布图计算滑窗步长
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="14501001")
        if not self.option("group"):
            raise OptionError("请设置group分组", code="14501002")
        else:
            try:
                group_samples = self.option("group").split(":")[1]
            except:
                raise OptionError("%s分组传入格式不对，应为分组名称:样本名", variables=(self.option("group")), code="14501003")
                # raise OptionError("{}分组传入格式不对，应为分组名称:样本名".format(self.option("group")))

    def run_single_group_compare(self):
        options = {
            "pop_final_vcf": self.option("pop_final_vcf"),
            "group": self.option("group"),
            "analysis_type": self.option("analysis_type"),
            "location": self.option("location"),
            "funtype": self.option("funtype"),
            "efftype": self.option("efftype"),
            "maf": self.option("maf"),
            "dep": self.option("dep"),
            "miss": self.option("miss"),
            "len": self.option("len"),
            "step": self.option("step"),
            "task_id": self.option("task_id"),
            "main_id": self.option("main_id"),
            "project_type": self.option("project_type")
        }
        self.single_group_compare = self.add_tool("wgs.single_group_compare")
        self.single_group_compare.set_options(options)
        self.single_group_compare.on("end", self.set_output, "single_group_compare")
        self.single_group_compare.run()

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            new = os.path.join(self.output_dir, f)
            old = os.path.join(obj.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        if self.option("main_id"):
            self.set_db()

    def set_db(self):
        self.logger.info("保存结果到mongo")
        compare_id = self.option("main_id")
        compare_api = self.api.api("wgs.snp_indel_compare")
        if self.option("project_type"):
            compare_api._project_type = self.option("project_type")
        download_path = self._sheet.output.rstrip('/') + "/diff.variant"
        types = self.option("analysis_type")
        compare_api.update_download_path(compare_id, types, download_path)
        self.end()

    def run(self):
        self.run_single_group_compare()
        super(SingleGroupCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SingleGroupCompareWorkflow, self).end()
