# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180612

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class DoubleGroupCompareWorkflow(Workflow):
    """
    交互分析：cnv的差异比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DoubleGroupCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pop_final_vcf", "type": "infile", 'format': 'bsa.vcf'},  # pop_final_vcf
            {"name": "group1", "type": "string"},  # group：sample1,sample2
            {"name": "group2", "type": "string"},  # group：sample1,sample2
            {"name": "analysis_type", "type": "string"},  # 分析类型，snp、indel
            {"name": "location", "type": "string"},  # 基因组区域, chr1,1,10000
            {"name": "funtype", "type": "string"},  # 变异位点功效,逗号分隔，如High,Low,Moderate,Modifier
            {"name": "efftype", "type": "string"},  # 变异位点类型,逗号分隔，如5_prime_UTR_variant,3_prime_UTR_variant
            {"name": "maf1", "type": "string"},  # group1的基因型频率,用-分隔开前后值
            {"name": "maf2", "type": "string"},  # group2的基因型频率,用-分隔开前后值
            {"name": "dep1", "type": "string"},  # group1的平均测序深度,用-分隔开前后值
            {"name": "dep2", "type": "string"},  # group2的平均测序深度,用-分隔开前后值
            {"name": "miss1", "type": "string"},  # group1的缺失率,用-分隔开前后值
            {"name": "miss2", "type": "string"},  # group2的缺失率,用-分隔开前后值
            {"name": "len1", "type": "string"},  # 变异长度，若为snp则固定为1，为indel时则是变异长度前面一个值
            {"name": "len2", "type": "string"},  # 变异长度，若为snp则固定为1，为indel时则是变异长度后面一个值
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="14500401")
        if not self.option("group1"):
            raise OptionError("请设置group1分组", code="14500402")
        else:
            try:
                group1_samples = self.option("group1").split(":")[1]
            except:
                # raise OptionError("{}分组传入格式不对，应为分组名称:样本名".format(self.option("group1")))
                raise OptionError("%s分组传入格式不对，应为分组名称:样本名", variables=(self.option("group1")), code="14500403")
        if not self.option("group2"):
            raise OptionError("请设置group2分组", code="14500404")
        else:
            try:
                group2_samples = self.option("group2").split(":")[1]
            except:
                raise OptionError("%s分组格式不对，应为分组名称:样本名", variables=(self.option("group2")), code="14500405")
                # raise OptionError("{}分组格式不对，应为分组名称:样本名".format(self.option("group2")))

    def run_double_group_compare(self):
        options = {
            "pop_final_vcf": self.option("pop_final_vcf").prop['path'],
            "group1": self.option("group1"),
            "group2": self.option("group2"),
            "analysis_type": self.option("analysis_type"),
            "location": self.option("location"),
            "funtype": self.option("funtype"),
            "efftype": self.option("efftype"),
            "maf1": self.option("maf1"),
            "maf2": self.option("maf2"),
            "dep1": self.option("dep1"),
            "dep2": self.option("dep2"),
            "miss1": self.option("miss1"),
            "miss2": self.option("miss2"),
            "len1": self.option("len1"),
            "len2": self.option("len2"),
            "task_id": self.option("task_id"),
            "main_id": self.option("main_id"),
            "project_type": self.option("project_type")
        }
        self.double_group_compare = self.add_tool("wgs.double_group_compare")
        self.double_group_compare.set_options(options)
        self.double_group_compare.on("end", self.set_output, "double_group_compare")
        self.double_group_compare.run()

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
        self.run_double_group_compare()
        super(DoubleGroupCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DoubleGroupCompareWorkflow, self).end()
