# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class DoubleGroupCompareAgent(Agent):
    """
    snp/indel 样本组与样本组间进行比较
    """
    def __init__(self, parent):
        super(DoubleGroupCompareAgent, self).__init__(parent)
        options = [
            {"name": "pop_final_vcf", "type": "string"},  # pop_final_vcf
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

    def check_options(self):
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="34502101")
        if not self.option("group1"):
            raise OptionError("请设置group1分组", code="34502102")
        else:
            try:
                group1_samples = self.option("group1").split(":")[1]
            except:
                raise OptionError("%s分组传入格式不对，应为分组名称:样本名",variables=(self.option("group1")), code="34502103")
        if not self.option("group2"):
            raise OptionError("请设置group2分组", code="34502104")
        else:
            try:
                group2_samples = self.option("group2").split(":")[1]
            except:
                raise OptionError("%s分组格式不对，应为分组名称:样本名",variables=(self.option("group2")), code="34502105")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(DoubleGroupCompareAgent, self).end()


class DoubleGroupCompareTool(Tool):
    def __init__(self, config):
        super(DoubleGroupCompareTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.double_group = self.config.PACKAGE_DIR + "/wgs/double_group_vcf.pl"

    def run_double_group_vcf(self):
        """
        double_group_vcf.pl
        """
        self.group1_name = self.option("group1").split(":")[0]
        self.group1_samples = self.option("group1").split(":")[1]
        self.group2_name = self.option("group2").split(":")[0]
        self.group2_samples = self.option("group2").split(":")[1]
        cmd = "{} {} -vcf {} ".format(self.perl_path, self.double_group, self.option("pop_final_vcf"))
        cmd += "-gro1 {} -gro2 {} ".format(self.group1_samples, self.group2_samples)
        cmd += "-g1name {} -g2name {} -out {}".format(self.group1_name, self.group2_name, self.work_dir + "/diff")
        if self.option("location"):
            cmd += " -pos {}".format(self.option("location"))
        if self.option("funtype"):
            cmd += " -funtype {}".format(self.option("funtype"))
        if self.option("efftype"):
            cmd += " -efftype {}".format(self.option("efftype"))
        if self.option("maf1"):
            cmd += " -maf1 {}".format(self.option("maf1"))
        if self.option("maf2"):
            cmd += " -maf2 {}".format(self.option("maf2"))
        if self.option("dep1"):
            cmd += " -dep1 {}".format(self.option("dep1"))
        if self.option("dep2"):
            cmd += " -dep2 {}".format(self.option("dep2"))
        if self.option("miss1"):
            cmd += " -miss1 {}".format(self.option("miss1"))
        if self.option("miss2"):
            cmd += " -miss2 {}".format(self.option("miss2"))
        if self.option("len1"):
            cmd += " -len1 '\{}'".format(self.option("len1"))
        if self.option("len2"):
            cmd += " -len2 '\{}'".format(self.option("len2"))
        command = self.add_command("double_group_vcf", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("double_group_vcf.pl运行完成")
        else:
            self.set_error("double_group_vcf.pl运行失败", code="34502101")

    def set_output(self):
        os.link(self.work_dir + "/diff/Ann.stat", self.output_dir + "/Ann.stat")
        os.link(self.work_dir + "/diff/Eff.stat", self.output_dir + "/Eff.stat")
        os.link(self.work_dir + "/diff/win.stat", self.output_dir + "/win.stat")
        if self.option("location"):
            os.link(self.work_dir + "/diff/pos.variant", self.output_dir + "/diff.variant")
        else:
            os.link(self.work_dir + "/diff/pop.variant", self.output_dir + "/diff.variant")

    def set_db(self):
        self.logger.info("保存结果到mongo")
        task_id = self.option("task_id")
        compare_id = self.option("main_id")
        compare_api = self.api.api("wgs.snp_indel_compare")
        if self.option("project_type"):
            compare_api._project_type = self.option("project_type")
        anno_file = self.output_dir + "/Ann.stat"
        eff_stat = self.output_dir + "/Eff.stat"
        variant_file = self.output_dir + "/diff.variant"
        win_stat = self.output_dir + "/win.stat"
        if self.option("analysis_type") == "snp":
            compare_api.add_sg_snp_indel_compare_stat(task_id, compare_id, anno_file, "snp")
            compare_api.add_sg_snp_indel_compare_eff_stat(task_id, compare_id, eff_stat, "snp")
            compare_api.add_sg_snp_indel_compare_detail(compare_id, variant_file, "snp", "two_group")
            compare_api.add_sg_manhattan(compare_id, variant_file, "snp")
            compare_api.add_sg_distribution(compare_id, win_stat, 100000, types="snp")
        else:
            compare_api.add_sg_snp_indel_compare_stat(task_id, compare_id, anno_file, "indel")
            compare_api.add_sg_snp_indel_compare_eff_stat(task_id, compare_id, eff_stat, "indel")
            compare_api.add_sg_snp_indel_compare_detail(compare_id, variant_file, "indel", "two_group")
            compare_api.add_sg_manhattan(compare_id, variant_file, "indel")
            compare_api.add_sg_distribution(compare_id, win_stat, 100000, types="indel")

    def run(self):
        super(DoubleGroupCompareTool, self).run()
        self.run_double_group_vcf()
        self.set_output()
        if self.option("main_id"):
            self.set_db()
        self.end()
