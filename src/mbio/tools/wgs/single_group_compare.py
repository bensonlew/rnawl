# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last modify 20180416

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SingleGroupCompareAgent(Agent):
    """
    snp/indel比较分析接口中用于组内比较分析, 所有参数多选的时候，都以英文逗号分隔
    """
    def __init__(self, parent):
        super(SingleGroupCompareAgent, self).__init__(parent)
        options = [
            {"name": "pop_final_vcf", "type": "infile", 'format': 'bsa.vcf'},  # pop_final_vcf
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
        # self.step.add_steps('cnvdiff')
        # self.on('start', self.step_start)
        # self.on('end', self.step_end)

    # def step_start(self):
    #     self.step.cnvdiff.start()
    #     self.step.update()
    #
    # def step_end(self):
    #     self.step.cnvdiff.finish()
    #     self.step.update()

    def check_options(self):
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="34505601")
        if not self.option("group"):
            raise OptionError("请设置group分组", code="34505602")
        else:
            try:
                group_samples = self.option("group").split(":")[1]
            except:
                raise OptionError("%s分组传入格式不对，应为分组名称:样本名", variables=(self.option("group")), code="34505603")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        super(SingleGroupCompareAgent, self).end()


class SingleGroupCompareTool(Tool):
    def __init__(self, config):
        super(SingleGroupCompareTool, self).__init__(config)
        self.script_path = self.config.PACKAGE_DIR + "/wgs/single_group_vcf.pl"
        self.perl_path = 'miniconda2/bin/perl '
        self.python = "miniconda2/bin/python"
        self.distribution_graph_path = self.config.PACKAGE_DIR + "/wgs/distribution_graph.py"

    def single_group_compare(self):
        group_samples = self.option("group").split(":")[1]
        cmd = "{} {} -gro {} -vcf {} -out {}".format(self.perl_path, self.script_path, group_samples,
                                                     self.option("pop_final_vcf").prop['path'], self.work_dir + "/diff")
        if self.option("location"):
            cmd += " -pos {}".format(self.option("location"))
        if self.option("funtype"):
            cmd += " -funtype {}".format(self.option("funtype"))
        if self.option("efftype"):
            cmd += " -efftype {}".format(self.option("efftype"))
        if self.option("maf"):
            cmd += " -maf {}".format(self.option("maf"))
        if self.option("dep"):
            cmd += " -dep {}".format(self.option("dep"))
        if self.option("miss"):
            cmd += " -miss {}".format(self.option("miss"))
        if self.option("len"):
            cmd += " -len '\{}'".format(self.option("len"))
        self.logger.info("开始进行single_group_compare")
        command = self.add_command("single_group_compare", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("single_group_compare完成！")
        else:
            self.set_error("single_group_compare出错！", code="34505601")
            self.set_error("single_group_compare出错！", code="34505607")

    def set_output(self):
        os.link(self.work_dir + "/diff/Ann.stat", self.output_dir + "/Ann.stat")
        os.link(self.work_dir + "/diff/Eff.stat", self.output_dir + "/Eff.stat")
        os.link(self.work_dir + "/diff/win.stat", self.output_dir + "/win.stat")
        os.link(self.work_dir + "/diff/diff.matrix", self.output_dir + "/diff.matrix")
        if self.option("location"):
            os.link(self.work_dir + "/diff/pos.variant", self.output_dir + "/diff.variant")
        else:
            os.link(self.work_dir + "/diff/pop.variant", self.output_dir + "/diff.variant")

    def distribution_graph(self):
        """
        将win.stat进行滑窗，用于画染色体分布图
        """
        win_stat = self.work_dir + "/diff/win.stat"
        step_win_stat = self.work_dir + "/win.stat.xls"
        cmd = "{} {} -step {} -i {} -o {}".format(self.python, self.distribution_graph_path, self.option("step"), win_stat, step_win_stat)
        self.logger.info("开始进行染色体分布图滑窗")
        command = self.add_command("distribution_graph", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("distribution_graph完成！")
        else:
            self.set_error("distribution_graph出错！", code="34505602")
            self.set_error("distribution_graph出错！", code="34505608")

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
        # win_stat = self.output_dir + "/win.stat"
        win_stat = self.work_dir + "/win.stat.xls"
        heatmap_file = self.output_dir + "/diff.matrix"
        if self.option("analysis_type") == "snp":
            compare_api.add_sg_snp_indel_compare_stat(task_id, compare_id, anno_file, "snp")
            compare_api.add_sg_snp_indel_compare_eff_stat(task_id, compare_id, eff_stat, "snp")
            compare_api.add_sg_snp_indel_compare_detail(compare_id, variant_file, "snp", "one_group")
            # compare_api.add_sg_distribution(compare_id, win_stat, "snp")
            compare_api.add_sg_distribution_new(compare_id, win_stat, "snp")
            compare_api.add_sg_heatmap(compare_id, heatmap_file, "snp")
        else:
            compare_api.add_sg_snp_indel_compare_stat(task_id, compare_id, anno_file, "indel")
            compare_api.add_sg_snp_indel_compare_eff_stat(task_id, compare_id, eff_stat, "indel")
            compare_api.add_sg_snp_indel_compare_detail(compare_id, variant_file, "indel", "one_group")
            # compare_api.add_sg_distribution(compare_id, win_stat, "indel")
            compare_api.add_sg_distribution_new(compare_id, win_stat, "indel")
            compare_api.add_sg_heatmap(compare_id, heatmap_file, "indel")

    def run(self):
        super(SingleGroupCompareTool, self).run()
        self.single_group_compare()
        self.distribution_graph()
        self.set_output()
        if self.option("main_id"):
            self.set_db()
        self.end()
