# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.workflow import Workflow
import os
from mbio.packages.metaasv.common_function import link_dir
from mbio.packages.alpha_diversity.group_file_split import group_file_spilt


class AlphaCompareWorkflow(Workflow):
    """
    metaasv alpha多样性指数组间差异比较
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(AlphaCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "est_table", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的多样性指数表
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入的group表
            {"name": "est_id", "type": "string"},##alpha多样性主表的ID
            {"name": "asv_id", "type": "string"},##ASV主表的ID 用于框架中的data传递
            {"name": "group_id", "type": "string"},
            {"name": "main_id", "type": "string"},##主表的ID
            {"name": "group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "analysis", "type": "string"},##分析是两组还是多组
            {"name": "test", "type": "string"}, ##检验方法
            {"name": "group_name", "type": "string"},##分组名称
            {"name": "methor", "type": "string"}, ##post_hoc检验方法
            {"name": "coverage", "type": "float"}, ##post_hoc检验方法的阈值，默认0.95
            {"name": "correction", "type": "string", "default": ""},##多重检验方法
            {"name": "category_name", "type": "string"},##分组名称
            {"name": "ci", "type": "float", "default": 0.05},##两组检验ci值
            {"name": "type", "type": "string", "default": "two.side"}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.group_file_dir = self.work_dir + '/two_group_output'

    def run_multiple(self):
        """
        运行多组 metastat
        :return:
        """
        self.multiple = self.add_tool("statistical.metastat")
        if self.option("test") == "anova":
            options = {
                "anova_input": self.option("est_table"),
                "anova_group": self.option("group_table"),
                "anova_correction": self.option("correction"),
                "test": self.option("test"),
                "anova_gname": self.option("group_name"),
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage"),
                "need_percents": "false",
                "percent_abund": "false"
            }
        else:
            options = {
                "kru_H_input": self.option("est_table"),
                "kru_H_group": self.option("group_table"),
                "kru_H_correction": self.option("correction"),
                "test": self.option("test"),
                "kru_H_gname": self.option("group_name"),
                "kru_H_methor": self.option("methor"),
                "kru_H_coverage": self.option("coverage"),
                "need_percents": "false",
                "percent_abund": "false"
            }
        self.multiple.set_options(options)
        self.multiple.on("end", self.set_db)
        self.multiple.run()

    def run_two_group(self):
        """
        运行两组 metastat
        :return:
        """
        group_name = group_file_spilt(self.option('group_table').prop['path'], self.group_file_dir)
        name_list = []
        for g in group_name:
            if g[0] > g[1]:
                gg = g[1]+'|'+g[0]
                name_list.append(gg)
            else:
                gg = g[0]+'|'+g[1]
                name_list.append(gg)
        self.group_name = ",".join(name_list)
        self.logger.info(self.group_name)
        self.two_group = self.add_tool("statistical.metastat")
        options = {
            "est_input": self.option("est_table"),
            "est_group": self.group_file_dir,
            "est_correction": self.option("correction"),##多重检验方法
            "est_type": self.option("type"), ##双尾检验
            "test": "estimator",
            "est_test_method": self.option("test"),##检验方法
        }
        self.two_group.set_options(options)
        self.two_group.on("end", self.set_db)
        self.two_group.run()

    def set_db(self):
        """
        链接结果文件和导入MongoDB
        """
        if self.option("analysis") in ["multiple"]:
            link_dir(self.multiple.output_dir, self.output_dir)
            api_multiple = self.api.api("metaasv.alpha_compare")
            stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
            boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
            # bar_path = self.multiple.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
            cifiles = []
            for r, d, f in os.walk(self.output_dir):
                for i in f:
                    if self.option("methor") in i:
                        ci_path = r + '/' + i
                        if not os.path.isfile(ci_path):
                            self.logger.error("找不到报告文件:{}".format(ci_path))
                            self.set_error("找不到报告文件")
                        cifiles.append(ci_path)
            api_multiple.add_species_difference_check_detail(statfile=stat_path, cifiles=cifiles, table_id=self.option('main_id'), major=False, posthoc=self.option("methor"), correlation_key="diff_id", coll_name="alpha_diversity_diff_table", main_coll="alpha_diversity_diff", type="multiple")
            api_multiple.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'), correlation_key="diff_id", coll_name="alpha_diversity_diff_box", main_coll="alpha_diversity_diff")
            api_multiple.add_species_difference_check_barplot(stat_path, self.option('main_id'), correlation_key="diff_id", coll_name="alpha_diversity_diff_bar", main_coll="alpha_diversity_diff")
            self.end()
        elif self.option("analysis") in ["two_group"]:
            link_dir(self.two_group.output_dir, self.output_dir)
            api_two_group = self.api.api("metaasv.alpha_compare")
            api_two_group.add_est_detail(statfile=self.output_dir,table_id=self.option('main_id'), correlation_key="diff_id", coll_name="alpha_diversity_diff_table", main_coll="alpha_diversity_diff", type="two_group")
            self.end()

    def run(self):
        if self.option("analysis") in ["multiple"]:
            self.run_multiple()
        elif self.option("analysis") in ["two_group"]:
            self.run_two_group()
        super(AlphaCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "alpha多样性指数检验结果目录", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*\.xls", "", "alpha多样性指数T检验结果表", 0, ""]
        ])
        super(AlphaCompareWorkflow, self).end()
