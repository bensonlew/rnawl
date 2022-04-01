# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import link_dir, link_file
from biocluster.workflow import Workflow
import os


class MultipleGroupWorkflow(Workflow):
    """
    metaasv 物种差异分析多组比较
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultipleGroupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string","default": ""},
            {"name": "update_info", "type": "string", "default": ""},
            {"name": "test", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": ""},
            {"name": "params", "type": "string","default": ""},
            {"name": "group_name", "type": "string", "default": ""},
            {"name": "methor", "type": "string"},
            {"name": "coverage", "type": "float"},
            {"name": "category_name", "type": "string", "default": ""},
            {"name": "update_info", "type": "string","default": ""},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.multiple = self.add_tool("statistical.metastat")

    def run_multiple(self):
        """
        运行多组比较 metastat
        """
        if self.option("test") == "anova":
            options = {
                "anova_input": self.option("otu_file"),
                "anova_group": self.option("group_file"),
                "anova_correction": self.option("correction"),
                "test": self.option("test"),
                "anova_gname": self.option("group_name"),
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage")
            }
        else:
            options = {
                "kru_H_input": self.option("otu_file"),
                "kru_H_group": self.option("group_file"),
                "kru_H_correction": self.option("correction"),
                "test": self.option("test"),
                "kru_H_gname": self.option("group_name"),
                "kru_H_methor": self.option("methor"),
                "kru_H_coverage": self.option("coverage")
            }
        self.multiple.set_options(options)
        self.multiple.on("end", self.set_db)
        self.multiple.run()

    def end(self):
        """
        上传结果文件和结束
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异多组比较结果目录", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, ""],
            [r".*(-).*", "xls", "组间差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值", 0, ""],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, ""]
        ])
        super(MultipleGroupWorkflow, self).end()

    def set_db(self):
        """
        连接结果文件和导入MongoDB
        """
        link_dir(self.multiple.output_dir, self.output_dir)
        api_multiple = self.api.api("metaasv.stat_test")
        stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
        bar_path = self.multiple.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
        if not os.path.isfile(stat_path):
            self.logger.error("找不到报告文件:{}".format(stat_path))
            self.set_error("找不到报告文件", code="")
        if not os.path.isfile(boxfile_path):
            self.logger.error("找不到报告文件:{}".format(boxfile_path))
            self.set_error("找不到报告文件", code="")
        cifiles = []
        for r, d, f in os.walk(self.output_dir):
            for i in f:
                if self.option("methor") in i:
                    ci_path = r + '/' + i
                    if not os.path.isfile(ci_path):
                        self.logger.error("找不到报告文件:{}".format(ci_path))
                        self.set_error("找不到报告文件")
                    cifiles.append(ci_path)
        api_multiple.add_species_difference_check_detail(statfile=stat_path, cifiles=cifiles, table_id=self.option('main_id'), level=self.option("level"), major=False, posthoc=self.option("methor"), correlation_key="multiple_id", coll_name="multiple_group_detail", main_coll="multiple_group")
        api_multiple.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'), correlation_key="multiple_id", coll_name="multiple_group_box", main_coll="multiple_group")
        api_multiple.add_species_difference_check_barplot(bar_path, self.option('main_id'), correlation_key="multiple_id", coll_name="multiple_group_bar", main_coll="multiple_group")
        self.add_name(stat_path)
        self.add_name(boxfile_path)
        self.end()

    def add_name(self,file):
        os.rename(file,file+".tmp")
        with open(file+".tmp","r") as f, open(file,"w") as t:
            data = f.readlines()
            if self.option("test") == "anova":
                statis = "F_statistic"
            else:
                statis = "H_statistic"
            t.write(""+"\t"+data[0].replace("statistic",statis))
            for i in data[1:]:
                t.write(i)
        os.remove(file+".tmp")

    def report_option(self):
        """
        打印参数
        :return:
        """
        self.logger.info('option otu_file : ' + self.option('otu_file').prop['path'])
        self.logger.info('option group_file : ' + self.option('group_file').prop['path'])
        # self.logger.info('option group_detail : ' + self.option('group_detail'))
        self.logger.info('option update_info : ' + self.option('update_info'))
        self.logger.info('option test : ' + self.option('test'))
        level = self.option('level')
        level = str(level)
        self.logger.info('option level : ' + level)
        self.logger.info('option correction : ' + self.option('correction'))
        self.logger.info('option params : ' + self.option('params'))
        self.logger.info('option group_name : ' + self.option('group_name'))
        self.logger.info('option methor : ' + self.option('methor'))
        coverage = self.option('coverage')
        coverage = str(coverage)
        self.logger.info('option coverage : ' + coverage)
        self.logger.info('option category_name : ' + self.option('category_name'))
        self.logger.info('option main_id : ' + self.option('main_id'))

    def run(self):
        self.run_multiple()
        self.report_option()
        super(MultipleGroupWorkflow, self).run()
