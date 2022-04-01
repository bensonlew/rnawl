# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import link_dir, link_file
from biocluster.workflow import Workflow
import os


class TwoGroupWorkflow(Workflow):
    """
    metaasv 物种差异分析两组比较
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoGroupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "group_name", "type": "string", "default": ""},
            {"name": "coverage", "type": "float"},
            {"name": "params", "type": "string", "default": ""},
            {"name": "category_name", "type": "string", "default": ""},
            {"name": "update_info", "type": "string", "default": ""},
            {"name": "main_id", "type": "string"},
            {"name": "signal_pair_file","type":"infile", "format":"meta.otu.otu_table"},
            {"name": "signal_pair_id", "type": "string", "default":""}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.two_group = self.add_tool("statistical.metastat")

    def run_two_group(self):
        """
        运行两组比较 metastat
        """
        if self.option("test") == "student":
            options = {
                "student_input": self.option("otu_file"),
                "student_group": self.option("group_file"),
                "student_ci": self.option("ci"),
                "student_correction": self.option("correction"),
                "student_type": self.option("type"),
                "test": self.option("test"),
                "student_gname": self.option("group_name"),
                "student_coverage": self.option("coverage")
            }
        elif self.option("test") == "mann":
            options = {
                "mann_input": self.option("otu_file"),
                "mann_ci": self.option("ci"),
                "mann_group": self.option("group_file"),
                "mann_correction": self.option("correction"),
                "mann_type": self.option("type"),
                "test": self.option("test"),
                "mann_gname": self.option("group_name"),
                "mann_coverage": self.option("coverage")
            }
        elif self.option("test") == "signal":
            options = {
                "signal_input": self.option("otu_file"),
                "signal_ci": self.option("ci"),
                "signal_group": self.option("group_file"),
                "signal_correction": self.option("correction"),
                "signal_type": self.option("type"),
                "test": self.option("test"),
                "signal_gname": self.option("group_name"),
                "signal_coverage": self.option("coverage"),
                "signal_pair_file" : self.option("signal_pair_file")
            }
        else:
            options = {
                "welch_input": self.option("otu_file"),
                "welch_ci": self.option("ci"),
                "welch_group": self.option("group_file"),
                "welch_correction": self.option("correction"),
                "welch_type": self.option("type"),
                "test": self.option("test"),
                "welch_gname": self.option("group_name"),
                "welch_coverage": self.option("coverage")
            }
        self.two_group.set_options(options)
        self.two_group.on("end", self.set_db)
        self.two_group.run()

    def end(self):
        """
        结束和上传结果文件
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两组比较结果目录", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, ""],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量", 0, ""],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, ""]
        ])
        super(TwoGroupWorkflow, self).end()

    def set_db(self):
        """
        连接文件目录和导入MongoDB数据
        """
        link_dir(self.two_group.output_dir, self.output_dir)
        api_two_group = self.api.api("metaasv.stat_test")
        stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
        ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
        bar_path = self.two_group.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
        if not os.path.isfile(stat_path):
            self.logger.error("找不到报告文件:{}".format(stat_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(boxfile_path):
            self.logger.error("找不到报告文件:{}".format(boxfile_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(ci_path):
            self.logger.error("找不到报告文件:{}".format(ci_path))
            self.set_error("找不到报告文件")
        # params = eval(self.option("params"))two_group_bar
        api_two_group.add_species_difference_check_detail(statfile=stat_path, cifiles=[ci_path], table_id=self.option('main_id'), level=self.option("level"), major=False, posthoc=None, correlation_key="compare_id", coll_name="two_group_detail", main_coll="two_group")
        api_two_group.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'), correlation_key="compare_id", coll_name="two_group_box", main_coll="two_group")
        api_two_group.add_species_difference_check_barplot(bar_path, self.option('main_id'), correlation_key="compare_id", coll_name="two_group_bar", main_coll="two_group")
        api_two_group.update_species_difference_check(self.option('main_id'), stat_path, ci_path, 'twogroup')
        self.add_name(stat_path)
        self.end()

    def add_name(self,file):
        os.rename(file,file+".tmp")
        count = 0
        for index, line in enumerate(open(self.option("group_file").prop["path"], 'r')):
            count += 1
        with open(file+".tmp","r") as f, open(file,"w") as t:
            data = f.readlines()
            if self.option("test") in ["mann","signal"]:
                t.write(""+"\t"+data[0].strip().replace("statistic", "W_statistic")+"\t"+"Z_value"+"\n")
                for i in data[1:]:
                    if i.strip().split("\t")[5] != "NA":
                        W = float(i.strip().split("\t")[5])
                        z_score = ((abs(W-float((count-1)*(count-1))/4))-0.5)/((float((count-1)*count*(2*count+1))/24)**0.5)
                        t.write(i.strip()+"\t"+str(z_score)+"\n")
                    else:
                        t.write(i.strip() + "\t" + "NA" + "\n")
            else:
                t.write(data[0].replace("statistic", "T_statistic"))
                for i in data[1:]:
                    t.write(i)
        os.remove(file + ".tmp")

    def run(self):
        self.run_two_group()
        super(TwoGroupWorkflow, self).run()
