# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from biocluster.workflow import Workflow
import os
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class FaprotaxGroupWorkflow(Workflow):
    """
    faprotax组间差异检验
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FaprotaxGroupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "faprotax_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "faprotax_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "group_method", "type": "string", "default": "two"},
            {"name": "test", "type": "string"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci_method", "type": "string", "default": "bootstrapr"},
            {"name": "coverage", "type": "float", "default": 0.95},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "methor", "type": "string"},
            {"name": "group_name", "type": "string"},
            {"name": "category_name", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "task_id", "type": "float"},
            {"name": "task_type", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_group = self.add_tool("statistical.metastat")
        self.multiple = self.add_tool("statistical.metastat")

    def run_two_group(self):
        if self.option("test") == "student":
            options = {
                "student_input": self.option("faprotax_table"),
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
                "mann_input": self.option("faprotax_table"),
                "mann_ci": self.option("ci"),
                "mann_group": self.option("group_file"),
                "mann_correction": self.option("correction"),
                "mann_type": self.option("type"),
                "test": self.option("test"),
                "mann_gname": self.option("group_name"),
                "mann_coverage": self.option("coverage")
            }
        else:
            options = {
                "welch_input": self.option("faprotax_table"),
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
        self.output_dir = self.two_group.output_dir
        self.two_group.run()

    def run_multiple(self):
        if self.option("test") == "anova":
            options = {
                "anova_input": self.option("faprotax_table"),
                "anova_group": self.option("group_file"),
                "anova_correction": self.option("correction"),
                "test": self.option("test"),
                "anova_gname": self.option("group_name"),
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage")
            }
        else:
            options = {
                "kru_H_input": self.option("faprotax_table"),
                "kru_H_group": self.option("group_file"),
                "kru_H_correction": self.option("correction"),
                "test": self.option("test"),
                "kru_H_gname": self.option("group_name"),
                "kru_H_methor": self.option("methor"),
                "kru_H_coverage": self.option("coverage")
            }
        self.multiple.set_options(options)
        self.multiple.on("end", self.set_db)
        self.output_dir = self.multiple.output_dir
        self.multiple.run()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
            if self.option("group_method") == "two":
                if os.path.exists(self.output_dir+"/多组功能组间差异检验柱形图.pdf"):
                    os.remove(self.output_dir+"/多组功能组间差异检验柱形图.pdf")
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("group_method") == "two":

            result_dir.add_relpath_rules([
                [".", "", "物种差异两组比较结果目录", 0, "110139"],
                ["./FAPROTAX两组差异检验柱形图.pdf", "pdf", "两组差异检验柱形图", 0, ""]
            ])
            result_dir.add_regexp_rules([
                [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, "110140"],
                [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量", 0, "110141"],
                [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "110142"]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "物种差异多组比较结果目录", 0, "110135"],
                ["./FAPROTAX多组差异检验柱形图.pdf", "pdf", "多组差异检验柱形图", 0, ""]
            ])
            result_dir.add_regexp_rules([
                [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, "110137"],
                # [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
                [r".*(-).*", "xls", "组间差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值", 0, "110136"],
                [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "110138"]
            ])
        super(FaprotaxGroupWorkflow, self).end()

    def set_db(self):
        """
        保存分析的结果表保存到mongo数据库中
        """
        if self.option("group_method") == "two":
            api_two_group = self.api.bugbase_group
            stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
            boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
            ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
            bar_path = self.two_group.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
            if not os.path.isfile(stat_path):
                self.logger.error("找不到报告文件:{}".format(stat_path))
                self.set_error("找不到报告文件", code="12704101")
            if not os.path.isfile(boxfile_path):
                self.logger.error("找不到报告文件:{}".format(boxfile_path))
                self.set_error("找不到报告文件", code="12704101")
            if not os.path.isfile(ci_path):
                self.logger.error("找不到报告文件:{}".format(ci_path))
                self.set_error("找不到报告文件", code="12704101")
            params = eval(self.option("params"))
            api_two_group.add_species_difference_check_detail(statfile=stat_path, cifiles=[ci_path],
                                                              table_id=self.option('main_id'),
                                                              major=False, posthoc=None,
                                                              correlation_key="faprotax_id",bar_name = "sg_faprotax_group_bar",
                                                              coll_name="sg_faprotax_group_detail", main_coll="sg_faprotax_group",task_type="two_group")
            api_two_group.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'),
                                                               correlation_key="faprotax_id", coll_name="sg_faprotax_group_box",
                                                               main_coll="sg_faprotax_group",task_type="two_group")
            #api_two_group.add_species_difference_check_barplot(bar_path, self.option('main_id'),
            #                                                   correlation_key="faprotax_id", coll_name="sg_faprotax_group_bar",
            #                                                   main_coll="sg_faprotax_group",task_type="two_group")
            api_two_group.update_species_difference_check(self.option('main_id'), stat_path, ci_path, 'sg_faprotax_group')
        else:
            api_multiple = self.api.bugbase_group
            stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
            stat_file = open(stat_path,"r")
            stat_data = stat_file.readlines()
            stat_file.close()
            stat_file2 = open(stat_path, "r")
            head = stat_file2.readline()
            if head.strip().split("\t")[0] == "name":
                pass
            else:
                stat_file3 = open(stat_path, "w")
                stat_file3.write("name" + "\t" + stat_data[0])
                for i in stat_data[1:]:
                    stat_file3.write(i)
                stat_file3.close()

            boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
            boxfile_file = open(boxfile_path, "r")
            boxfile_data = boxfile_file.readlines()
            boxfile_file.close()
            boxfile_file2 = open(boxfile_path, "r")
            head = boxfile_file2.readline()
            if head.strip().split("\t")[0] == "name":
                pass
            else:
                boxfile_file3 = open(boxfile_path, "w")
                boxfile_file3.write("name" + "\t" + boxfile_data[0])
                for i in boxfile_data[1:]:
                    boxfile_file3.write(i)
                boxfile_file3.close()

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
            api_multiple.add_species_difference_check_detail(statfile=stat_path, cifiles=cifiles,
                                                             table_id=self.option('main_id'),
                                                             major=False,
                                                             posthoc=self.option("methor"),
                                                             correlation_key="faprotax_id",
                                                             coll_name="sg_faprotax_group_detail",
                                                             main_coll="sg_faprotax_group",
                                                             bar_name="sg_faprotax_group_bar",
                                                             task_type="multiple_group")
            api_multiple.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'),
                                                              correlation_key="faprotax_id",
                                                              coll_name="sg_faprotax_group_box",
                                                              main_coll="sg_faprotax_group",
                                                              task_type="multiple_group")
            #api_multiple.add_species_difference_check_barplot(bar_path, self.option('main_id'),
            #                                                  correlation_key="faprotax_id",
            #                                                  coll_name="sg_faprotax_group_bar",
            #                                                  main_coll="sg_faprotax_group",
            #                                                  task_type="multiple_group")
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_faprotax_group")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "faprotax_group",
                "interaction": 1,
                "main_table": "sg_faprotax_group",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        if self.option("group_method") == "two":
            self.run_two_group()
        else:
            self.run_multiple()
        super(FaprotaxGroupWorkflow, self).run()
