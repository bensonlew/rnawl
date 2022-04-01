# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""组间差异性两组比较检验分析"""

from biocluster.workflow import Workflow
import os
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.metaasv.common_function import link_dir, link_file
from mbio.packages.meta.save_params import save_params


class TwoGroupWorkflow(Workflow):
    """
    报告中调用组间差异性分析检验时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoGroupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "type", "type": "string", "default": "two.side"},
            # {"name": "update_info", "type": "string"},
            {"name": "test", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "group_name", "type": "string"},
            {"name": "coverage", "type": "float"},
            {"name": "params", "type": "string"},
            {"name": "category_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "signal_pair_file","type":"infile", "format":"meta.otu.otu_table"},
            {"name": "signal_pair_id", "type": "string", "default":""}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_group = self.add_tool("statistical.metastat")

    def run_two_group(self):
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
        #self.output_dir = self.two_group.output_dir
        self.two_group.run()

    def end(self):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两组比较结果目录", 0, "110139"],
            [r"./两组比较差异检验柱形图.pdf", "pdf", "两组比较中具有显著差异的丰度排行前15的物种差异检验柱形图", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, "110140"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量", 0, "110141"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "110142"]
        ])
        super(TwoGroupWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        link_dir(self.two_group.output_dir, self.output_dir)
        api_two_group = self.api.stat_test
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
        api_two_group.add_species_difference_check_detail(statfile=stat_path, cifiles=[ci_path], table_id=self.option('main_id'), level=self.option("level"), check_type='two_group', params=self.option("params"), category_name=self.option('category_name'), group_id=params["group_id"], from_otu_table=params["otu_id"], major=False, posthoc=None)
        api_two_group.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'))
        print bar_path
        api_two_group.add_species_difference_check_barplot(bar_path, self.option('main_id'))
        api_two_group.update_species_difference_check(self.option('main_id'), stat_path, ci_path, 'twogroup')
        self.add_name(stat_path)
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_species_difference_check")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "species_difference_two_group",
                "interaction": 1,
                "main_table": "sg_species_difference_check",
            })
            self.figsave.run()
        else:
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
