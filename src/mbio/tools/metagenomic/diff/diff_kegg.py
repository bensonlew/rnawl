# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2018/12/10'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.tools.statistical.metastat import MetastatAgent,MetastatTool
from mbio.packages.statistical.metastat import diff_kegg_test

class DiffKeggAgent(Agent):
    """
    version 1.0
    """

    def __init__(self, parent):
        super(DiffKeggAgent, self).__init__(parent)
        options = [
            {"name": "test", "type": "string"},  # 两组比较检验方法
            {"name": "kegg_path", "type": "string"}, # kegg注释结果文件路径
            # {"name": "xml_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "img_path", "type": "string"}, # 图片生成的结果路径
            {"name": "test_result", "type": "outfile", "format": "sequence.profile_table"}, # 差异检验结果
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option('test') or self.option('test') not in ["mann", "signal", "welch", "student"]:
            raise OptionError("所输入的检验名称不对", code="32800101")
        if not self.option('kegg_path'):
            raise OptionError("需要设置kegg结果路径", code="32800102")
        if not self.option("group").is_set:
            raise OptionError("需要设置分组表", code="32800103")
        if self.option("correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                                                          "none"]:
            raise OptionError("多重检验校正的方法不被支持", code="32800104")
        if self.option("ci") <= 0 or self.option("ci") >= 1:
            raise OptionError("所输入的显著水平不在范围值内", code="32800105")
        if self.option("type") not in ["two.side", "greater", "less"]:
            raise OptionError("所输入的类型不在范围内", code="32800106")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'


class DiffKeggTool(Tool):
    def __init__(self, config):
        super(DiffKeggTool, self).__init__(config)
        self.r_path = '/program/R-3.3.1/bin/Rscript'
        self.pic_path = self.config.PACKAGE_DIR + "/annotation/mg_annotation/kegg_pathway_img.py"
        self.python_path = "/program/Python/bin/python"
        self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"

    def run_pic(self):
        """
        description
        :return:
        """
        self.logger.info("start output kegg img")
        pathway_file = self.option("kegg_path") + "/kegg_pathway_eachmap.xls"
        enzyme_file = self.option("kegg_path") + "/kegg_enzyme_profile.xls"
        cmd = "{} {} -o {} -p {} -ko {} -KO {} -png_file {} -html {} -enzyme_file {}".format(
            self.python_path, self.pic_path, self.output_dir, pathway_file, "Pathway",
            "KO_list", "True", self.html_path, enzyme_file)
        command = self.add_command("output_kegg_pathway_img", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("output_kegg_pathway_img succeed")
        else:
            self.set_error("output kegg pathway img failed", code="32800101")

    def run_diff(self):
        self.logger.info("start two group test")
        protein_profile_file = self.option("kegg_path") + "/kegg_KO_profile.xls"
        enzyme_profile_file = self.option("kegg_path") + "/kegg_enzyme_profile.xls"
        filter_file = self.output_dir + '/filter.xls'
        output_file = self.output_dir + '/test_result.xls'
        diff_kegg_test(protein_profile_file, enzyme_profile_file, filter_file, self.option("group").path, output_file,
                       choose_test=self.option("test"), ci=self.option("ci"), test_type=self.option("type"),
                       mul_test=self.option("correction"))
        cmd = self.r_path + " run_%s_test.r" % self.option("test")
        lines = []
        with open(filter_file, 'r') as r:
            for l in r:
                l.strip() and lines.append(l)
        if not lines:
            self.set_error("无有效的KEGG注释用于差异分析")

        self.logger.info("开始运行差异检验")
        command = self.add_command("diff_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("差异检验运行完成！")
        else:
            self.set_error("差异检验运行出错！", code="32800102")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option("img_path", self.output_dir + "/pathway_img")
        self.option("test_result", self.output_dir + "/test_result.xls")
        self.logger.debug(self.option("img_path"))
        self.logger.debug(self.option("test_result").path)
        self.logger.debug(type(self.option("img_path")))

    def run(self):
        super(DiffKeggTool, self).run()
        self.run_pic()
        self.run_diff()
        self.set_output()
        self.end()