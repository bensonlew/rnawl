# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import shutil
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import gevent
import pandas as pd
from mbio.packages.metaasv.common_function import link_dir,link_file


class CompositionAnalysisModule(Module):
    """
    metaasv 群落barpie和heatmap
    """
    def __init__(self, work_id):
        super(CompositionAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string", "default": "barpie,heatmap"},###分析类型
            {"name": "abundtable", "type": "infile", "format": "meta.otu.otu_table"},  # 物种/功能丰度表格
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},###分组表
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "method", "type": "string", "default": "average"},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": "average"},  # 样本层次聚类方式，默认不聚类
            {"name": "species_number", "type": "string", "default": "50"},  # 物种数目，默认top50
            {"name": "others", "type": "float", "default": 0.05},  # 组成分析中用于将丰度小于0.05/其它的物种归为others
            {"name": "level_color", "type": "string", "default": "3"},
            {"name": "sample_distance", "type": "string", "default": "bray_curtis"},  ##样本计算距离
            {"name": "species_distance", "type": "string", "default": "bray_curtis"},  ##物种计算距离
            {"name": "fill_zero", "type": "string", "default": "true"},  # 是否对0全表用最小值进行填充
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
        ]
        self.add_option(options)
        self.analysis = self.option("analysis").split(",")
        self.analysis_list = []

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        for i in self.analysis:
            if i in ['barpie', 'heatmap']:
                print i
                break
            else:
                print i
                raise OptionError('没有选择任何分析或者分析类型选择错误：%s', variables=(self.option('analysis')))
        if not self.option("abundtable").is_set:
            raise OptionError("请传入物种/功能/基因丰度表格！")
        if not self.option("group").is_set:
            raise OptionError("请传入分组表格！")

    def run_bar_sort_samples(self):
        """
        运行barpie
        :return:
        """
        self.sort_samples = self.add_tool("metaasv.sort_samples")
        self.sort_samples.set_options({
            "in_otu_table": self.option("abundtable"),
            "group_table": self.option("group"),
            "method": self.option("add_Algorithm"),
            "others": self.option("others"),
            "fill_zero": "false",
            "abundance_method": "all"
        })
        self.analysis_list.append(self.sort_samples)

    def run_heatmap(self):
        """
        运行heatmap
        :return:
        """
        input_table = self.change_otuname(self.option("abundtable").prop['path'], "asv")
        self.heatmap = self.add_module("metaasv.heatmap")
        opts = ({
            "abundtable": input_table,
            "group": self.option("group"),
            "species_number": self.option("species_number"),
            "method": self.option("method"),
            "sample_method": self.option("sample_method"),
            "add_Algorithm": self.option("add_Algorithm"),
            "sample_distance": self.option("sample_distance"),
            "species_distance": self.option("species_distance"),
            "fill_zero": self.option("fill_zero")
            })
        if self.option("phy_newick").is_set:
            opts["phy_newick"] = self.option("phy_newick")
        self.heatmap.set_options(opts)
        self.analysis_list.append(self.heatmap)

    def run_analysis(self):
        if "heatmap" in self.analysis:
            self.run_heatmap()
        if "barpie" in self.analysis:
            self.run_bar_sort_samples()
        self.on_rely(self.analysis_list, self.set_output)
        for module in self.analysis_list:
            module.run()
            gevent.sleep(0)

    def change_otuname(self, tablepath, file_name):
        """
        改换asv名称
        :param tablepath:
        :param file_name:
        :return:
        """
        self.rename_dict = {}
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        with open(tablepath, "r") as f, open(newtable, "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                origin_name = lines[0]
                specimen = re.subn("^.*;", "", lines[0])[0]
                if specimen not in self.rename_dict:
                    self.rename_dict[specimen] = origin_name
                g.write(specimen + "\t" + lines[1])
        return newtable

    def convert_format(self, input_file, outout_file):
        """
        metaasv 筛选leve_color
        :return:
        """
        level_color = self.option('level_color')
        self.logger.info("+++++++++{}".format(self.rename_dict))
        level_color_number = int(level_color) - 1
        with open(input_file, 'r') as f, open(outout_file, "w") as w:
            lines = f.readlines()
            w.write(lines[0].strip() + "\tlevel_color\n")
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[0] in self.rename_dict.keys():
                    asv_name = self.rename_dict[line[0]]
                    asv_name_list = asv_name.strip().split(";")
                    level_name = asv_name_list[level_color_number]
                    w.write("\t".join(line) + "\t{}\n".format(level_name))

    def set_output(self):
        """
        设置结果文件
        :return:
        """
        if 'barpie' in self.analysis:
            bar_path = os.path.join(self.output_dir, "CommunityBarPie")
            link_dir(self.sort_samples.output_dir, bar_path)
            if os.path.exists(os.path.join(self.output_dir, "CommunityBarPie", "taxa.percents.mongo.xls")):
                os.remove(os.path.join(self.output_dir, "CommunityBarPie", "taxa.percents.mongo.xls"))
        if 'heatmap' in self.analysis:
            heatmap_path = os.path.join(self.output_dir, "CommunityHeatmap")
            link_dir(self.heatmap.output_dir, heatmap_path)
            if os.path.exists(os.path.join(self.output_dir, "CommunityHeatmap", "taxa.percents.mongo.xls")):
                os.remove(os.path.join(self.output_dir, "CommunityHeatmap", "taxa.percents.mongo.xls"))
            insert_asv_table = os.path.join(self.output_dir, "CommunityHeatmap", "taxa.table.xls")
            output_asv_table = os.path.join(self.work_dir, "taxa.table.mongo.xls")
            insert_asv_percents_table = os.path.join(self.output_dir, "CommunityHeatmap", "taxa.percents.table.xls")
            output_asv_percents_table = os.path.join(self.work_dir, "taxa.percents.mongo.xls")
            self.convert_format(insert_asv_table, output_asv_table)
            self.convert_format(insert_asv_percents_table, output_asv_percents_table)
            link_file(output_asv_percents_table, insert_asv_percents_table)
            link_file(output_asv_table, insert_asv_table)
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(CompositionAnalysisModule, self).run()
        self.run_analysis()

    def end(self):
        """
        结束
        :return:
        """
        super(CompositionAnalysisModule, self).end()
