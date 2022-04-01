# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify:2017.08.22
#last_modify:2018.10.24--qingchen.zhang

"""bar / circos / heatmap/bubble"""
import os
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class CompositionAnalysisModule(Module):
    def __init__(self, work_id):
        super(CompositionAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string", "default": "bar,heatmap,circos, bubble"},
            {"name": "abundtable", "type": "infile", "format": "meta.otu.otu_table"},  # 物种/功能丰度表格
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "species_number", "type": "string", "default": "50"},  # 物种数目，默认top50
            {"name": "others", "type": "float", "default": 0.01},  # 组成分析中用于将丰度小于0.01/其它的物种归为others
            {"name": "level_type", "type": "string", "default": ""},    #用于heatmap/bubble分析
            {"name": "anno_type", "type": "string", "default": ""},    #用于bubble和heatmap分析的nr的结果分离
            {"name": "project_name", "type":"string", "default": ""},   #用于区分meta和metagenomic
            {"name": "level_color", "type": "string", "default": ""},
            {"name": "normalization", "type": "string", "default": ""},  # 对数据进行标准化, row, col, 或 “”
        ]
        self.add_option(options)
        self.analysis = []

    def check_options(self):
        analysis = self.option('analysis').split(',')
        for i in analysis:
            if i in ['bar', 'circos', 'heatmap', 'bubble']:
                print(i)
                break
            else:
                print(i)
                raise OptionError('没有选择任何分析或者分析类型选择错误：%s', variables=(self.option('analysis')), code="22700301")
        if not self.option("abundtable").is_set:
            raise OptionError("请传入物种/功能/基因丰度表格！", code="22700302")
        if not self.option("group").is_set:
            raise OptionError("请传入分组表格！", code="22700303")

    def run_bar_sort_samples(self):
        if 'bar' not in self.option('analysis'):
            return
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples.set_options({
            "in_otu_table": self.option("abundtable"),
            "group_table": self.option("group"),
            "method": self.option("add_Algorithm"),
            "others": self.option("others")
        })
        self.analysis.append(self.sort_samples)

    def run_circos_sort_samples(self):
        if 'circos' not in self.option('analysis'):
            return
        self.sort_samples2 = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples2.set_options({
            "in_otu_table": self.option("abundtable"),
            "group_table": self.option("group"),
            "method": self.option("add_Algorithm"),
            "others": self.option("others")
        })
        self.analysis.append(self.sort_samples2)

    def run_heatmap(self):
        if 'heatmap' not in self.option('analysis'):
            return
        self.heatmap = self.add_module("meta.composition.heatmap")
        self.heatmap.set_options({
            "abundtable": self.option("abundtable"),
            "group": self.option("group"),
            "species_number": self.option("species_number"),
            "method": self.option("method"),
            "sample_method": self.option("sample_method"),
            "add_Algorithm": self.option("add_Algorithm"),
            "level_type": self.option("level_type"),
            "anno_type": self.option('anno_type'),
            "level_color": self.option('level_color'),
            "normalization": self.option('normalization')
            })
        self.analysis.append(self.heatmap)

    def run_bubble(self):
        """
        bubble分析
        :return:
        """
        if 'bubble' not in self.option('analysis'):
            return
        self.bubble = self.add_module("meta.composition.heatmap")
        self.bubble.set_options({
            "abundtable": self.option("abundtable"),
            "group": self.option("group"),
            "species_number": self.option("species_number"),
            "add_Algorithm": self.option("add_Algorithm"),
            "level_type": self.option("level_type"),
            "anno_type": self.option('anno_type'),
            "level_color": self.option('level_color'),
            "analysis": self.option('analysis'),
            "others": self.option('others'),

        })
        self.analysis.append(self.bubble)

    def run_venn(self):
        if 'venn' not in self.option('analysis'):
            return
        no_zero_otu = os.path.join(self.work_dir, "otu.nozero")
        if not self.option("group").is_set:
            return
        self.samples = []
        self.groups = []
        with open(self.option("group").path) as f:
            data = f.readlines()
            for i in data[1:]:
                self.samples.append(i.strip().split("\t")[0])
                self.groups.append(i.strip().split("\t")[1])
        if len(list(set(self.groups))) < 2:
            return
        self.option("abundtable").sub_otu_sample(self.samples, no_zero_otu)
        num_lines = sum(1 for line in open(no_zero_otu))
        if num_lines < 11:
            return
        self.venn = self.add_tool("graph.venn_table")
        self.venn.set_options({
            "otu_table": no_zero_otu,
            "group_table": self.option("group")
        })
        self.analysis.append(self.venn)

    def run_analysis(self):
        self.run_heatmap()
        self.run_bubble()
        self.run_circos_sort_samples()
        self.run_bar_sort_samples()
        self.run_venn()
        self.on_rely(self.analysis, self.set_output)
        for module in self.analysis:
            module.run()

    def set_output(self):
        if 'bar' in self.option('analysis'):
            self.linkdir(self.sort_samples.output_dir, 'bar')
        if 'circos' in self.option('analysis'):
            self.linkdir(self.sort_samples2.output_dir, 'circos')
        if 'heatmap' in self.option('analysis'):
            self.linkdir(self.heatmap.output_dir, 'heatmap')
        if 'bubble' in self.option('analysis'):
            self.linkdir(self.bubble.output_dir, 'bubble')
        if 'barpie' in self.option('analysis'):
            self.linkdir(self.sort_samples.output_dir, 'barpie')
        self.end()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(CompositionAnalysisModule, self).run()
        self.run_analysis()
