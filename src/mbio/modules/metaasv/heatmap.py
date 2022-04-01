# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir


class HeatmapModule(Module):
    """
    metaasv heatmap
    """
    def __init__(self, work_id):
        super(HeatmapModule, self).__init__(work_id)
        options = [
            {"name": "abundtable", "type": "infile", "format": "meta.otu.otu_table"},  # 物种/功能丰度表格
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "species_number", "type": "string", "default": "50"},  # 物种数目，默认top50
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "analysis","type": "string", "default": "heatmap"},  #作为bubble和heatmap分析的参数
            {"name": "others", "type": "float", "default": ""},
            {"name": "sample_distance", "type": "string", "default": "bray_curtis"}, ##样本计算距离
            {"name": "species_distance", "type": "string", "default": "bray_curtis"},##物种计算距离
            {"name": "fill_zero", "type": "string", "default": "true"}, #  是否对0全表用最小值进行填充
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("abundtable").is_set:
            raise OptionError("请传入物种/功能/基因丰度表格！")
        if not self.option("group").is_set:
            raise OptionError("请传入分组表格！")
        if self.option('method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的物种层次聚类方式：%s')
        if self.option('sample_method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s')
        if self.option('add_Algorithm') not in ['sum', 'average', 'middle', ""]:
            raise OptionError('错误的样本求和方式：%s')
        if self.option("method") != "":
            if self.option("species_number") != "":
                if int(self.option("species_number")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s')

    def run_sort_samples(self):
        self.sort_samples = self.add_tool("metaasv.sort_samples")
        if self.option("species_number") != "":
            number = int(self.option("species_number"))
        else:
            number = 0
        if self.option('analysis') == 'bubble':
            self.sort_samples.set_options({
                "in_otu_table": self.option("abundtable"),
                "group_table": self.option("group"),
                "method": self.option("add_Algorithm"),
                "top": number,
                "others": self.option('others'),
                "abundance_method": "all",
                "fill_zero": self.option("fill_zero")
            })
        else:
            self.sort_samples.set_options({
                "in_otu_table": self.option("abundtable"),
                "group_table": self.option("group"),
                "method": self.option("add_Algorithm"),
                "top": number,
                "abundance_method": "all",
                "fill_zero": self.option("fill_zero")
            })
        if self.option("sample_method") == "":
            if self.option("method") == "":
                self.sort_samples.on("end", self.end)
            else:
                self.sort_samples.on("end", self.run_species_matrix)
        else:
            self.sort_samples.on("end", self.run_matrix)
        self.sort_samples.run()

    def run_matrix(self):
        """
        运行计算距离矩阵
        :return:
        """
        new_otu_file_path = self.sort_samples.option("out_otu_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        sample = content[0].strip().split('\t')
        if len(sample) < 3 :
            raise OptionError('样本聚类的个数不能为1', code="22700407")
        self.logger.info("正在进行样本距离计算")
        self.matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.matrix.set_options({'method': self.option("sample_distance"),
                                'otutable': new_otu_file_path
                                 })
        self.matrix.on('end', self.run_hcluster)
        self.matrix.run()

    def run_hcluster(self):
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.hcluster.set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'linkage': self.option("sample_method"),
        })
        if self.option("method") == "":
            self.hcluster.on('end', self.end)
        else:
            self.hcluster.on('end', self.run_species_matrix)
        self.hcluster.run()

    def transposition(self,old, path):
        """
        转置一个丰度表
        """
        file_ = list()
        with open(old, 'rb') as r:
            linelist = [l.strip('\r\n') for l in r.readlines()]
        for row in linelist:
            row = re.split("\t", row)
            file_.append(row)
        zip_line = zip(*file_)
        with open(path, 'wb') as w:
            for my_l in zip_line:
                w.write("\t".join(my_l) + "\n")

    def run_species_matrix(self):
        new_otu_file_path = self.sort_samples.option("out_otu_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        if len(content) == 2:
            raise OptionError('物种/功能的个数小于2，无法进行物种/功能聚类分析')
        self.logger.info("正在进行物种/功能/基因距离计算")
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        self.transposition(new_otu_file_path, trans_otu)
        self.species_matrix = self.add_tool("meta.beta_diversity.distance_calc")
        opts = ({
            "method": self.option("species_distance"),
            "otutable": trans_otu
        })
        if self.option("phy_newick").is_set:
            opts['newicktree'] = self.option('phy_newick')
        self.species_matrix.set_options(opts)
        self.species_matrix.on('end', self.run_species_cluster)
        self.species_matrix.run()

    def run_species_cluster(self):
        self.species_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.species_hcluster.set_options({
            'dis_matrix': self.species_matrix.option('dis_matrix'),
            'linkage': self.option("method")
        })
        self.species_hcluster.on('end', self.end)
        self.species_hcluster.run()

    def run(self):
        super(HeatmapModule, self).run()
        self.run_sort_samples()

    def end(self):
        """
        生成结果文件目录
        :return:
        """
        if self.option('analysis')== 'heatmap':
            link_dir(self.sort_samples.output_dir, self.output_dir)
        if self.option("sample_method") != "":
            if os.path.exists(self.output_dir + "/specimen_hcluster.tre"):
                os.remove(self.output_dir + "/specimen_hcluster.tre")
            os.link(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + "/specimen_hcluster.tre")
        if self.option("method") != "":
            if os.path.exists(self.output_dir + "/species_hcluster.tre"):
                os.remove(self.output_dir + "/species_hcluster.tre")
            os.link(self.species_hcluster.output_dir+"/hcluster.tre", self.output_dir + "/species_hcluster.tre")
        super(HeatmapModule, self).end()
