# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class PearsonsCorrelationModule(Module):
    """
    metaasv pearsons_correlation分析
    """
    def __init__(self, work_id):
        super(PearsonsCorrelationModule, self).__init__(work_id)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},  # 物种/功能丰度表格
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"}, #环境因子表
            {"name": "method", "type": "string", "default": "pearsonr"},##选择的相关性计算方法
            {"name": "env_cluster", "type": "string", "default": ""},  # 默认不聚类 环境因子
            {"name": "species_cluster", "type": "string", "default": ""},  # 默认不聚类 物种
            {"name": "env_distance", "type": "string", "default": "bray_curtis"}, ##环境因子计算距离
            {"name": "species_distance", "type": "string", "default": "bray_curtis"},##物种计算距离
            {"name": "top_species", "type": "int", "default": 0},
        ]
        self.add_option(options)
        self.correlation = self.add_tool('statistical.pearsons_correlation')

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError("请传入物种/功能/基因丰度表格！")
        if not self.option("envtable").is_set:
            raise OptionError("请传入环境因子分组表格！")
        if self.option('method') not in ['pearsonr', 'spearmanr', 'kendalltau', ""]:
            raise OptionError('错误的物种层次聚类方式：%s' % self.option('method'))
        if self.option('env_cluster') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s'%self.option('env_cluster'))
        if self.option('species_cluster') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s'%self.option('species_cluster'))
        if self.option("species_cluster") != "":
            if self.option("top_species") != "":
                if int(self.option("top_species")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s')

    def run_correlation(self):
        """
        计算相关性
        :return:
        """
        env_cluster = self.option("env_cluster")
        species_cluster = self.option("species_cluster")
        if env_cluster == "":
            env_cluster = "average"
        if species_cluster == "":
            species_cluster = "average"
        options = {
            'otutable': self.option('otutable'),
            'envtable': self.option('envtable'),
            "method": self.option('method'),
            "env_cluster": env_cluster,
            "species_cluster": species_cluster,
            "top_species": self.option('top_species'),
            }

        self.correlation.set_options(options)
        if self.option("env_cluster") == "":
            if self.option("species_cluster") == "":
                self.correlation.on("end", self.end)
            else:
                self.correlation.on("end", self.run_species_matrix)
        else:
            self.correlation.on("end", self.run_matrix)
        self.correlation.run()

    def run_matrix(self):
        """
        运行计算距离矩阵
        :return:
        """
        new_otu_file_path = self.correlation.option("cor_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        sample = content[0].strip().split('\t')
        if len(sample) < 3 :
            raise OptionError('样本聚类的个数不能为1')
        self.logger.info("正在进行样本距离计算")
        self.matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.matrix.set_options({'method': self.option("env_distance"),
                                'otutable': new_otu_file_path
                                 })
        self.matrix.on('end', self.run_hcluster)
        self.matrix.run()

    def run_hcluster(self):
        """
        构建树
        :return:
        """
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.hcluster.set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'linkage': self.option("env_cluster"),
        })
        if self.option("species_cluster") == "":
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
        """
        计算物种的距离
        :return:
        """
        new_otu_file_path = self.correlation.option("cor_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        if len(content) == 2:
            raise OptionError('物种/功能的个数小于2，无法进行物种/功能聚类分析')
        self.logger.info("正在进行物种/功能/基因距离计算")
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        self.transposition(new_otu_file_path, trans_otu)
        self.species_matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.species_matrix.set_options({
            "method": self.option("species_distance"),
            "otutable": trans_otu
        })
        self.species_matrix.on('end', self.run_species_cluster)
        self.species_matrix.run()

    def run_species_cluster(self):
        """
        构建物种树
        :return:
        """
        self.species_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.species_hcluster.set_options({
            'dis_matrix': self.species_matrix.option('dis_matrix'),
            'linkage': self.option("species_cluster")
        })
        self.species_hcluster.on('end', self.end)
        self.species_hcluster.run()

    def run(self):
        super(PearsonsCorrelationModule, self).run()
        self.run_correlation()

    def end(self):
        """
        生成结果文件目录
        :return:
        """
        link_dir(self.correlation.output_dir, self.output_dir)
        if self.option("env_cluster") != "":
            if os.path.exists(self.output_dir + "/env_tree.tre"):
                os.remove(self.output_dir + "/env_tree.tre")
            os.link(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + "/env_tree.tre")
        if self.option("species_cluster") != "":
            if os.path.exists(self.output_dir + "/species_tree.tre"):
                os.remove(self.output_dir + "/species_tree.tre")
            os.link(self.species_hcluster.output_dir+"/hcluster.tre", self.output_dir + "/species_tree.tre")
        super(PearsonsCorrelationModule, self).end()
