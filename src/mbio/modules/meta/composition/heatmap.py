# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 20180102
#last_modify: 20181107--qingchen.zhang：将最后的taxa_table表格做一步处理之后再上传


import os
import re
import pandas as pd
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class HeatmapModule(Module):
    """
    群落heatmap/bubble
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
            {"name": "level_type", "type": "string"},  #作为split_samples分析的输入文件的参数
            {"name": "anno_type", "type":"string", "default": ""},  #主要是用于整理nr分析内容的结果用
            {"name": "level_color", "type": "string", "default": ""},
            {"name": "others", "type": "float", "default": ""},
            {"name": "normalization", "type": "string", "default": ""},  # 对数据进行标准化, row, col, 或 “”
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("abundtable").is_set:
            raise OptionError("请传入物种/功能/基因丰度表格！", code="22700401")
        if not self.option("group").is_set:
            raise OptionError("请传入分组表格！", code="22700402")
        if self.option('method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的物种层次聚类方式：%s', variables=(self.option('method')), code="22700403")
        if self.option('sample_method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s', variables=(self.option('sample_method')), code="22700404")
        if self.option('add_Algorithm') not in ['sum', 'average', 'middle', ""]:
            raise OptionError('错误的样本求和方式：%s', variables=(self.option('add_Algorithm')), code="22700405")
        if self.option("method") != "":
            if self.option("species_number") != "0":
                if int(self.option("species_number")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s', variables=(self.option('species_number')), code="22700406")

    def run_sort_samples(self):
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        number = int(self.option("species_number"))
        if self.option('analysis') == 'bubble':
            self.sort_samples.set_options({
                "in_otu_table": self.option("abundtable"),
                "group_table": self.option("group"),
                "method": self.option("add_Algorithm"),
                "top": number,
                "others": self.option('others'),
            })
        else:
            self.sort_samples.set_options({
                "in_otu_table": self.option("abundtable"),
                "group_table": self.option("group"),
                "method": self.option("add_Algorithm"),
                "top": number,
                "normalization": self.option('normalization'),
            })
        #self.sort_samples.on('end',self.get_species)
        if self.option("sample_method") == "":
            if self.option("method") == "":
                self.sort_samples.on("end", self.convert_format)
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
        self.matrix.set_options({'method': "bray_curtis",
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
            self.hcluster.on('end', self.convert_format)
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
            raise OptionError('物种/功能的个数小于2，无法进行物种/功能聚类分析', code="22700408")
        self.logger.info("正在进行物种/功能/基因距离计算")
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        self.transposition(new_otu_file_path, trans_otu)
        self.species_matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.species_matrix.set_options({
            "method": "bray_curtis",
            "otutable": trans_otu
        })
        self.species_matrix.on('end', self.run_species_cluster)
        self.species_matrix.run()

    def run_species_cluster(self):
        self.species_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.species_hcluster.set_options({
            'dis_matrix': self.species_matrix.option('dis_matrix'),
            'linkage': self.option("method")
        })
        self.species_hcluster.on('end', self.convert_format)
        self.species_hcluster.run()

    def run(self):
        super(HeatmapModule, self).run()
        self.run_sort_samples()

    def convert_format(self):
        """
        宏基因组对bubble和heatmap转level_color
        :return:
        """
        if self.option('analysis')== 'bubble':
            level_color = self.option('level_color')
            taxa_table_path = os.path.join(self.sort_samples.output_dir, 'taxa.percents.table.xls')
            taxa_table = pd.read_table(taxa_table_path, sep='\t', header=0)
            split_table= pd.DataFrame(taxa_table)
            new = split_table.iloc[-1]
            new_list = pd.DataFrame(new)
            well = new_list.T
            new_split = pd.concat([well,split_table],axis=0)
            new_split = new_split.reset_index()
            del new_split['index']
            split_table = new_split.drop(list(new_split.index)[-1:])
            if level_color != "":
                split_table[self.option('level_type')], split_table[self.option('level_color')] = split_table[self.option('level_type')].str.split('|', 1).str
                split_table.to_csv(os.path.join(self.sort_samples.output_dir, 'taxa.percents.table.xls'), index=False, sep='\t', encoding="utf-8")
            else:
                pass
        elif self.option('analysis')== 'heatmap':
            level_color = self.option('level_color')
            taxa_table_path = os.path.join(self.sort_samples.output_dir, 'taxa.table.xls')
            taxa_table = pd.read_table(taxa_table_path, sep='\t', header=0)
            split_table= pd.DataFrame(taxa_table)
            if level_color != "":
                if self.option("sample_method") == "" and self.option("method") == "":
                    split_table[self.option('level_type')], split_table[self.option('level_color')] = split_table[self.option('level_type')].str.split('|', 1).str
                    split_table = split_table.sort_values(by=[self.option('level_color')], ascending=True)
                    split_table.to_csv(os.path.join(self.sort_samples.output_dir, 'taxa.table.xls'), index=False, sep='\t', encoding="utf-8")
                else:
                    split_table[self.option('level_type')], split_table[self.option('level_color')] = split_table[self.option('level_type')].str.split('|', 1).str
                    split_table.to_csv(os.path.join(self.sort_samples.output_dir, 'taxa.table.xls'), index=False, sep='\t', encoding="utf-8")
            else:
                pass
        self.end()

    def end(self):
        #anno_type = self.option('anno_type')
        if self.option('analysis')== 'bubble':
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):  # 先判断是否已有结果，如果已存在，先删除再链接  modified by guhaidong 20180102
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.sort_samples.output_dir + "/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            if os.path.exists(self.output_dir + "/taxa.table.xls"):  # 先判断是否已有结果，如果已存在，先删除再链接  modified by guhaidong 20180102
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.sort_samples.output_dir + "/taxa.table.xls", self.output_dir + "/taxa.table.xls")
        else:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):  # 先判断是否已有结果，如果已存在，先删除再链接  modified by guhaidong 20180102
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.sort_samples.output_dir + "/taxa.table.xls", self.output_dir + "/taxa.table.xls")
        if self.option("sample_method") != "":
            if os.path.exists(self.output_dir + "/specimen_hcluster.tre"):
                os.remove(self.output_dir + "/specimen_hcluster.tre")
            os.link(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + "/specimen_hcluster.tre")
        if self.option("method") != "":
            if os.path.exists(self.output_dir + "/species_hcluster.tre"):
                os.remove(self.output_dir + "/species_hcluster.tre")
            os.link(self.species_hcluster.output_dir+"/hcluster.tre", self.output_dir + "/species_hcluster.tre")
        super(HeatmapModule, self).end()
