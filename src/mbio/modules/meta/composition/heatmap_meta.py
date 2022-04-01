# -*- coding: utf-8 -*-

import os
import re
import pandas as pd
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import shutil


class HeatmapMetaModule(Module):
    """
    和交互分析一样的群落heatmap图
    """
    def __init__(self, work_id):
        super(HeatmapMetaModule, self).__init__(work_id)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "level", "type": "string", "default": "7"},  # 输入的OTU level
            {"name": "species_number", "type": "string", "default": ""},  # 物种数目，默认全部物种
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "level_color", "type": "string", "default": "3"}
        ]
        self.add_option(options)
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.sample_matrix = self.add_tool("meta.beta_diversity.distance_calc") # 20161206 2 lines
        self.sample_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        group_table_path = os.path.join(self.work_dir, "group_table.xls")

        self.s2_file_path = ""
        self.color_dict = {}

    def check_options(self):
        if self.option('method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的物种层次聚类方式：%s', variables=(self.option('method')), code="12702001")
        if self.option('sample_method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s', variables=(self.option('sample_method')), code="12702002")
        if self.option('add_Algorithm') not in ['sum', 'average', 'middle', ""]:
            raise OptionError('错误的样本求和方式：%s', variables=(self.option('add_Algorithm')), code="12702003")
        # if (self.option("method") != "" and self.option("species_number") != ("" and "all")):
        if self.option("method") != "":
            if (self.option("species_number") != "" and self.option("species_number") != "0"):
                if int(self.option("species_number")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s', variables=(self.option('species_number')), code="12702004")

    def run_sort_samples(self):  # modify by zhujuan 20171127 先挑选样品，后取top物种（get_species）
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path,
            "method": self.option("add_Algorithm"),
            "top_n": self.option("species_number") or 0
        })
        self.sort_samples.on('end', self.get_species)
        self.sort_samples.run()

    def get_species(self):
        old_otu_file_path = self.sort_samples.option("out_otu_table").prop['path']
        if self.option("sample_method") != "":  # 判断是否能做聚类
            print(self.group_table_path)
            sample_n = open(self.group_table_path, "r")
            content = sample_n.readlines()
            if len(content) == 2:
                raise OptionError('样本聚类的个数不能为：1', code="12702005")
            if self.option("add_Algorithm") != "":
                sample_ = []
                sample_n2 = open(self.group_table_path, "r")
                content = sample_n2.readlines()
                for f in content:
                    f = f.strip("\n")
                    arr = f.strip().split("\t")
                    if arr[0] != "#sample":
                        if arr[1] not in sample_:
                            sample_.append(arr[1])
                if len(sample_) == 1:
                    raise OptionError('当计算分组丰度并且进行样本聚类的分析时,样本分组不能为1', code="12702006")
        self.s2_file_path = os.path.join(self.work_dir, "s_otu_table.xls")
        otu_relative_path = self.sort_samples.option('level_otu_table').path
        self.otu_relative = os.path.join(self.output_dir, 'heatmap.taxa.relative.xls')
        with open(old_otu_file_path, 'r') as r, open(self.s2_file_path, 'w') as w:
            for n in r:
                n = n.strip("\n")
                brr = n.strip().split(";")
                if brr[0].startswith("OTU"):
                    w.write(brr[0] + "\n")
                else:
                    w.write(brr[-1] + "\n")
                    # by houshuang 20191010 分类水平选择domain时没有颜色水平>>>
                    if self.option("level_color") != "":
                        name = re.split("\t", brr[-1])
                        self.color_dict[name[0]] = brr[int(self.option("level_color")) - 1].strip(";").strip()
                    # <<<
            self.new_otu_file_path = self.s2_file_path

        with open(otu_relative_path, 'r') as r, open(self.otu_relative, 'w') as w:
            for n in r:
                n = n.strip("\n")
                brr = n.strip().split(" ")
                if brr[0] == "OTU":
                    w.write(brr[0] + " " + brr[1] + "\n")
                else:
                    w.write(brr[-1] + "\n")

        if self.option("method") == "":
            if self.option("sample_method") == "":
                self.end()
            else:
                self.run_sample_matrix()
        else:
            self.run_matrix()

    def transposition(self, old, path):
        """
        转置一个otu表
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

    def run_matrix(self):
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        self.transposition(self.new_otu_file_path, trans_otu)
        self.matrix.set_options({
            "method": "bray_curtis",
            "otutable": trans_otu
        })
        self.matrix.on('end', self.run_cluster)
        self.matrix.run()

    def run_cluster(self):
        options = {
            "dis_matrix": self.matrix.option('dis_matrix'),
            "linkage": self.option("method")
        }
        self.hcluster.set_options(options)
        if self.option("sample_method") != "":
            self.hcluster.on('end', self.run_sample_matrix)
        else:
            self.hcluster.on('end', self.end)
        self.hcluster.run()

    def run_sample_matrix(self):  # 20161206
        self.logger.info("正在进行样本距离计算")
        options = {
            "method": "bray_curtis",
            "otutable": self.new_otu_file_path
        }
        self.sample_matrix.set_options(options)
        self.sample_matrix.on('end', self.run_sample_hcluster)
        self.logger.info("样本距离计算结束")
        self.sample_matrix.run()

    def run_sample_hcluster(self):
        self.logger.info("正在进行样本聚类计算")
        options = {
            "dis_matrix": self.sample_matrix.option('dis_matrix'),
            "linkage": self.option("sample_method")
        }
        self.sample_hcluster.set_options(options)
        self.sample_hcluster.on('end', self.end)
        self.logger.info("样本聚类计算结束")
        self.sample_hcluster.run()

    def replace_0_new_otu_file_path(self, infile):
        """
        将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一
        :return:
        """
        self.logger.info("正在将结果文件中的丰度值0替换")
        min_list = []  # add modify by qingchen.zhang@20190306 用于增加
        second_list = []
        real_zero_otu_new = self.work_dir + "/real_otu_new.xls"
        with open(infile, 'r') as f, open(real_zero_otu_new, 'w') as w:
            lines = f.readlines()
            w.write("{}".format(lines[0]))
            for line in lines[1:]:
                line = line.strip().split('\t')
                global line, min
                for i in range(1, len(line)):
                    if line[i] == "":
                        continue
                    if float(line[i]) != 0.0:
                        line_min = float(line[i])
                        min_list.append(line_min)
                        min_line = min(min_list)
                second_list.append(min_line)
            min_table = float(min(second_list))
            for line in lines[1:]:
                line = line.strip().split('\t')
                data_list = []
                for i in range(1, len(line)):
                    if line[i] == "":
                        continue
                    if float(line[i]) == 0.0:
                        line[i] = float(line[i]) + float(min_table / 10)
                    else:
                        line[i] = float(line[i])
                    data_list.append(line[i])
                w.write("{}\t{}\n".format(line[0], "\t".join(str(i) for i in data_list)))
        os.rename(real_zero_otu_new, infile)

    def replace_0_new_otu_file_path2(self, infile, outfile):
        """
        将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一
        :return:
        """
        self.logger.info("正在将结果文件中的丰度值0替换")
        min_list = []  # add modify by qingchen.zhang@20190306 用于增加
        second_list = []
        real_zero_otu_new = self.work_dir + "/real_otu_new.xls"
        with open(infile, 'r') as f, open(outfile, 'w') as w:
            lines = f.readlines()
            w.write("{}".format(lines[0]))
            for line in lines[1:]:
                line = line.strip().split('\t')
                global line, min
                for i in range(1, len(line)):
                    if line[i] == "":
                        continue
                    if float(line[i]) != 0.0:
                        line_min = float(line[i])
                        min_list.append(line_min)
                        min_line = min(min_list)
                second_list.append(min_line)
            min_table = float(min(second_list))
            for line in lines[1:]:
                line = line.strip().split('\t')
                data_list = []
                for i in range(1, len(line)):
                    if line[i] == "":
                        continue
                    if float(line[i]) == 0.0:
                        line[i] = float(line[i]) + float(min_table / 10)
                    else:
                        line[i] = float(line[i])
                    data_list.append(line[i])
                w.write("{}\t{}\n".format(line[0], "\t".join(str(i) for i in data_list)))

    def run(self):
        super(HeatmapMetaModule, self).run()
        self.group_table_path = self.option('group').prop["path"]
        self.run_sort_samples()

    def end(self):
        #anno_type = self.option('anno_type')
        self.replace_0_new_otu_file_path(self.otu_relative)
        self.replace_0_new_otu_file_path2(self.sort_samples.option("out_otu_table").prop['path'], self.output_dir + "/sort_sample.xls")
        if os.path.exists(os.path.join(self.output_dir, "heatmap.taxa.table.xls")):
            os.remove(os.path.join(self.output_dir, "heatmap.taxa.table.xls"))
        os.link(self.new_otu_file_path, os.path.join(self.output_dir, "heatmap.taxa.table.xls"))
        if os.path.exists(self.sample_hcluster.output_dir + "/hcluster.tre"):
            if os.path.exists(self.output_dir + "/sample_hcluster.tre"):
                os.remove(self.output_dir + "/sample_hcluster.tre")
            shutil.copy(self.sample_hcluster.output_dir + "/hcluster.tre", self.output_dir + "/sample_hcluster.tre")
        if os.path.exists(self.hcluster.output_dir + "/hcluster.tre"):
            if os.path.exists(self.output_dir + "/species_hcluster.tre"):
                os.remove(self.output_dir + "/species_hcluster.tre")
            shutil.copy(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + "/species_hcluster.tre")
        super(HeatmapMetaModule, self).end()
