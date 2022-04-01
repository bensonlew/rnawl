# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

"""聚类heatmap图模块"""
import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class HcHeatmapModule(Module):
    """
    小工具： 聚类热图
    参照多样性-群落heatmap图
    """
    def __init__(self, work_id):
        super(HcHeatmapModule, self).__init__(work_id)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入的group表
            {"name": "top_number", "type": "int", "default": ""},  # 物种数目，默认全部物种
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "group_method", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "is_standard", "type": "bool", "default": False}, # 是否进行标准化“对0进行替换”， 默认不做标准化
            {"name": "distance", "type": "string", "default": "bray_curtis"}, ## 计算物种的距离，默认为欧氏距离
            {"name": "specimen_distance", "type": "string", "default": "bray_curtis"} ## 计算样本的距离，默认为欧氏距离
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if self.option('method') not in ['average', 'single', 'complete', "", "none"]:
            raise OptionError('错误的物种层次聚类方式：%s', self.option('method'), code="24400101")
        if self.option('sample_method') not in ['average', 'single', 'complete', "", "none"]:
            raise OptionError('错误的样本层次聚类方式：%s', self.option('sample_method'), code="24400102")
        if self.option('group_method') not in ['sum', 'average', 'middle', "", "none"]:
            raise OptionError('错误的样本求和方式：%s', self.option('group_method'), code="24400103")
        if int(self.option("top_number")) == 1:
            raise OptionError('物种聚类的个数不能为：%s', self.option('top_number'), code="24400104")
        
    def run_sort_samples(self):  
        """
        功能：对样本进行合并和排序
        :return:
        """
        self.logger.info("正在对样本进行合并和排序！")
        self.sort_samples = self.add_tool("meta.sort_samples_toolapps")
        opts = {
            "in_otu_table": self.option("otu_table"),
            "group_table": self.option("group_table"),
            "top": int(self.option("top_number"))
        }
        if self.option("group_method") in ['','none']:
            opts["method"] = ''
        else:
            opts["method"] = self.option("group_method")
        self.sort_samples.set_options(opts)
        if self.option("sample_method") in ["", 'none']:
            if self.option("method") in ["", 'none']:
                self.sort_samples.on("end", self.set_output)
            else:
                self.sort_samples.on("end", self.run_species_matrix)
        else:
            self.sort_samples.on('end',self.run_matrix)
        self.sort_samples.run()

    def run_matrix(self):
        """
        功能：计算距离矩阵
        :return:
        """
        new_otu_file_path = self.sort_samples.option("out_otu_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        sample = content[0].strip().split('\t')
        if len(sample) < 3 :
            raise OptionError('样本聚类的个数不能为1', code="24400105")
        self.logger.info("正在进行样本距离计算!")
        self.matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.matrix.set_options({
            'method': self.option("specimen_distance"),
            'otutable': new_otu_file_path
            })
        self.matrix.on('end', self.run_hcluster)
        self.matrix.run()

    def run_hcluster(self):
        """
        功能：对矩阵进行聚类
        :return:
        """
        self.logger.info("正在进行样本聚类计算!")
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        opts = {
            'dis_matrix': self.matrix.option('dis_matrix'),
            }
        if self.option("sample_method") in ['','none']:
            opts["linkage"] = ''
        else:
            opts["linkage"] = self.option("sample_method")
        self.hcluster.set_options(opts)
        if self.option("method") in ["", 'none']:
            self.hcluster.on('end', self.set_output)
        else:
            self.hcluster.on('end', self.run_species_matrix)
        self.hcluster.run()

    def transposition(self,old, path):
        """
        转置丰度表
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
        功能：计算物种或者行的矩阵
        :return:
        """
        new_otu_file_path = self.sort_samples.option("out_otu_table").prop['path']
        sample_n = open(new_otu_file_path, "r")
        content = sample_n.readlines()
        if len(content) == 2:
            raise OptionError('物种/功能的个数小于2，无法进行物种/功能聚类分析', code="24400106")
        self.logger.info("正在进行物种/功能/基因距离计算")
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        self.transposition(new_otu_file_path, trans_otu)
        self.species_matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.species_matrix.set_options({
            "method": self.option("distance"),
            "otutable": trans_otu
        })
        self.species_matrix.on('end', self.run_species_cluster)
        self.species_matrix.run()

    def run_species_cluster(self):
        """
        功能：对矩阵进行聚类
        :return:
        """
        self.logger.info("正在进行物种聚类计算!")
        self.species_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        opts = {
            'dis_matrix': self.species_matrix.option('dis_matrix')
            }
        if self.option("method") in ['','none']:
            opts["linkage"] = ''
        else:
            opts["linkage"] = self.option("method")
        self.species_hcluster.set_options(opts)
        self.species_hcluster.on('end', self.set_output)
        self.species_hcluster.run()

    def replace_0_new_otu_file_path(self,in_table, out_table):
        """
        将丰度表中的所有绝对丰度为0的结果在此基础上加上整张丰度表最小值的十分之一
        :return:
        """
        self.logger.info("正在对otu结果进行标准化")
        min_list = []
        all_min_list = []
        new_otu_file_path = in_table
        real_zero_otu_new = out_table
        with open(new_otu_file_path, 'r') as f, open(real_zero_otu_new, 'w') as w:
            lines = f.readlines()
            w.write("{}".format(lines[0]))
            for line in lines[1:]:
                line = line.strip().split('\t')
                for i in range(1, len(line)):
                    if float(line[i]) != 0.0:
                        line_min = line[i]
                        if line_min not in min_list:
                            min_list.append(line_min)
                min_number = min(min_list)
                all_min_list.append(min_number)
            min_table = float(min(all_min_list))
            for lin in lines[1:]:
                lin = lin.strip().split('\t')
                data_list = []
                for i in range(1, len(lin)):
                    if float(lin[i]) == 0.0:
                        lin[i] = float(lin[i]) + float(min_table)/10
                    else:
                        lin[i] = float(lin[i])
                    data_list.append(lin[i])
                w.write("{}\t{}\n".format(lin[0], "\t".join(str(i) for i in data_list)))
        self.logger.info("标准化完成!")
        return real_zero_otu_new

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        if self.option("is_standard"):
            new_otu_file_path = self.replace_0_new_otu_file_path(os.path.join(self.sort_samples.output_dir,"taxa.table.xls"), os.path.join(self.work_dir, "taxa.table.xls"))
            new_otu_file_persent_path = self.replace_0_new_otu_file_path(os.path.join(self.sort_samples.output_dir,"taxa.percents.table.xls"), os.path.join(self.work_dir, "taxa.percents.table.xls"))
        else:
            new_otu_file_path = os.path.join(self.sort_samples.output_dir,"taxa.table.xls")
            new_otu_file_persent_path = os.path.join(self.sort_samples.output_dir,"taxa.percents.table.xls")
        absolute = os.path.join(self.output_dir, "heatmap.taxa.table.xls")
        if os.path.exists(absolute):
            os.remove(absolute)
        os.link(new_otu_file_path, absolute)
        persent = os.path.join(self.output_dir, "heatmap.taxa.percents.table.xls")
        if os.path.exists(persent):
            os.remove(persent)
        os.link(new_otu_file_persent_path, persent)
        if self.option("sample_method") not in ['','none']:
            if os.path.exists(os.path.join(self.output_dir, "specimen_hcluster.tre")):
                os.remove(os.path.join(self.output_dir, "specimen_hcluster.tre"))
            os.link(os.path.join(self.hcluster.output_dir, "hcluster.tre"), os.path.join(self.output_dir, "specimen_hcluster.tre"))
        if self.option("method") not in ['','none']:
            if os.path.exists(os.path.join(self.output_dir, "species_hcluster.tre")):
                os.remove(os.path.join(self.output_dir, "species_hcluster.tre"))
            os.link(os.path.join(self.species_hcluster.output_dir, "hcluster.tre"), os.path.join(self.output_dir, "species_hcluster.tre"))
        self.logger.info("设置结果文件目录完成")
        self.end()

    def end(self):
        """
        结束
        :return:
        """
        self.logger.info("开始上传结果文件")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "群落Heatmap分析结果输出目录", 0, ""],
            ["./heatmap.taxa.table.xls", "xls", "群落Heatmap分析绝对丰度数据表", 0, ""],
            ["./heatmap.taxa.percents.table.xls", "xls", "群落Heatmap分析百分比数据表", 0, ""],
            ["./sample_hcluster.tre", "tre", "样本聚类树", 0, ""],
            ["./species_hcluster.tre", "tre", "物种聚类树", 0, ""]
        ])
        super(HcHeatmapModule, self).end()

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行module！")
        super(HcHeatmapModule, self).run()
        self.run_sort_samples()
