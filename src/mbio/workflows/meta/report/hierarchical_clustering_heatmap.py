# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan '

"""聚类heatmap图模块"""
import os
import json
import datetime
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mainapp.models.mongo.public.meta.meta import Meta
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params

class HierarchicalClusteringHeatmapWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HierarchicalClusteringHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "level", "type": "string", "default": "9"},  # 输入的OTU level
            {"name": "group_detail", "type": "string"},  # 输入的group_detail 示例如下
            # {"A":["578da2fba4e1af34596b04ce","578da2fba4e1af34596b04cf","578da2fba4e1af34596b04d0"],"B":["578da2fba4e1af34596b04d1","578da2fba4e1af34596b04d3","578da2fba4e1af34596b04d5"],"C":["578da2fba4e1af34596b04d2","578da2fba4e1af34596b04d4","578da2fba4e1af34596b04d6"]}
            {"name": "species_number", "type": "string", "default": ""},  # 物种数目，默认全部物种
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level_color", "type": "string", "default": ""}  # by houshuang 20190820
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.matrix = self.add_tool("meta.beta_diversity.distance_calc")
        self.sample_matrix = self.add_tool("meta.beta_diversity.distance_calc") # 20161206 2 lines
        self.sample_hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        group_table_path = os.path.join(self.work_dir, "group_table.xls")
        meta = Meta()
        meta._config = self.config  # 兼容不同mongo库版本
        self.group_table_path = meta.group_detail_to_table(self.option("group_detail"), group_table_path)
        self.s2_file_path = ""  #
        self.color_dict = {}  # by houshuang 20191010

    def check_options(self):
        if self.option('method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的物种层次聚类方式：%s', variables=(self.option('method')), code="12702001")
        if self.option('sample_method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s', variables=(self.option('sample_method')), code="12702002")
        if self.option('add_Algorithm') not in ['sum', 'average', 'middle', ""]:
            raise OptionError('错误的样本求和方式：%s', variables=(self.option('add_Algorithm')), code="12702003")
        # if (self.option("method") != "" and self.option("species_number") != ("" and "all")):
        if self.option("method") != "" :
            if (self.option("species_number") != "" and self.option("species_number") != "0"):
                if int(self.option("species_number")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s', variables=(self.option('species_number')), code="12702004")
        # if self.option("sample_method") != "":
        #     print(self.group_table_path)
        #     sample_n = open(self.group_table_path, "r")
        #     content = sample_n.readlines()
        #     if len(content) == 2:
        #         raise OptionError('样本聚类的个数不能为：1' )
        #     if self.option("add_Algorithm") != "":
        #         sample_ = []
        #         sample_n2 = open(self.group_table_path, "r")
        #         content = sample_n2.readlines()
        #         for f in content:
        #             f = f.strip("\n")
        #             arr = f.strip().split("\t")
        #             if arr[0] != "#sample":
        #                 if arr[1] not in sample_:
        #                     sample_.append(arr[1])
        #         if len(sample_) == 1:
        #             raise OptionError('当计算分组丰度并且进行样本聚类的分析时,样本分组不能为1' )

    def run_sort_samples(self):  # modify by zhujuan 20171127 先挑选样品，后取top物种（get_species）
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            # "group_table": self.option("group_detail"),
            "group_table": self.group_table_path,
            "method": self.option("add_Algorithm"),
            "top_n": self.option("species_number") or 0
        })
        self.sort_samples.on('end',self.get_species)
        self.sort_samples.run()

    def get_species(self):
        global new_otu_file_path
        old_otu_file_path =self.sort_samples.option("out_otu_table").prop['path']
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

        # if self.option("species_number") == "all":
        #     self.s2_file_path = os.path.join(self.work_dir, "s_otu_table.xls")
        #     s1 = open(old_otu_file_path, 'r')
        #     contento = s1.readlines()
        #     for n in contento:
        #         n = n.strip("\n")
        #         brr = n.strip().split(" ")
        #         if brr[0] == "OTU":
        #             with open(self.s2_file_path, "a") as w:
        #                 w.write(brr[0] + " " + brr[1] + "\n")
        #         else:
        #             with open(self.s2_file_path, "a") as w:
        #                 w.write(brr[-1] + "\n")
        #     new_otu_file_path = self.s2_file_path

        self.s2_file_path = os.path.join(self.work_dir, "s_otu_table.xls")
        otu_relative_path = self.sort_samples.option('level_otu_table').path
        self.otu_relative = os.path.join(self.output_dir, 'heatmap.taxa.relative.xls')
        with open(old_otu_file_path, 'r') as r, open(self.s2_file_path, 'w') as w:
            for n in r:
                n = n.strip("\n")
                brr = n.strip().split(" ")
                if brr[0] == "OTU":
                    w.write(brr[0] + " " + brr[1] + "\n")
                else:
                    w.write(brr[-1] + "\n")
                    # by houshuang 20191010 分类水平选择domain时没有颜色水平>>>
                    if self.option("level_color") != "":
                        name = re.split("\t", brr[-1])
                        self.color_dict[name[0]] = brr[int(self.option("level_color"))-1].strip(";").strip()
                    # <<<
            new_otu_file_path = self.s2_file_path
        
        with open(otu_relative_path, 'r') as r, open(self.otu_relative, 'w') as w:
            for n in r:
                n = n.strip("\n")
                brr = n.strip().split(" ")
                if brr[0] == "OTU":
                    w.write(brr[0] + " " + brr[1] + "\n")
                else:
                    w.write(brr[-1] + "\n")

        # list2 = [] # 放入sum值 为重复做准备
        # if (self.option("species_number") != "0" and self.option("species_number") != ""):
        #     species_nu = int(self.option("species_number"))
        #     if species_nu >= all:
        #         new_otu_file_path = self.s2_file_path
        #     else:
        #         new_otu_file_path = os.path.join(self.work_dir, "new_otu_table.xls")
        #         middle_file = open(middle_otu_file_path, "r")
        #         content2 = middle_file.readlines()
        #         with open(new_otu_file_path, "a") as w:
        #             w.write(first_line + "\n")

        #         t = 0
        #         for i in range(species_nu):  # print(list1[i])
        #             if list1[i] in list2:  # print(list1[i])
        #                 continue
        #             for line in content2:
        #                 line = line.strip("\n")
        #                 arr_2 = line.strip().split(":")
        #                 if arr_2[1] == str(list1[i]):
        #                     f = arr_2[0].strip().split(" ")
        #                     linecontent_2 = f[-1] + "\n"
        #                     # linecontent_2 = arr_2[0] + "\n"
        #                     list2.append(list1[i])  # 为重复的sum判断做准备
        #                     with open(new_otu_file_path, "a") as w:
        #                         w.write(linecontent_2)
        #                         t = t + 1
        #                         if t == species_nu:
        #                             break
        if self.option("method") == "":
            if self.option("sample_method") == "":
                self.set_db()
            else:
                self.run_sample_matrix()
        else:
            self.run_matrix()

    def transposition(self,old, path):
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
        #self.new_otu_file_path = os.path.join(self.work_dir, "new_otu_table.xls")
        self.transposition(new_otu_file_path, trans_otu)
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
            self.hcluster.on('end', self.set_db)
        self.hcluster.run()

    def run_sample_matrix(self): #20161206
        self.logger.info("正在进行样本距离计算")
        options = {
            "method": "bray_curtis",
            "otutable": new_otu_file_path
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
        self.sample_hcluster.on('end', self.set_db)
        self.logger.info("样本聚类计算结束")
        self.sample_hcluster.run()

    def sort_new_otu_file_path(self):     # by guanqing.zou 20180403
        sortfile = os.path.join(self.work_dir, "sort_otu_table.xls")
        self.logger.info("\n\n\n"+new_otu_file_path)
        s1 = open(new_otu_file_path, 'r')
        contento = s1.readlines()
        retlist=[]
        for n in contento[1:]:
            n = n.strip("\n")
            brr = n.strip().split("\t")
            sum=0
            for i in range(1,len(brr)):
                sum+=float(brr[i])
            retlist.append([n,sum])
        retlist_s=sorted(retlist,key=lambda b:b[1],reverse=True)
        with open(sortfile, "a") as w:
            w.write(contento[0])
            for oneline,t in retlist_s:
                w.write(oneline + "\n")
        s1.close()
        global new_otu_file_path
        new_otu_file_path = sortfile

    def replace_0_new_otu_file_path(self, infile):
        """
        将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一
        :return:
        """
        self.logger.info("正在将结果文件中的丰度值0替换")
        min_list = []   #add modify by qingchen.zhang@20190306 用于增加
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
                        line[i] = float(line[i]) + float(min_table/10)
                    else:
                        line[i] = float(line[i])
                    data_list.append(line[i])
                w.write("{}\t{}\n".format(line[0], "\t".join(str(i) for i in data_list)))
        # global new_otu_file_path
        # new_otu_file_path = real_zero_otu_new
        os.rename(real_zero_otu_new, infile)

    def set_db(self):
        # self.sort_new_otu_file_path()   #by guanqing.zou 20180403
        self.replace_0_new_otu_file_path(new_otu_file_path)  #by qingchen.zhang@20190307
        self.replace_0_new_otu_file_path(self.otu_relative)
        sample_tree = ""
        species_tree = ""
        sample_list = []
        species_list = []
        self.logger.info("正在写入mongo数据库")
        # myParams = json.loads(self.sheet.params)
        if self.option("method") != "":
            species_tree_path = self.hcluster.option("newicktree").prop['path']
            if os.path.exists(species_tree_path):
                with open(species_tree_path, "r") as f:
                    species_tree = f.readline().strip()
                    raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                    species_list = [i[1] for i in raw_samp]
        if self.option("sample_method") != "":
            sample_tree_path = self.sample_hcluster.option("newicktree").prop['path']
            if os.path.exists(sample_tree_path):
                with open(sample_tree_path, "r") as f:
                    sample_tree = f.readline().strip()
                    raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', sample_tree)
                    sample_list = [i[1] for i in raw_samp]
        # else:
            # if self.option("add_Algorithm") != "":
        if (self.option("add_Algorithm") != "" and self.option("sample_method") == ""):
            sample_name = open(self.group_table_path, "r")
            content = sample_name.readlines()
            for f in content:
                f = f.strip("\n")
                arr = f.strip().split("\t")
                if arr[0] != "#sample":
                    if arr[1] not in sample_list:
                        sample_list.append(arr[1])
        api_otu = self.api.hierarchical_clustering_heatmap
        # new_otu_id = api_otu.add_sg_hc_heatmap(self.sheet.params, self.option("input_otu_id"), None,
        #                                        sample_tree = sample_tree, sample_list = sample_list,
        #                                        species_tree = species_tree, species_list = species_list)
        api_otu.add_sg_hc_heatmap_detail(new_otu_file_path, self.color_dict, self.option("main_id"), self.option("input_otu_id"),
                                         sample_tree=sample_tree, sample_list=sample_list,
                                         species_tree = species_tree, species_list = species_list,
                                         otu_relative=self.otu_relative)  # add color_dict by houshuang 20191010
        # self.add_return_mongo_id("sg_hc_heatmap", self.option("main_id"))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"),"sg_hc_heatmap")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "otu_analysis_heatmap",
                "interaction": 1,
                "main_table": "sg_hc_heatmap",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        # shutil.copy(os.path.join(self.work_dir, "sort_otu_table.xls"), self.output_dir + "/heatmap.taxa.table.xls")
        os.link(new_otu_file_path, os.path.join(self.output_dir, "heatmap.taxa.table.xls"))
        if os.path.exists(self.sample_hcluster.output_dir + "/hcluster.tre"):
            shutil.copy(self.sample_hcluster.output_dir + "/hcluster.tre", self.output_dir + "/sample_hcluster.tre")
        if os.path.exists(self.hcluster.output_dir + "/hcluster.tre"):
            shutil.copy(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + "/species_hcluster.tre")
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "群落Heatmap分析结果输出目录", 0, "110077"],
            ["./heatmap.taxa.table.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, "110079"],
            ["./heatmap.taxa.relative.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, ""],
            ["./sample_hcluster.tre", "tre", "样本聚类树", 0, "110078"],
            ["./species_hcluster.tre", "tre", "物种聚类树", 0, "110080"],
            ["./群落Heatmap图.pdf", "pdf", "物种群落Heatmap图", 0, ""]
        ])
        super(HierarchicalClusteringHeatmapWorkflow, self).end()

    def run(self):
        self.run_sort_samples()
        super(HierarchicalClusteringHeatmapWorkflow, self).run()
