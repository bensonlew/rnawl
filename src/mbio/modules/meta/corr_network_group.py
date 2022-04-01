# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
from mainapp.models.mongo.public.meta.meta import Meta

class CorrNetworkGroupModule(Module):

    def __init__(self, work_id):
        super(CorrNetworkGroupModule, self).__init__(work_id)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "grouptable", "type": "string"}, # group 文件的目录
            {"name": "lable", "type": "float", "default": 0.03},  # 设定lable（距离水平）的值
            {"name": "method", "type": "string", "default": "pearson"},  # 设定计算相似性的算法
            {"name": "coefficient", "type": "float", "default": 0.08},  # 相关性系数阈值
            {"name": "significance", "type": "float", "default": 0.05},
            {"name": "color_level", "type": "string", "default": "3"}, # 颜色显示水平
            {"name": "abundance", "type": "int", "default": 50},  #设定物种总丰度值前50的物种信息
        ]
        self.add_option(options)
        self.cor_network_list = []

    def run_cor_network(self):
        newtable = self.change_otuname(self.option('otutable').prop['path'])
        for file in os.listdir(self.option("grouptable")):
            all_sample = []
            with open(self.option("grouptable") + "/" + file) as f:
                for i in f.readlines()[1:]:
                    all_sample.append(i.strip().split("\t")[0])
            if not len(all_sample) < 3:
                self.corr_network_analysis = self.add_module('meta.otu.corr_network_analysis')
                options = {
                    'otutable': newtable,
                    'grouptable': self.option("grouptable") + "/" + file,
                    'lable': self.option('lable'),
                    'method': self.option('method'),
                    'coefficient': self.option("coefficient"),
                    'significance': self.option("significance"),
                    'group_name': os.path.basename(file).strip(".group.txt"),
                }
                self.corr_network_analysis.set_options(options)
                self.cor_network_list.append(self.corr_network_analysis)
        if len(self.cor_network_list) > 1:
            self.on_rely(self.cor_network_list, self.set_output)
        else:
            self.cor_network_list[0].on('end', self.set_output)
        for tool in self.cor_network_list:
            tool.run()

    def change_otuname(self, tablepath):
        '''
        处理原始的otu表的otu名字过长，然后对otu表中的物种的总丰度值进行统计，然后筛选留下丰度值前x的物种，默认是50
        :param tablepath:
        :return:
        '''
        newtable = os.path.join(self.work_dir, 'otutable1.xls')  # 保存只留最后的物种名字文件
        newtable4 = os.path.join(self.work_dir, 'species_color.txt')  # 保存物种与门对应的数据，用于在后面的物种丰度中添加门数据
        x = self.option('abundance')
        if x == 0:
            x = 500
        else:
            pass
        f2 = open(newtable, 'w+')
        f4 = open(newtable4, 'w+')
        color_index = int(self.option('color_level')) - 1  ###
        with open(tablepath, 'r') as f:
            i = 0
            for line in f:
                if i == 0:
                    i = 1
                    f2.write(line)
                else:
                    line = line.strip().split('\t')
                    line_data = line[0].strip().split(' ')
                    line_he = "".join(line_data)
                    # tmp1 = line_he.strip().split(";")[-1:][0]  # 输出最后一个物种名
                    tmp2 = line_he.strip().split(";")[color_index]  # 输出门分类

                    line[0] = line_he
                    f4.write(line_he + '\t' + tmp2 + "\n")
                    for i in range(0, len(line)):
                        if i == len(line) - 1:
                            f2.write("%s\n" % (str(int(float(line[i])))))
                        elif i == 0:
                            f2.write("%s\t" % (line[i]))
                        else:
                            f2.write("%s\t" % (str(int(float(line[i])))))
        f2.close()
        f4.close()
        newtable1 = os.path.join(self.work_dir, 'species_abundance.txt')
        newtable3 = os.path.join(self.work_dir, 'otutable2.xls')
        f1 = open(newtable1, 'w+')  # 报存物种与总丰度的对应值，画网络图需要
        f3 = open(newtable3, 'w+')  # 报存物种总丰度前50的物种
        with open(newtable, "r") as r:
            data_line = r.readlines()[1:]
            list_sum = []
            list_otu = []
            new_dict = {}
            for line in data_line:
                line = line.strip().split("\t")
                list_otu.append(line[0])
                list1 = []
                for i in range(1, len(line)):
                    list1.append(int(float(line[i])))
                    a = sum(list1)
                list_sum.append(a)
            new_dict = dict(zip(list_otu, list_sum))
            new_dict1 = sorted(new_dict.iteritems(), key=lambda d: d[1], reverse=True)[0:int(x)]
            list_result_otu = []
            f1.write("species" + "\t" + "abundance" + "\t" + "color_level" + "\n")
            for i in new_dict1:
                if int(i[1]) != 0:
                    list_result_otu.append(str(i[0]))
                with open(newtable4, "r") as x:
                    data = x.readlines()
                    for line in data:
                        line = line.strip().split("\t")
                        if str(line[0]) == str(i[0]) and str(i[0]) in list_result_otu:
                            f1.write(str(i[0]).strip().split(';')[-1:][0] + "\t" + str(i[1]) + "\t" + str(
                                line[1]) + "\n")
        with open(newtable, "r") as m:
            i = 0
            for line1 in m:
                if i == 0:
                    i = 1
                    f3.write(line1)
                else:
                    line1 = line1.strip().split("\t")
                    if line1[0] in list_result_otu:
                        line1_data = line1[0].strip().split(';')[-1:][0]
                        line1[0] = line1_data
                        for i in range(0, len(line1)):
                            if i == len(line1) - 1:
                                f3.write("%s\n" % (line1[i]))
                            else:
                                f3.write("%s\t" % (line1[i]))
        f1.close()
        f3.close()
        return newtable3

    def set_output(self):
        if not os.path.exists(self.output_dir + "/CorrNetworkSpearmanGenus/"):
            os.mkdir(self.output_dir + "/CorrNetworkSpearmanGenus/")
        for tool in self.cor_network_list:
            if os.path.exists(tool.output_dir):
                try:
                    file_name = os.listdir(tool.output_dir)[0]
                    if (os.path.exists(tool.output_dir + '/' + file_name)) and (len(os.listdir(tool.output_dir + '/' + file_name)) > 0):
                        os.system('cp -r %s %s' % (tool.output_dir + '/' + file_name, self.output_dir + "/CorrNetworkSpearmanGenus/"))
                except:
                    self.logger.info("{} hava none result".format(tool.output_dir))
        os.link(self.work_dir + '/species_abundance.txt',self.output_dir + "/CorrNetworkSpearmanGenus/species_abundance.txt")
        os.link(self.option('otutable').prop["path"], self.output_dir  + '/CorrNetworkSpearmanGenus/newtable.txt')
        self.end()

    def run(self):
        super(CorrNetworkGroupModule, self).run()
        self.run_cor_network()


    def end(self):
        super(CorrNetworkGroupModule, self).end()

