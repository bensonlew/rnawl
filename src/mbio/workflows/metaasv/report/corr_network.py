# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.workflow import Workflow
import os
from mbio.packages.metaasv.common_function import link_dir


class CorrNetworkWorkflow(Workflow):
    """
    metaasv多样性 单因素相关性网络
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CorrNetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},##选择的otu表
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},##group表
            {"name": "group_detail", "type": "string"},##传入的group_detail
            {"name": "level", "type": 'int', "default": 9},##选择的分类水平
            {"name": "lable", "type": "float", "default": 0.03},##固定参数设置的label水平
            {"name": "method", "type": "string", "default": "pearson"},##相关性系数计算方法
            {"name": "coefficient", "type": "float", "default": 0.4},##此系数不在页面显示，相关性系数
            {"name": "significance", "type": "float", "default": 0.05},#此系数不在页面显示，显著性水平
            {"name": "abundance", "type": "int", "default": 50},  #设定物种总丰度值前50的物种信息
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},##主表ID
            {"name": "color_level", "type": "string", "default": "4"},# 颜色显示水平
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        # self.corr_network_analysis = self.add_module('meta.otu.corr_network_analysis')
        self.otu_association = self.add_tool("meta.otu.otu_association")
        self.corr_network_calc = self.add_tool("metaasv.corr_network_calc")

    def change_otuname(self, tablepath):
        '''
        处理原始的otu表的otu名字过长，然后对otu表中的物种的总丰度值进行统计，然后筛选留下丰度值前x的物种，默认是50
        :param tablepath:
        :return:
        '''
        newtable = os.path.join(self.work_dir, 'otutable1.xls')   #保存只留最后的物种名字文件
        newtable4 = os.path.join(self.work_dir, 'species_color.txt')  #保存物种与门对应的数据，用于在后面的物种丰度中添加门数据
        f2 = open(newtable, 'w+')
        f4 = open(newtable4, 'w+')
        color_index = int(self.option('color_level')) -1   ###
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
                    #tmp1 = line_he.strip().split(";")[-1:][0]  # 输出最后一个物种名
                    tmp2 = line_he.strip().split(";")[color_index]  # 输出门分类

                    line[0] = line_he
                    f4.write(line_he + '\t' + tmp2 + "\n")
                    for i in range(0, len(line)):
                        if i == len(line) - 1:
                            f2.write("%s\n" % (line[i]))
                        else:
                            f2.write("%s\t" % (line[i]))
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
                    list1.append(int(line[i]))
                    a = sum(list1)
                list_sum.append(a)
            new_dict = dict(zip(list_otu, list_sum))
            new_dict1 = sorted(new_dict.iteritems(), key=lambda d: d[1], reverse=True)
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
                            f1.write(str(i[0]).strip().split(';')[-1:][0] + "\t" + str(i[1]) + "\t" + str(line[1]) + "\n")
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

    def otuassociation_run(self):
        self.otu_association.set_options({
            "otutable": self.newtable,
            "lable": self.option("lable"),
            "method": self.option("method")
        })
        self.otu_association.on('end', self.corrnetworkcalc_run)
        self.otu_association.run()

    def corrnetworkcalc_run(self):
        lable = self.option("lable")
        method = self.option("method")
        corr_name = "shared." + str(lable) + "." + str(method) + ".corr"
        corr_table = os.path.join(self.otu_association.output_dir, corr_name)
        print corr_table
        self.corr_network_calc.set_options({
            "corr_table": corr_table,
            "coefficient": self.option("coefficient"),
            "significance": self.option("significance")
        })
        self.corr_network_calc.on('end', self.set_db)
        self.corr_network_calc.run()

    def end(self):
        repaths = [
            [".", "", "物种相关性网络分析结果目录", 0, ""],
            ["./association", "", "物种相关性计算结果输出目录", 0, ""],
            ["./corr_network_calc", "", "物种相关性网络分析结果输出目录", 0, ""],
            ["./association/shared.txt", "", "shared文件", 0, ""],
            ["./association/shared.0.03.spearman.corr", "", "物种相关性计算结果表", 0, ""],
            ["./corr_network_calc/corr_network_attributes.txt", "", "网络的单值属性表", 0, ""],
            ["./corr_network_calc/corr_network_by_cut.txt", "", "相关系数筛选后网络边文件", 0, ""],
            ["./corr_network_calc/corr_network_centrality.txt", "", "网络节点的中心系数表", 0, ""],
            ["./corr_network_calc/corr_network_clustering.txt", "", "网络节点的聚类系数表", 0, ""],
            ["./corr_network_calc/corr_network_degree_distribution.", "txt", "网络节点的度分布表", 0, ""],
            ["./corr_network_calc/corr_network_node_degree.txt", "", "网络节点的度统计表", 0, ""]
        ]
        # regexps = [
        #     [r"^otu_association/shared", "", "物种相似性网络边文件夹", 0, ""]
        # ]
        os.rename(self.output_dir+"/otu_association", self.output_dir+"/association")
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        # sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))
        super(CorrNetworkWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        link_dir(self.otu_association.output_dir, self.output_dir+"/otu_association")
        link_dir(self.corr_network_calc.output_dir, self.output_dir+"/corr_network_calc")
        self.logger.info("链接结果文件成功！")
        api_corrnetwork = self.api.api("metaasv.corr_network")
        node_links_path = self.output_dir + '/corr_network_calc/corr_network_by_cut.txt'
        node_abundance_path = self.work_dir + '/species_abundance.txt'
        network_clustering_path = self.output_dir + '/corr_network_calc/corr_network_clustering.txt'
        network_degree_path = self.output_dir + '/corr_network_calc/corr_network_node_degree.txt'
        network_centrality_path = self.output_dir + '/corr_network_calc/corr_network_centrality.txt'
        network_attributes_path = self.output_dir + '/corr_network_calc/corr_network_attributes.txt'
        network_degree_distribution = self.output_dir + '/corr_network_calc/corr_network_degree_distribution.txt'
        if not os.path.isfile(node_links_path):
            self.set_error("没有符合相应条件的物种相关关系存在，请将相关系数阈值调小点，或者请选择足够多的样本与物种，"
                           "或者选择更低的分类水平！")
        if not os.path.isfile(node_abundance_path):
            self.logger.error("找不到报告文件:{}".format(node_abundance_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(network_clustering_path):
            self.logger.error("找不到报告文件:{}".format(network_clustering_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(network_degree_distribution):
            self.logger.error("找不到报告文件:{}".format(network_degree_distribution))
            self.set_error("找不到报告文件")
        if not os.path.isfile(network_centrality_path):
            self.logger.error("找不到报告文件:{}".format(network_centrality_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(network_attributes_path):
            self.logger.error("找不到报告文件:{}".format(network_attributes_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(network_degree_path):
            self.logger.error("找不到报告文件:{}".format(network_degree_path))
            self.set_error("找不到报告文件")
        main_id = self.option("main_id")
        profile1 = self.newtable
        attributes_file = self.output_dir + "/corr_network_calc/corr_network_attributes.txt"
        api_corrnetwork.add_network_corfd(attributes_file, main=False, main_table_id=main_id)
        api_corrnetwork.add_network_corfd_link(main_id, self.output_dir + '/corr_network_calc')
        api_corrnetwork.add_network_corfd_degree(main_id, self.output_dir+ '/corr_network_calc')
        api_corrnetwork.add_network_corfd_node(main_id, self.output_dir+ '/corr_network_calc', profile1,
                                               int(self.option('color_level')))
        self.end()

    def run(self):
        self.newtable = self.change_otuname(self.option('otutable').prop['path'])
        self.otuassociation_run()
        super(CorrNetworkWorkflow, self).run()
