# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
from biocluster.workflow import Workflow
import os
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class CorrNetworkWorkflow(Workflow):
    """
    报告中进行物种相关性网络构建与分析时使用
    多样性 单因素网络分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CorrNetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "level", "type": 'int', "default": 9},
            {"name": "lable", "type": "float", "default": 0.03},
            {"name": "method", "type": "string", "default": "pearson"},
            {"name": "coefficient", "type": "float", "default": 0.04},
            {"name": "significance", "type": "float", "default": 0.05},
            {"name": "abundance", "type": "int", "default": 50},  #设定物种总丰度值前50的物种信息
            {"name": "update_info", "type": "string"},
            {"name": "corr_network_id", "type": "string"},
            {"name": "color_level", "type": "string", "default": "4"},           # 颜色显示水平  v4
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.corr_network_analysis = self.add_module('meta.otu.corr_network_analysis')

    def change_otuname(self, tablepath):
        '''
        处理原始的otu表的otu名字过长，然后对otu表中的物种的总丰度值进行统计，然后筛选留下丰度值前x的物种，默认是50
        :param tablepath:
        :return:
        '''
        newtable = os.path.join(self.work_dir, 'otutable1.xls')   #保存只留最后的物种名字文件
        newtable4 = os.path.join(self.work_dir, 'species_color.txt')  #保存物种与门对应的数据，用于在后面的物种丰度中添加门数据
        x = self.option('abundance')
        if x == 0:
            x = 500
        else:
            pass
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


    def run_corrnetworkcalc(self):
        newtable = self.change_otuname(self.option('otutable').prop['path'])
        print newtable
        options = {
            'otutable': newtable,
            'grouptable': self.option('grouptable'),
            'lable': self.option('lable'),
            'method': self.option('method'),
            'coefficient': self.option("coefficient"),
            'significance': self.option("significance")
            }
        self.corr_network_analysis.set_options(options)
        self.corr_network_analysis.on('end', self.set_db)
        self.output_dir = self.corr_network_analysis.output_dir
        self.corr_network_analysis.run()


    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        repaths = [
            [".", "", "物种相关性网络分析结果目录", 0, "110193"],
            ["./association", "", "物种相关性计算结果输出目录", 0, "110201"],
            ["./corr_network_calc", "", "物种相关性网络分析结果输出目录", 0, "110194"],
            ["./association/shared.txt", "txt", "shared文件", 0, "110203"],
            ["./corr_network_calc/corr_network_attributes.txt", "txt", "网络的单值属性表", 0, "110197"],
            ["./corr_network_calc/corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件", 0, "110196"],
            ["./corr_network_calc/corr_network_centrality.txt", "txt", "网络节点的中心系数表", 0, "110195"],
            ["./corr_network_calc/corr_network_clustering.txt", "txt", "网络节点的聚类系数表", 0, "110200"],
            ["./corr_network_calc/corr_network_degree_distribution.txt", "txt", "网络节点的度分布表", 0, "110199"],
            ["./corr_network_calc/corr_network_node_degree.txt", "txt", "网络节点的度统计表", 0, "110198"],
            ["./corr_network_calc/单因素相关性网络图.pdf", "pdf", "物种与物种间的相关性网络图", 0, ""],
            ["./heatmap/物种相关性Heatmap图.pdf", "pdf", "物种与物种间的相关性Heatmap图", 0, ""]
        ]
        regexps = [
            [r"^association/shared", "corr", "物种相似性网络边文件", 0, "110202"]
        ]
        os.rename(self.output_dir+"/otu_association", self.output_dir+"/association")
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        for i in self.get_upload_files():
            self.logger.info('upload file:{}'.format(str(i)))
        super(CorrNetworkWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_corrnetwork = self.api.corr_network
        # lable = self.option("lable")
        # method = self.option("method")
        # corr_name = "shared." + str(lable) + "." + str(method) + ".otu.corr"
        # node_links_path = os.path.join(self.output_dir, "otu_association/", corr_name)
        # print node_links_path
        node_links_path = self.output_dir + '/corr_network_calc/corr_network_by_cut.txt'
        node_abundance_path = self.work_dir + '/species_abundance.txt'
        network_clustering_path = self.output_dir + '/corr_network_calc/corr_network_clustering.txt'
        network_degree_path = self.output_dir + '/corr_network_calc/corr_network_node_degree.txt'
        network_centrality_path = self.output_dir + '/corr_network_calc/corr_network_centrality.txt'
        network_attributes_path = self.output_dir + '/corr_network_calc/corr_network_attributes.txt'
        network_degree_distribution = self.output_dir + '/corr_network_calc/corr_network_degree_distribution.txt'
        if not os.path.isfile(node_links_path):
            self.set_error("没有符合相应条件的物种相关关系存在，请将相关系数阈值调小点，或者请选择足够多的样本与物种，"
                           "或者选择更低的分类水平！", code="12700901")
        if not os.path.isfile(node_abundance_path):
            self.logger.error("找不到报告文件:{}".format(node_abundance_path))
            self.set_error("找不到报告文件", code="12700902")
        if not os.path.isfile(network_clustering_path):
            self.logger.error("找不到报告文件:{}".format(network_clustering_path))
            self.set_error("找不到报告文件", code="12700902")
        if not os.path.isfile(network_degree_distribution):
            self.logger.error("找不到报告文件:{}".format(network_degree_distribution))
            self.set_error("找不到报告文件", code="12700902")
        if not os.path.isfile(network_centrality_path):
            self.logger.error("找不到报告文件:{}".format(network_centrality_path))
            self.set_error("找不到报告文件", code="12700902")
        if not os.path.isfile(network_attributes_path):
            self.logger.error("找不到报告文件:{}".format(network_attributes_path))
            self.set_error("找不到报告文件", code="12700902")
        if not os.path.isfile(network_degree_path):
            self.logger.error("找不到报告文件:{}".format(network_degree_path))
            self.set_error("找不到报告文件", code="12700902")

        api_corrnetwork.add_network_links_table(file_path=node_links_path, node_id_file=network_degree_path ,table_id=self.option("corr_network_id"))
        api_corrnetwork.add_network_abundance_table(file_path=node_abundance_path, table_id=self.option("corr_network_id"))
        api_corrnetwork.add_network_cluster_degree(file1_path=network_degree_path, file2_path=network_clustering_path, table_id=self.option("corr_network_id"))
        api_corrnetwork.add_network_centrality(file_path=network_centrality_path, table_id=self.option("corr_network_id"))
        api_corrnetwork.add_network_attributes(file_path=network_attributes_path, table_id=self.option("corr_network_id"))
        api_corrnetwork.add_network_degree_distribution(file_path=network_degree_distribution, table_id=self.option("corr_network_id"))
        corr_file =  self.output_dir + '/heatmap/corr.xls'
        p_file = self.output_dir + '/heatmap/pvalue.xls'
        tree_file = self.output_dir + '/heatmap/corr.cluster_tree.xls'
        if not os.path.exists(tree_file):
            tree_file = None
        api_corrnetwork.add_heatmap_corr_detail(self.option("corr_network_id"), corr_file , p_file, tree_file=tree_file)


        # corr_network_id_tab = self.option("corr_network_id")
        # self.add_return_mongo_id('sg_corr_network', corr_network_id_tab)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("corr_network_id"), "sg_corr_network")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("corr_network_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "corr_network_analysis",
                "interaction": 1,
                "main_table": "sg_corr_network",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_corrnetworkcalc()
        super(CorrNetworkWorkflow, self).run()
