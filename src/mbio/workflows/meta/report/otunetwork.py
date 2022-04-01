# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
from biocluster.workflow import Workflow
import os
from collections import defaultdict
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class OtunetworkWorkflow(Workflow):
    """
    报告中进行OTU网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OtunetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "network_id", "type": "string"},
            {"name": "add_Algorithm", "type": "string", "default": ""} # 分组样本求和算法，默认不求和
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.otunetwork = self.add_tool('meta.otu.otunetwork')

    def change_otuname(self, tablepath):
        """
        这一步骤只是将otu名字中的空格去掉
        :param tablepath:
        :return:
        """
        newtable = os.path.join(self.work_dir, 'otutable1.xls')
        f2 = open(newtable, 'w+')
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
                    line[0] = line_he
                    for i in range(0, len(line)):
                        if i == len(line) - 1:
                            f2.write("%s\n" % (line[i]))
                        else:
                            f2.write("%s\t" % (line[i]))
        f2.close()
        return newtable

    def cat_samples(self, otu, method, grouptable):
        """
        合并同一分组的样本，可进行求和，求平均，求中位数
        :param method:
        :return:
        """
        # grouptable = os.path.join(self.work_dir, "grouptable_input.group.xls")
        cat_otu_path = os.path.join(self.work_dir, "cat_otu.xls")
        sample_group = dict()  # 一个样本是属于哪个group的
        index_sample = dict()  # 一个OTU表中第几列属于哪个样本
        group_sample_num = defaultdict(int)  # 一个分组里面有多少的样本
        with open(grouptable, "rb") as r:
            line = r.next()
            for line in r:
                line = line.rstrip().split("\t")
                sample_group[line[0]] = line[1]
                group_sample_num[line[1]] += 1

        with open(otu, "rb") as r, open(cat_otu_path, 'wb') as w:
            group_list = list()
            for v in sample_group.values():
                group_list.append(v)
                group_list = list(set(group_list))
            # print group_list
            # l = len(group_list) #zx

            line = r.next().rstrip().split("\t")
            # print line
            for i in range(len(line)):
                index_sample[i] = line[i]
            # print index_sample
            w.write(index_sample[0] + "\t")
            w.write("\t".join(group_list) + "\n")
            for line in r:
                line = line.rstrip().split("\t")
                num = defaultdict(int)
                middle_num = defaultdict(int)
                tmp = list()
                list1 = []
                mid_num = dict()
                w.write(line[0] + "\t")
                for i in range(1, len(line)):
                    num[sample_group[index_sample[i]]] += int(line[i])
                for m in group_list:
                    for i in range(1, len(line)):
                        if sample_group[index_sample[i]] == m:
                            list1.append(int(line[i]))
                            if len(list1) == group_sample_num[m]:
                                list1.sort()
                                yu = int(group_sample_num[m]) % 2
                                index = int(int(group_sample_num[m]) / 2)
                                if yu == 0:
                                    mid_num[m] = int(round((int(list1[index - 1]) + int(list1[index])) / 2))
                                    list1 = []
                                else:
                                    mid_num[m] = list1[index]
                                    list1 = []
                if method == "sum":
                    for g in group_list:
                        tmp.append(str(num[g]))
                if method == "average":
                    for g in group_list:
                        avg = round(round(num[g],5) / (group_sample_num[g]), 5)
                        tmp.append(str(avg))
                if method == "middle":
                    for g in group_list:
                        tmp.append(str(mid_num[g]))
                w.write("\t".join(tmp))
                w.write("\n")
        return cat_otu_path

    def set_group_file(self):
        """
        tofile文件当group_id为all时，grouptable_input.group.xls为空，直接在这里处理all的分组情况
        :return:
        """
        grouptable = os.path.join(self.work_dir, "grouptable_input.group.xls")
        newtable = os.path.join(self.work_dir, 'otutable.xls')
        with open(newtable, "rb") as r, open(grouptable, "wb") as w:
            data1 = r.readlines()[0]
            w.write("#sample\tgroup" + "\n")
            data = data1.strip().split("\t")
            for i in range(1, len(data)):
                w.write(data[i] + "\tAll\n")
        return grouptable


    def run_otunetwork(self):
        self.logger.info(self.option('grouptable'))
        if self.option('add_Algorithm') in ["sum", "average", "middle"]:
            if self.option("group_id") == "all":
                grouptable = self.set_group_file()
            else:
                grouptable = os.path.join(self.work_dir, "grouptable_input.group.xls")
            otu = os.path.join(self.work_dir, 'otutable.xls')
            otu_hebing = self.cat_samples(otu, self.option('add_Algorithm'), grouptable)
            newtable = self.change_otuname(otu_hebing)
            options = {'otutable': newtable}
        else:
            newtable = self.change_otuname(self.option('otutable').prop['path'])
            if self.option("group_id") == 'all':
                options = {'otutable': newtable}
            else:
                options = {
                    'otutable': newtable,
                    'grouptable': self.option('grouptable')
                }
        self.otunetwork.set_options(options)
        self.otunetwork.on('end', self.set_db)
        self.output_dir = self.otunetwork.output_dir
        self.otunetwork.run()


    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "共现性网络分析结果目录", 0, "110184"],
            ["./real_node_table.txt", "txt", "网络节点属性表", 0, "110190"],
            ["./real_edge_table.txt", "txt", "网络边的属性表", 0,"110186"],
            ["./real_dc_otu_degree.txt", "txt", "网络物种节点度分布表", 0, "110187"],
            ["./real_dc_sample_degree.txt", "txt", "网络sample节点度分布表", 0, "110189"],
            ["./real_dc_sample_otu_degree.txt", "txt", "网络所有节点度分布表", 0, "110192"],
            ["./network_centrality.txt", "txt", "网络中心系数表", 0, "110188"],
            ["./network_attributes.txt", "txt", "网络单值属性表", 0, "110191"],
            ["./network_degree.txt", "txt", "OTU网络度统计总表", 0, "110185"],
            ["./共线性网络图.pdf", "pdf", "物种在样本中的共存在关系结果图", 0, ""]
        ])
        super(OtunetworkWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_otunetwork = self.api.otunetwork
        node_table_path = self.output_dir + '/real_node_table.txt'
        edge_table_path = self.output_dir + '/real_edge_table.txt'
        otu_degree_path = self.output_dir + '/real_dc_otu_degree.txt'
        sample_degree_path = self.output_dir + '/real_dc_sample_degree.txt'
        sample_otu_degree_path = self.output_dir + '/real_dc_sample_otu_degree.txt'
        network_centrality_path = self.output_dir + '/network_centrality.txt'
        network_attributes_path = self.output_dir + '/network_attributes.txt'
        network_degree_path = self.output_dir + '/network_degree.txt'
        if not os.path.isfile(node_table_path):
            self.logger.error("找不到报告文件:{}".format(node_table_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(edge_table_path):
            self.logger.error("找不到报告文件:{}".format(edge_table_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(otu_degree_path):
            self.logger.error("找不到报告文件:{}".format(otu_degree_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(sample_degree_path):
            self.logger.error("找不到报告文件:{}".format(sample_degree_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(sample_otu_degree_path):
            self.logger.error("找不到报告文件:{}".format(sample_otu_degree_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(network_centrality_path):
            self.logger.error("找不到报告文件:{}".format(network_centrality_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(network_attributes_path):
            self.logger.error("找不到报告文件:{}".format(network_attributes_path))
            self.set_error("找不到报告文件", code="12703101")
        if not os.path.isfile(network_degree_path):
            self.logger.error("找不到报告文件:{}".format(network_degree_path))
            self.set_error("找不到报告文件", code="12703101")

        api_otunetwork.add_node_table(file_path=node_table_path, table_id=self.option("network_id"))
        # api_otunetwork.add_node_table_group(file_path=node_table_path, table_id=self.option("network_id"))
        api_otunetwork.add_edge_table(file_path=edge_table_path, table_id=self.option("network_id"))
        api_otunetwork.add_network_attributes(file_path=network_attributes_path, table_id=self.option("network_id"))
        api_otunetwork.add_network_degree(file1_path=otu_degree_path, file2_path=sample_degree_path,
                                          file3_path=sample_otu_degree_path, table_id=self.option("network_id"))
        api_otunetwork.add_network_centrality(file_path=network_centrality_path, table_id=self.option("network_id"))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("network_id"), "sg_network")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("network_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "network_analysis",
                "interaction": 1,
                "main_table": "sg_network",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_otunetwork()
        super(OtunetworkWorkflow, self).run()
