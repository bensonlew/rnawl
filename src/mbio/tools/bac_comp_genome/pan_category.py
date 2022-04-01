# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,json
import subprocess
from biocluster.core.exceptions import OptionError
import shutil
import pandas as pd
import math
import re
from mbio.packages.bac_comp_genome.common_function import format_category,add_signature


class PanCategoryAgent(Agent):
    """
    version v1.0
    author: qingchen.zhang
    last modified:2019.10.16
    细菌比较基因组用其根据选择的分类方案统计category信息
    整数型  {core : {min: 20, max: 20}, dis:{min :2 , max : 19}, unique : {min:1, max:1}}
    百分比 {core : {min: 0.95, max: 1}, dis:{min :0.05 , max : 0.95}, unique : {min:0, max:0.05}}
    两种类型不同的处理，否则会逻辑混乱
    """

    def __init__(self, parent):
        super(PanCategoryAgent, self).__init__(parent)
        options = [
            {"name": "cluster", "type": "infile", "format": "sequence.profile_table"},  # 输入fasta文件
            {"name": "category", "type": "string"}, #输入分类方案类型{'core'：{"min":1, "max": 1}, 'dispensable'：{"min":0.15, "max": 0.95}, 'soft_core'{"min":0.95, "max": 1}, 'unique'{"min":0, "max": 0.15}}min能取到，max取不到
            {"name": "percent", "type":"bool", "default": True}, #false 表示按个数，true表示按百分比
            {"name": "category_name", "type":"string"}, #输入分组方案的名称
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},#输入group表选择哪些样本进行导表和计算
            # {"name": "step", "type":"string", "default": "1"}, #输入步长
        ]
        self.add_option(options)
        self.step.add_steps('cluster')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cluster.start()
        self.step.update()

    def step_end(self):
        self.step.cluster.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("cluster").is_set:
            raise OptionError("必须设置参数cluster")
        # if not self.option("percent"):
        #     raise OptionError("必须指定方案是按照百分比计算还是个数计算")
        if not self.option("category") :
            raise OptionError("必须设置分组方案")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PanCategoryAgent, self).end()


class PanCategoryTool(Tool):
    def __init__(self, config):
        super(PanCategoryTool, self).__init__(config)
        self._version = '1.0'
        self.cluster = self.option("cluster").prop["path"]

    def run(self):
        """
        开始运行
        :return:
        """
        self.logger.info('开始运行')
        super(PanCategoryTool, self).run()
        self.run_category()
        self.set_output()
        self.end()

    def run_category(self):
        """
        根据cluster聚类结果统计不同分组方案的list
        :return:
        """
        cluster_data = pd.read_table(self.cluster, sep="\t", header=0)
        columns = cluster_data.columns
        if self.option("group_table").is_set:
            all_sample_list = self.get_group()
            columns_number = len(all_sample_list) #得到样本数目
        else:
            all_sample_list = list(columns[3:])
            columns_number = len(all_sample_list) #得到样本数目
        all_sample_list.sort()
        for sample in all_sample_list:
            if sample not in columns:
                self.set_error("传入样本名称错误")
        try:
            all_category = json.dumps(self.option("category"))
            all_category = json.loads(all_category)
            all_category = eval(all_category)
        except Exception, e:
            self.logger.info("传入的category不是正确的json格式")
        category_list = self.option("category_name").split(",")
        category_list.sort()
        self.logger.info("all_category : {}".format(all_category))
        self.category = format_category(self.option("percent"),all_category, columns_number, category_list)

        self.logger.info(">>>>>>>>>>>>>%s"%self.category)
        # pop_data = cluster_data['Cluster_ID']
        cluster_data["Category"] = cluster_data['Cluster_ID']
        categ_index = {}
        self.logger.info(">>>>>>>>>>>>>%s"%columns_number)

        cluster_data = add_signature(cluster_data,self.category, category_list,self.option("percent"),columns_number)

        sample_dict = {}
        self.logger.info(cluster_data.head())
        cluster_data.to_csv("table.xls", sep='\t',index=0)
        cluster_table = cluster_data.set_index('Cluster_ID')
        sample_cluster_dict = {}
        for sample in all_sample_list: ##将样本的基因替换为cluster，方便统计cluster
            data_cluster = cluster_table[[sample, 'Category']]
            data_cluster = data_cluster[data_cluster[sample] != "-"]
            data_cluster[sample] = data_cluster[sample].index
            # self.logger.info("bug1%s完事了！"% data_cluster)
            data_cluster.reset_index()
            # self.logger.info("bug2%s完事了！"% data_cluster)
            #data_cluster = data_cluster.drop('Cluster_ID', axis=1)  ##gaohao
            sample_cluster = data_cluster.groupby(data_cluster['Category']).apply(lambda x:[';'.join(x[sample])])###gaohao
            self.logger.info("++++++++++++++++++++++%s--------------------"%sample_cluster.head())
            new_cluster_dict = {}
            for x in sample_cluster.index:
                new_cluster_dict[x] = ";".join(sample_cluster[x]).split(";")
            sample_cluster_dict[sample] = new_cluster_dict

        total_cluster_dict = {}
        for sample in all_sample_list: #获取所有样本总的cluster的list
            new_data_dict = sample_cluster_dict[sample]
            all_sample_cluster_list = []
            for cate in category_list:
                if cate in new_data_dict:
                    all_sample_cluster_list += new_data_dict[cate]
            all_sample_cluster_list = list(set(all_sample_cluster_list))
            # for data_j in new_data_list:
            #     if data_j != "-":
            #         if data_j not in all_sample_cluster_list:
            #             all_sample_cluster_list.append(data_j)
            total_cluster_dict[sample] = all_sample_cluster_list

        cluster_data = cluster_data.set_index('Category')
        for sample in all_sample_list:#获取所有样本的不同分类水平的list，这里的index为category，根据category找对应的gene
            newdata = cluster_data[sample]
            newdata = newdata[newdata != "-"]
            file_total = newdata.groupby(newdata.index).apply(';'.join)
            new_class_dict = {}
            for x in file_total.index:
                new_class_dict[x] = file_total[x]
            sample_dict[sample] = new_class_dict
        # self.logger.info(sample_dict["GCF_000216375.1"])

        total_gene_dict = {}
        for sample in all_sample_list: #获取所有每个样本的total_gene list
            new_level_dict = sample_dict[sample]
            all_sample_gene_list = []
            for cate in category_list:
                if cate in new_level_dict:
                    data_list = new_level_dict[cate].split(";")
                    for x in data_list:
                        all_sample_gene_list += [y.split("|")[1] for y in x.split(",")]
            # for data_i in new_data_list:
            #     if data_i != "-":
            #         data_gene_name_list = data_i.split(",")
            #         for g_n_j in data_gene_name_list:
            #             if g_n_j not in all_sample_gene_list:
            #                 all_sample_gene_list.append(g_n_j.split("|")[1])
            total_gene_dict[sample] = all_sample_gene_list

        stat_cluster = self.work_dir + "/pangenome_genes.xls" #获得各分类单元的数据和统计基因数目
        cluster_result = self.work_dir + "/pangenome_clusters.xls" #获得各分类单元的数据和统计cluster数目
        with open(stat_cluster, 'w') as w, open(cluster_result, "w") as outf:
            category_list_number = [(x + "_num") for x in category_list]
            total_categ_list = list(set(category_list).union(set(category_list_number)))
            total_categ_list.sort()
            w.write("Sample_name\tTotal_Gene\tTotal_num\t{}\n".format('\t'.join(total_categ_list)))
            outf.write("Sample_name\tTotal_Cluster\tTotal_num\t{}\n".format('\t'.join(total_categ_list)))
            for sample in all_sample_list:
                sample_name = sample
                total = len(total_gene_dict[sample])
                # total_list = []
                # for i in total_gene_dict[sample]:
                #     if i != '-':
                #         sample_gene_list = i.split(",")
                #         for j in sample_gene_list:
                #             gene_name = j.split("|")[1]
                #             total_list.append(gene_name)
                w.write("{}\t{}\t{}".format(sample_name, ','.join(total_gene_dict[sample]), total))

                total_cluster = len(total_cluster_dict[sample])## 计算总的数目
                # total_cluster_list = []
                # for i in total_cluster_dict[sample]:## 计算总的cluster的数目
                #     if i not in total_cluster_list:
                #         total_cluster_list.append(i)
                outf.write("{}\t{}\t{}".format(sample_name, ','.join(total_cluster_dict[sample]), total_cluster))

                for categ in category_list:
                    # self.logger.info(sample_dict[sample].keys())
                    if categ in sample_dict[sample]:
                        categ_name = sample_dict[sample][categ]
                        self.logger.info("bug3%s完事了！"% categ_name)
                        categ_list = categ_name.split(";")
                        new_categ_list = []
                        for gene in categ_list:
                            if gene != '-':
                                new_gene_name = gene.split(",")
                                for gene_i in new_gene_name:
                                    name = gene_i.split("|")[1]
                                    if name not in new_categ_list:
                                        new_categ_list.append(name)
                        categ_num = len(new_categ_list)
                        categ_index_name = ','.join(new_categ_list)
                        # self.logger.info(categ_index_name)
                        w.write("\t{}\t{}".format(categ_index_name, categ_num))
                    else:
                        w.write("\t{}\t{}".format("-", 0))

                    sample_categ_dict = sample_cluster_dict[sample]
                    if categ in sample_categ_dict:
                        new_categ_list2 = sample_categ_dict[categ]
                        # self.logger.info("bug4%s完事了！"% categ_name2)
                        # categ_list2 = categ_name2[0]
                        # new_categ_list2 = []
                        # for cluster in categ_list2:
                        #     if cluster not in new_categ_list2:
                        #         new_categ_list2.append(cluster)
                        categ_num2 = len(new_categ_list2)
                        categ_index_name2 = ','.join(new_categ_list2)
                        outf.write("\t{}\t{}".format(categ_index_name2, categ_num2))
                    else:
                        outf.write("\t{}\t{}".format("-", 0))

                w.write('\n')
                outf.write('\n')

        distribution = self.work_dir + "/pangenome_distribution.xls" #获得分布直方图的数据
        with open(distribution, 'w') as outw:
            outw.write("Sample_num\tCluster_num\n")
            for i in range(1, len(all_sample_list)+1):
                sample_count_list = list(cluster_data['Sample_number'])
                sample_count = sample_count_list.count(i)
                outw.write("{}\t{}\n".format(i, sample_count))

    def get_group(self):
        """
        根据传入的group表获取字段名称
        :return:
        """
        all_list = []
        group_file = self.option("group_table")
        with open(group_file, 'r') as f:
            for line in f:
                line = line.strip().split("\t")
                group_name = line[0]
                sample_name = line[1]
                if sample_name not in all_list:
                    all_list.append(sample_name)
        return all_list

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始链接文件到结果文件夹！")
        oldfiles = os.path.join(self.work_dir, "pangenome_genes.xls")
        newfiles = os.path.join(self.output_dir, "pangenome_genes.xls")
        if os.path.exists(newfiles):
            os.remove(newfiles)
        os.link(oldfiles,newfiles)
        oldfiles2 = os.path.join(self.work_dir, "pangenome_clusters.xls")
        newfiles2 = os.path.join(self.output_dir, "pangenome_clusters.xls")
        if os.path.exists(newfiles2):
            os.remove(newfiles2)
        os.link(oldfiles2,newfiles2)

        pangenome_distribution = os.path.join(self.output_dir, "pangenome_distribution.xls")
        if os.path.exists(pangenome_distribution):
            os.remove(pangenome_distribution)
        os.link(os.path.join(self.work_dir, "pangenome_distribution.xls"), pangenome_distribution)
        self.logger.info("链接文件完成！")
