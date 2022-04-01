# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import re
from bson.objectid import ObjectId
import pandas as pd
import numpy as np
import math
import json

class GenesetClusterWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "express_file", "type": "string", "default": "none"},  # 输入文件，差异基因表达量矩阵
            {"name": "samples_distance_method", "type": "string","default":"complete"},  # 计算距离的方式
            {"name": "genes_distance_method", "type": "string","defalut":"complete"},  # 计算距离的方式
            {"name": "samples_distance_algorithm", "type": "string", "default": "pearson"},# 基因距离算法，只对hclust，默认是euclidean
            {"name": "genes_distance_algorithm", "type": "string", "default": "euclidean"},# 样本聚类算法，只对hclust, 默认是pearson
            {"name": "log", "type": "int", "default": None},  # 画热图时对原始表进行取对数处理，底数为10或2
            {"name": "method", "type": "string", "default": "hclust"},  # 聚类方法选择
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "type","type":"string","default":"gene"}, #gene/transcript 给to_file传递参数
            {"name": "express_method","type":"string","default":"rsem"},#rsem/featurecounts 给to_file传递参数参数
            {"name": "sub_num", "type": "int", "default": 0},  # 子聚类的数目
            {"name": "geneset_cluster_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name":"class_code_type","type":"string","default":"express_diff"},
            # {"name":"gene_cluster",'type':bool}, #是否基因聚类
            # {"name":"sample_cluster",'type':bool}, #是否样本聚类
            {"name":"class_code","type":"string"}, #class_code信息
            {"name": "gene_list", "type": "string"},  #输出gene_list
            {"name": "use_group", "type": "string"},  # added by gdq
            {"name": "express_level", "type": "string"}# 对应 tpm/fpkm字段, 传给workflow
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cluster = self.add_tool("rna.cluster")
        self.output_dir = self.cluster.output_dir
        with open(self.option('group_id'), 'r+') as f1:
            f1.readline()
            if not f1.readline():
                self.group_id = 'all'
            else:
                self.group_id = self.option('group_id')

    def get_samples(self): #add by khl 20170504
        edger_group_path = self.option("group_id")
        self.logger.info(edger_group_path)
        samples = []
        with open(edger_group_path,'r+') as f1:
            f1.readline()
            for lines in f1:
                line=lines.strip().split("\t")
                samples.append(line[0])
        print(samples)
        return samples

    def get_group_dict(self):
        group_file = self.option("group_id")
        group_dict = dict()
        with open(group_file) as f:
            _ = f.readline()
            for line in f:
                sample, group = line.strip().split()
                group_dict.setdefault(group, list())
                group_dict[group].append(sample)
        return group_dict
    
    def fpkm(self,samples):  #add by khl 20170504
        fpkm_path = self.option("express_file").split(",")[0]
        if self.group_id.lower() == "all":
            return fpkm_path
        fpkm = pd.read_table(fpkm_path, sep="\t")
        print(fpkm.columns)
        no_samp = []
        sample_total = fpkm.columns[1:]
        for sam in sample_total:
            try:
                if sam not in samples:
                    no_samp.append(sam)
            except Exception:
                pass
        if no_samp:
            new_fpkm = fpkm.drop(no_samp, axis=1)
            # print new_fpkm.columns
            self.new_fpkm = self.cluster.work_dir + "/ref_fpkm"
            header=['']
            header.extend(samples)
            new_fpkm.columns=header
            new_fpkm.to_csv(self.new_fpkm, sep="\t",index=False)
            return self.new_fpkm
        else:
            return fpkm_path
    
    def run_cluster(self):
        if self.group_id in ['all', 'All', 'ALL']:
            new_fpkm = self.option("express_file").split(",")[0]
        else:
            specimen = self.get_samples()
            new_fpkm = self.fpkm(specimen)
        if self.option('use_group') == 'yes' and (self.group_id not in ['all', 'All', 'ALL']):
            group_dict = self.get_group_dict()
            exp_pd = self.process_exp_matrix(new_fpkm, log_base=None, group_dict=group_dict)
            exp_pd.to_csv(new_fpkm, header=True, index=True, sep='\t')
        self.logger.info(self.option("method"))
        options = {
            "sub_num": self.option("sub_num"),
            "method": self.option("method"),
            "log": self.option("log")
        }
        if self.option("genes_distance_method") == "":
            genes_distance_method ="complete"
            options["genes_distance_method"]=genes_distance_method
        if self.option("method") == "hclust":
            options['genes_distance_algorithm'] = self.option("genes_distance_algorithm")
            options['samples_distance_algorithm'] = self.option("samples_distance_algorithm")
            if self.option("samples_distance_method") == "":
                samples_distance_method = "complete"
                options["samples_distance_method"]=samples_distance_method

        gene_list = self.option("gene_list")
        self.logger.info("gene_list文件路径为:{}".format(self.option("gene_list")))
        gene_list_id = []
        with open(gene_list,'r+') as f1:
            for lines in f1:
                gene_list_id.append(lines.strip())
        self.logger.info("打印gene_list_id")
        # self.logger.info(gene_list_id)
        out_fpkm_path = self.cluster.work_dir + "/new_fpkm"
        self.filter_file(new_fpkm, gene_list_id, out_fpkm_path)  #过滤掉samples 和gene_list 剩下的表达量矩阵做聚类分析
        options["diff_fpkm"] = out_fpkm_path
        self.logger.info("生成fpkm文件成功!")
        self.logger.info(out_fpkm_path)
        
        self.cluster.set_options(options)
        self.cluster.on("end", self.set_db)
        self.cluster.run()

    def filter_file(self, infile, seq_list, output):
        with open(infile, 'rb') as r, open(output, 'wb') as w:
            w.write(r.readline())
            self.logger.info("开始筛选基因id")
            i = 0
            for line in r:
                # self.logger.info(line.split('\t')[0])
                if line.split('\t')[0] in seq_list:
                    i += 1
                    w.write(line)
            if i <= 1:
                self.logger.error("基因集中的基因/转录本均不属于选择的表达量表")

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_cluster = self.api.denovo_cluster  # #不确定,增加一个database 
        if self.option("method") == "hclust":
            self.logger.info("开始导mongo表！")
            hclust_path = os.path.join(self.output_dir, "hclust")
            sub_clusters = os.listdir(hclust_path)
            with open(self.cluster.work_dir + '/hc_gene_order') as r:
                genes = [i.strip('\n') for i in r.readlines()]
            with open(self.cluster.work_dir + '/hc_sample_order') as r:
                specimen = [i.strip('\n') for i in r.readlines()]
            for sub_cluster in sub_clusters:
                # if re.search(r'hclust_heatmap.xls', sub_cluster):  # 在输出文件夹里找到hclust文件夹里面的热图表
                #     sub_cluster_path = os.path.join(hclust_path, sub_cluster)
                #     with open(sub_cluster_path, 'rb') as heatmap:
                #         specimen = heatmap.readline()  # 第一行信息
                #         heat_lst = heatmap.readlines()[1:]
                #         for gene_num in heat_lst:
                #             gene = gene_num.split('\t')[0]
                #             genes.append(gene)  # 获取第一列从第二行开始的信息，返回列表

                if re.match('subcluster', sub_cluster):  # 找到子聚类的文件进行迭代
                    sub = sub_cluster.split("_")[1]
                    sub_path = os.path.join(hclust_path, sub_cluster)
                    api_cluster.add_cluster_detail(cluster_id=self.option("geneset_cluster_id"), sub=sub, class_code =self.option('class_code'),sub_path=sub_path,project='ref',workflow=False,query_type=self.option("type"))
                    self.logger.info("开始导子聚类函数！")
                if re.search('samples_tree', sub_cluster):  # 找到sample_tree
                    sample_tree = os.path.join(hclust_path, sub_cluster)
                    self.logger.info("sample_tree产生")
                    
                if re.search('genes_tree', sub_cluster):  # 找到gene_tree
                    gene_tree = os.path.join(hclust_path, sub_cluster)
                    self.logger.info("gene_tree产生")
                    
            self.update_cluster(table_id=self.option("geneset_cluster_id"), genes=genes, sample_tree=sample_tree, gene_tree=gene_tree, specimen=specimen)
            self.logger.info(gene_tree)
            self.logger.info(genes)
            self.logger.info(specimen)
            self.logger.info("更新geneset_cluster主表成功！")
        else:
            kmeans_path = os.path.join(self.output_dir, "kmeans")
            sub_clusters = os.listdir(kmeans_path)
            genes = []
            for sub_cluster in sub_clusters:
                if re.match('subcluster', sub_cluster):
                    sub = sub_cluster.split("_")[1]
                    sub_path = os.path.join(kmeans_path, sub_cluster)
                    api_cluster.add_cluster_detail(cluster_id=self.option("geneset_cluster_id"), sub=sub, class_code =self.option('class_code'), sub_path=sub_path,workflow=False,query_type=self.option("type"))
                if re.search(r'kmeans_heatmap.xls', sub_cluster):
                    sub_cluster_path = os.path.join(kmeans_path, sub_cluster)
                    with open(sub_cluster_path, 'rb') as heatmap:
                        specimen = heatmap.readline().strip().split("\t")  # 第一行信息
                        
                        heat_lst = heatmap.readlines()
                        for gene_num in heat_lst:
                            gene = gene_num.split('\t')[0]
                            genes.append(gene)  # 获取
            self.update_cluster(table_id=self.option("geneset_cluster_id"), genes=genes, specimen=specimen, gene_tree=None, sample_tree=None)
        self.end()

    def update_cluster(self, table_id, genes, specimen, gene_tree, sample_tree):
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_ref_rna'
        collection = db['sg_geneset_cluster']
        if gene_tree:
            with open(gene_tree, 'rb') as g:
                gene_tree = g.readlines()[0].strip('\n')
        if sample_tree:
            with open(sample_tree, 'rb') as s:
                sample_tree = s.readlines()[0].strip('\n')
        collection.update({'_id': ObjectId(table_id)}, {'$set': {'specimen': specimen, 'genes': genes, 'gene_tree': gene_tree, 'sample_tree': sample_tree}})

    def run(self):
        self.run_cluster()
        super(GenesetClusterWorkflow, self).run()

    def end(self):
        output1_dir = self.output_dir
        result = self.add_upload_dir(output1_dir)
        result.add_relpath_rules([[".", "", "基因集聚类分析结果文件"], ])
        super(GenesetClusterWorkflow, self).end()

    def process_exp_matrix(self, exp_matrix, log_base=None, group_dict=None):
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
        else:
            print(exp_matrix, 'is assumed to be a pandas DataFrame Object')
            all_exp_pd = exp_matrix
        all_exp_pd.index.name = 'seq_id'

        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)

        if log_base:
            if log_base == math.e:
                all_exp_pd = np.log(all_exp_pd + 1)
            elif log_base == 2:
                all_exp_pd = np.log2(all_exp_pd + 1)
            elif log_base == 10:
                all_exp_pd = np.log10(all_exp_pd + 1)
            else:
                raise Exception('log base of {} is not supported'.format(log_base))
        # return
        return all_exp_pd