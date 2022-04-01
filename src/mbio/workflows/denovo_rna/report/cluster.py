# -*- coding: utf-8 -*-
# __author__ = "chenyanyan"
# last_modify:2016.10.12

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import re
from bson.objectid import ObjectId


class ClusterWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "express_file", "type": "string", "default": "none"},  # 输入文件，差异基因表达量矩阵
            {"name": "distance_method", "type": "string", "default": "euclidean"},  # 计算距离的算法
            {"name": "log", "type": "int", "default": 10},  # 画热图时对原始表进行取对数处理，底数为10或2
            {"name": "cluster_method", "type": "string", "default": "hclust"},  # 聚类方法选择
            {"name": "sub_num", "type": "int", "default": 0},  # 子聚类的数目
            {"name": "cluster_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "gene_list", "type": "string"}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cluster = self.add_tool("denovo_rna.express.cluster")
        self.output_dir = self.cluster.output_dir

    def run_cluster(self):
        self.logger.info(self.option("cluster_method"))
        options = {
            "distance_method": self.option("distance_method"),
            "sub_num": self.option("sub_num"),
            "method": self.option("cluster_method"),
            "log": self.option("log")
        }
        if self.option("gene_list") == 'all':
            options["diff_fpkm"] = self.option("express_file").split(',')[0]
        else:
            gene_list = self.option('gene_list').split(',')
            self.filter_file(self.option("express_file").split(',')[0], gene_list, './new_fpkm')
            options["diff_fpkm"] = './new_fpkm'
        self.logger.info(self.option("express_file").split(',')[0])
        self.cluster.set_options(options)
        self.logger.info(self.cluster.option('method'))
        self.cluster.on("end", self.set_db)
        self.cluster.run()

    def filter_file(self, infile, seq_list, output):
        with open(infile, 'rb') as r, open(output, 'wb') as w:
            w.write(r.readline())
            for line in r:
                if line.split('\t')[0] in seq_list:
                    w.write(line)

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_cluster = self.api.denovo_cluster  # #不确定,增加一个database 
        if self.option("cluster_method") == "hclust":
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
                    api_cluster.add_cluster_detail(cluster_id=self.option("cluster_id"), sub=sub, sub_path=sub_path)

                if re.match('samples_tree', sub_cluster):  # 找到sample_tree
                    sample_tree = os.path.join(hclust_path, sub_cluster)
                if re.match('genes_tree', sub_cluster):  # 找到gene_tree
                    gene_tree = os.path.join(hclust_path, sub_cluster)

            self.update_cluster(table_id=self.option("cluster_id"), genes=genes, sample_tree=sample_tree, gene_tree=gene_tree, specimen=specimen)
        else:
            kmeans_path = os.path.join(self.output_dir, "kmeans")
            sub_clusters = os.listdir(kmeans_path)
            genes = []
            for sub_cluster in sub_clusters:
                if re.match('subcluster', sub_cluster):
                    sub = sub_cluster.split("_")[1]
                    sub_path = os.path.join(kmeans_path, sub_cluster)
                    api_cluster.add_cluster_detail(cluster_id=self.option("cluster_id"), sub=sub, sub_path=sub_path)
                if re.search(r'kmeans_heatmap.xls', sub_cluster):
                    sub_cluster_path = os.path.join(kmeans_path, sub_cluster)
                    with open(sub_cluster_path, 'rb') as heatmap:
                        specimen = heatmap.readline()  # 第一行信息
                        heat_lst = heatmap.readlines()[1:]
                        for gene_num in heat_lst:
                            gene = gene_num.split('\t')[0]
                            genes.append(gene)  # 获取
            self.update_cluster(table_id=self.option("cluster_id"), genes=genes, specimen=specimen, gene_tree=None, sample_tree=None)
        self.end()

    def update_cluster(self, table_id, genes, specimen, gene_tree, sample_tree):
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_rna'
        collection = db['sg_denovo_cluster']
        if gene_tree:
            with open(gene_tree, 'rb') as g:
                gene_tree = g.readlines()[0].strip('\n')
        if sample_tree:
            with open(sample_tree, 'rb') as s:
                sample_tree = s.readlines()[0].strip('\n')
        collection.update({'_id': ObjectId(table_id)}, {'$set': {'specimen': specimen, 'genes': genes, 'gene_tree': gene_tree, 'sample_tree': sample_tree}})

    def run(self):
        self.run_cluster()
        super(ClusterWorkflow, self).run()
