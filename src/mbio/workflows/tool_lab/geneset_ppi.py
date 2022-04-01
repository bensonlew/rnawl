# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last update
from biocluster.workflow import Workflow
import json
import pandas as pd
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import os,re
import shutil
from collections import OrderedDict
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
import glob
from biocluster.config import Config
from biocluster.file import getsize, exists
from biocluster.file import download
import time
from bson.objectid import ObjectId


class GenesetPpiWorkflow(Workflow):
    """
    报告中进行ppi网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_list", "type": "infile", "format": "ref_rna_v2.common"},  # 基因列表
            {'name': 'species_name', 'type': 'string'},
            {"name": "species", "type": "string", "default": 'Homo sapiens(9606)'},
            {"name": "combine_score", "type": "int", "default": 1000},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "score", "type": "float", "default": 0.4},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.ppinetwork = self.add_module('tool_lab.ppinetwork_analysis')

    def run_ppinetwork(self):
        self.species = re.findall(r'\((.*?)\)', self.option('species'))[0]
        options = {
            'diff_exp_gene': self.geneset_file,
            'species': int(self.species),
            'score': self.option('score'),
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.ppinetwork.run()


    def end(self):
        os.rename(os.path.join(self.ppinetwork.output_dir,"ppinetwork_topology","gene_interaction_network_centrality.txt"),os.path.join(self.ppinetwork.output_dir,"ppinetwork_topology","PPI_centrality.txt"))
        os.rename(os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "gene_interaction_network_degree_distribution.txt"), os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "PPI_degree_distribution.txt"))
        os.rename(os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "gene_interaction_network_node_degree.txt"), os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "PPI_node_degree.txt"))
        result_dir = self.add_upload_dir(self.ppinetwork.output_dir)

        self.inter_dirs = [
            ["06 Advanced_Analysis", "", "高级分析结果目录",0],
            ["06 Advanced_Analysis/05 PPI_Analysis", "", "蛋白互作网络分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白互作网络分析", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ["./ppinetwork_map", "", "基因id比对到string数据库文件目录", 0, "211082"], 
            ["./map", "", "与String数据库比对结果目录",0,"211560"],
            ["./map/seq_mapped.txt", "txt", "比对到数据库的序列信息表",0],
            ["./map/unmapped_seq.txt", "txt", "未必对到数据库用的序列信息表",0],
            ["./map/unmapped_db.txt", "txt", "数据库中该物种本次没有比对到的序列信息表",0],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件", 0],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件", 0,0,"211564"],
            ["./ppinetwork_map/diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因列表", 0, "211085"],
            ["./ppinetwork_predict/all_nodes.txt", "txt", "PPI网络节点信息列表", 0, "211086"],
            ["./ppinetwork_predict/network_stats.txt", "txt", "PPI网络统计结果表", 0, "211087"],
            # ["./ppinetwork_predict/interaction.txt", "txt", "PPI网络边信息列表", 0, "211088"],
            ["./ppinetwork_predict/interaction_detail.txt", "txt", "PPI网络边信息列表", 0, "211089"],
            # ["./ppinetwork_predict/gene_gene.txt", "txt", "基因id与蛋白质对应表", 0, "211090"],
            ["./ppinetwork_topology/PPI_centrality.txt", "txt", "PPI网络中心系数表", 0, "211091"],
            ["./ppinetwork_topology/*.pdf", "txt", "PPI网络统计图", 0],
            # ["./ppinetwork_topology/gene_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表", 0, "211092"],
            # ["./ppinetwork_topology/gene_interaction_network_transitivity.txt", "txt", "PPI网络传递性", 0, "211093"],
            ["./ppinetwork_topology/PPI_degree_distribution.txt", "txt", "PPI网络度分布表", 0, "211094"],
            # ["./ppinetwork_topology/gene_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络", 0, "211095"],
            ["./ppinetwork_topology/PPI_node_degree.txt", "txt", "PPI网络节点度属性表", 0, "211096"]
        ])
        super(GenesetPpiWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中

        """
        api_ppinetwork = self.api.api('tool_lab.geneset_ppi')
        all_nodes_path = self.ppinetwork.output_dir + '/ppinetwork_predict/all_nodes.txt'   # 画图节点属性文件
        interaction_path = self.ppinetwork.output_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件

        # nodes_df = pd.read_table(all_nodes_path)
        # nodes_df.sort_values(['degree'], ascending=[False], inplace=True)
        # nodes_df['node_id'] = range(0, len(nodes_df))
        # self.node2acc = dict(zip(nodes_df['STRING_id'], nodes_df['accession_id']))
        # self.node2id = dict(zip(nodes_df['STRING_id'], nodes_df['node_id']))
        # self.acc2id = dict(zip(nodes_df['accession_id'], nodes_df['node_id']))
        # nodes_df = nodes_df[["method", "accession_id"]]
        # nodes_df = nodes_df.rename({"method": "group", "accession_id": "id"}, axis=1)
        # edge_df = pd.read_table(interaction_path)
        # edge_df['to_accession'] = edge_df['to'].map(lambda x: self.node2acc[x])
        # edge_df['from_id'] = edge_df['from'].map(lambda x: self.node2id[x])
        # edge_df['to_id'] = edge_df['to'].map(lambda x: self.node2id[x])
        # edge_df = edge_df.rename({"from_id": "source", "to_id": "target","combined_score":"distance"}, axis=1)
        # self.logger.info("{}".format(str(edge_df.columns)))
        # edge_df = edge_df[["source", "target",'distance']]
        # with open(os.path.join(self.ppinetwork.output_dir,"ppinetwork_predict/network.json"), "w") as f:
        #     json.dump(dict(nodes=nodes_df.to_dict("records"),links=edge_df.to_dict("records")), f)
        network_stats_path = self.ppinetwork.output_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
        network_centrality_path = self.ppinetwork.output_dir + '/ppinetwork_topology/gene_interaction_network_centrality.txt'
        # network_clustering_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_clustering.txt'
        degree_distribution_path = self.ppinetwork.output_dir + '/ppinetwork_topology/gene_interaction_network_degree_distribution.txt'
        network_node_degree_path = self.ppinetwork.output_dir + '/ppinetwork_topology/gene_interaction_network_node_degree.txt'
        if not os.path.isfile(all_nodes_path):
            self.set_error("找不到报告文件:%s", variables = (all_nodes_path), code = "13701801")
        if not os.path.isfile(interaction_path):
            self.set_error("找不到报告文件:%s", variables = (interaction_path), code = "13701802")
        if not os.path.isfile(network_stats_path):
            self.set_error("找不到报告文件:%s", variables = (network_stats_path), code = "13701803")
        # if not os.path.isfile(network_clustering_path):
        #     self.set_error("找不到报告文件:{}".format(network_clustering_path))
        if not os.path.isfile(network_centrality_path):
            self.set_error("找不到报告文件:%s", variables = (network_centrality_path), code = "13701804")
        if not os.path.isfile(degree_distribution_path):
            self.set_error("找不到报告文件:%s", variables = (degree_distribution_path), code = "13701805")
        if not os.path.isfile(network_node_degree_path):
            self.set_error("找不到报告文件:%s", variables = (network_node_degree_path), code = "13701806")

        print('start insert')
        api_ppinetwork.add_node_table(file_path=all_nodes_path, table_id=self.option("main_id"))   # 节点的属性文件（画网络图用）
        api_ppinetwork.add_edge_table(file_path=interaction_path, table_id=self.option("main_id"))  # 边信息
        api_ppinetwork.update_db_record('sg_geneset_ppi', self.option('main_id'),
                                   status="end")
        print('end insert')
        self.end()

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            self.geneset_file = self.download_s3_file(self.file_path, 'geneset_table.txt')
        if self.option('source') == 'tool_lab':
            self.geneset_file = self.option('geneset_list').prop['path']
        self.run_ppinetwork()
        super(GenesetPpiWorkflow, self).run()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        return file_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path
