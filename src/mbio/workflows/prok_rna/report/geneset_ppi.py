# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last update
from biocluster.workflow import Workflow
import os

class ProteinsetPpiWorkflow(Workflow):
    """
    报告中进行ppi网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_list", "type": "infile", "format": "itraq_and_tmt.common"},  # 基因列表
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 1000},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "update_info", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "seq", type="infile", "format": "prok_rna.common"},
            {"name": "ppi_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ppinetwork = self.add_module('itraq_and_tmt.ppinetwork_analysis')

    def run_ppinetwork(self):
        # self.logger.info(self.option("diff_exp_gene").path)
        diff_exp_gene = os.path.join(self.work_dir, "proteinset_list.txt")
        options = {
            'diff_exp_gene': self.option('proteinset_list'),
            'species': self.option('species'),
            'seq': self.option('seq').prop['path'],
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.output_dir = self.ppinetwork.output_dir
        self.ppinetwork.run()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PPI网络分析结果输出目录"],
            ["./ppinetwork_map", "", "基因id比对到string数据库文件目录"],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件目录"],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件目录"],
            ["./ppinetwork_map/diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因列表"],
            ["./ppinetwork_predict/all_nodes.txt", "txt", "PPI网络节点信息列表"],
            ["./ppinetwork_predict/network_stats.txt", "txt", "PPI网络统计结果表"],
            ["./ppinetwork_predict/interaction.txt", "txt", "PPI网络边信息列表"],
            ["./ppinetwork_predict/interaction_detail.txt", "txt", "PPI网络边信息列表"],
            ["./ppinetwork_predict/gene_protein.txt", "txt", "基因id与蛋白质对应表"],
            ["./ppinetwork_topology/protein_interaction_network_centrality.txt", "txt", "PPI网络中心系数表"],
            ["./ppinetwork_topology/protein_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表"],
            ["./ppinetwork_topology/protein_interaction_network_transitivity.txt", "txt", "PPI网络传递性"],
            ["./ppinetwork_topology/protein_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表"],
            ["./ppinetwork_topology/protein_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络"],
            ["./ppinetwork_topology/protein_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表"]
        ])
        super(ProteinsetPpiWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_ppinetwork = self.api.api('itraq_and_tmt.proteinset_ppi')
        all_nodes_path = self.output_dir + '/ppinetwork_predict/all_nodes.txt'   # 画图节点属性文件
        interaction_path = self.output_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件
        network_stats_path = self.output_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
        network_centrality_path = self.output_dir + '/ppinetwork_topology/protein_interaction_network_centrality.txt'
        # network_clustering_path = self.output_dir + '/ppinetwork_topology/protein_interaction_network_clustering.txt'
        degree_distribution_path = self.output_dir + '/ppinetwork_topology/protein_interaction_network_degree_distribution.txt'
        network_node_degree_path = self.output_dir + '/ppinetwork_topology/protein_interaction_network_node_degree.txt'
        if not os.path.isfile(all_nodes_path):
            raise OptionError("找不到报告文件：%s", variables = (all_nodes_path), code = "15000501")
        if not os.path.isfile(interaction_path):
            raise OptionError("找不到报告文件：%s", variables = (interaction_path), code = "15000502")
        if not os.path.isfile(network_stats_path):
            raise OptionError("找不到报告文件：%s", variables = (network_stats_path), code = "15000503")
        # if not os.path.isfile(network_clustering_path):
        #     raise Exception("找不到报告文件:{}".format(network_clustering_path))
        if not os.path.isfile(network_centrality_path):
            raise OptionError("找不到报告文件：%s", variables = (network_centrality_path), code = "15000504")
        if not os.path.isfile(degree_distribution_path):
            raise OptionError("找不到报告文件：%s", variables = (degree_distribution_path), code = "15000505")
        if not os.path.isfile(network_node_degree_path):
            raise OptionError("找不到报告文件：%s", variables = (network_node_degree_path), code = "15000506")

        print 'start insert'
        api_ppinetwork.add_node_table(file_path=all_nodes_path, table_id=self.option("ppi_id"))   # 节点的属性文件（画网络图用）
        api_ppinetwork.add_edge_table(file_path=interaction_path, table_id=self.option("ppi_id"))  # 边信息
        api_ppinetwork.add_network_attributes(file2_path=network_stats_path,
                                              table_id=self.option("ppi_id"))  # 网络全局属性
        # api_ppinetwork.add_network_cluster_degree(file1_path=network_node_degree_path,
                                                  # file3_path=all_nodes_path,
                                                  # table_id=self.option("ppi_id"))
        # 节点的聚类与degree，画折线图
        api_ppinetwork.add_network_centrality(file_path=network_centrality_path, file1_path=all_nodes_path,
                                              table_id=self.option("ppi_id"))  # 中心信息
        api_ppinetwork.add_degree_distribution(file_path=degree_distribution_path, table_id=self.option("ppi_id"))
        print 'end insert'
        self.end()

    def run(self):
        self.run_ppinetwork()
        super(ProteinsetPpiWorkflow, self).run()
