# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last update
from biocluster.workflow import Workflow
import os
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json
import glob

class ProteinsetPpiWorkflow(Workflow):
    """
    报告中进行ppi网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_list", "type": "infile", "format": "itraq_and_tmt.common"},  # 蛋白列表
            {"name": "species", "type": "int", "default": 9606},
            {"name": "proteinset_id", "type": "string"},
            {"name": "combine_score", "type": "int", "default": 1000},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "score", "type": "float", "default": 0.4},
            {"name": "update_info", "type": "string"},
            {"name": "seq", "type": "infile", "format": "itraq_and_tmt.common"},
            {"name": "ppi_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ppinetwork = self.add_module('itraq_and_tmt.ppinetwork_analysis')
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/05_PPI')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(ProteinsetPpiWorkflow, self).send_log(data)

    def run_ppinetwork(self):
        # self.logger.info(self.option("diff_exp_gene").path)
        diff_exp_gene = os.path.join(self.work_dir, "proteinset_list.txt")
        options = {
            'diff_exp_gene': self.option('proteinset_list').prop['path'],
            'species': self.option('species'),
            'seq': self.option('seq').prop['path'],
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.output_dir = self.ppinetwork.output_dir
        self.ppinetwork.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        
        network_centrality = os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "protein_interaction_network_centrality.txt")
        network_degree = os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", "protein_interaction_network_degree_distribution.txt")
        chart.chart_ppi_centrality_degree_web(network_centrality, network_degree)

        # "ppi_centrality_all":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_all/ppinetwork_topology/protein_interaction_network_centrality.txt",
        # "ppi_centrality_up":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_up/ppinetwork_topology/protein_interaction_network_centrality.txt",
        # "ppi_centrality_down":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_down/ppinetwork_topology/protein_interaction_network_centrality.txt",
        # "ppi_degree_all":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_all/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",
        # "ppi_degree_up":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_up/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",
        # "ppi_degree_down":"/mnt/ilustre/users/sanger-dev/workspace/20210430/Itraqtmt_s5jn_s0ieqm5rf9dfelc1t88ajh/DiffPpi/output/{group_name}_down/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",

        # if "ppi_centrality_all" in chart_json and "ppi_centrality_up" in chart_json and "ppi_centrality_down" in chart_json:
        #     ppi_centrality_all_files = [chart_json["ppi_centrality_all"].format(group_name=group) for group in chart_json["diff_group"]]
        #     ppi_centrality_up_files = [chart_json["ppi_centrality_up"].format(group_name=group) for group in chart_json["diff_group"]]
        #     ppi_centrality_down_files = [chart_json["ppi_centrality_down"].format(group_name=group) for group in chart_json["diff_group"]]
        #     if files_all_exit(ppi_centrality_all_files) and files_all_exit(ppi_centrality_up_files) and files_all_exit(ppi_centrality_down_files):
        #         self.chart_ppi_centrality(chart_json["diff_group"], ppi_centrality_all_files, ppi_centrality_up_files, ppi_centrality_down_files)

        # if "ppi_degree_all" in chart_json and "ppi_degree_up" in chart_json and "ppi_degree_down" in chart_json:
        #     ppi_degree_all_files = [chart_json["ppi_degree_all"].format(group_name=group) for group in chart_json["diff_group"]]
        #     ppi_degree_up_files = [chart_json["ppi_degree_up"].format(group_name=group) for group in chart_json["diff_group"]]
        #     ppi_degree_down_files = [chart_json["ppi_degree_down"].format(group_name=group) for group in chart_json["diff_group"]]
        #     if files_all_exit(ppi_degree_all_files) and files_all_exit(ppi_degree_up_files) and files_all_exit(ppi_degree_down_files):
        #         self.chart_ppi_degree(chart_json["diff_group"], ppi_degree_all_files, ppi_degree_up_files, ppi_degree_down_files)

        # chart.chart_sample_corr(sample_corr_matrix, sample_corr_tree)
        chart.to_pdf()

        # move pdf to result dir
        for i in ["ppi.degree.line.showCurve.pdf", "ppi.centrality.line.showCurve.pdf"]:
            if os.path.exists(os.path.join(self.work_dir, i)):
                os.link(os.path.join(self.work_dir, i), os.path.join(self.ppinetwork.output_dir, "ppinetwork_topology", i[:-13]+'pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/05_PPI", "", "蛋白互作网络分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "PPI网络分析结果输出目录", 0, "220083"],
            ["./ppinetwork_map", "", "accession_id比对到string数据库文件目录", 0, "220084"],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件目录", 0, "220085"],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件目录", 0, "220086"],
            ["./ppinetwork_map/diff_exp_mapped.txt", "txt", "含有STRINGid的差异蛋白列表", 0, "220087"],
            ["./ppinetwork_predict/all_nodes.txt", "txt", "PPI网络节点信息列表", 0, "220088"],
            ["./ppinetwork_predict/network_stats.txt", "txt", "PPI网络统计结果表", 0, "220089"],
            ["./ppinetwork_predict/interaction.txt", "txt", "PPI网络边信息列表", 0, "220090"],
            ["./ppinetwork_predict/interaction_detail.txt", "txt", "PPI网络边信息列表", 0, "220091"],
            ["./ppinetwork_predict/gene_protein.txt", "txt", "accession_idid与蛋白质对应表", 0, "220092"],
            ["./ppinetwork_topology/protein_interaction_network_centrality.txt", "txt", "PPI网络中心系数表", 0, "220093"],
            ["./ppinetwork_topology/protein_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表", 0, "220094"],
            ["./ppinetwork_topology/protein_interaction_network_transitivity.txt", "txt", "PPI网络传递性", 0, "220095"],
            ["./ppinetwork_topology/protein_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表", 0, "220096"],
            ["./ppinetwork_topology/protein_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络", 0, "220097"],
            ["./ppinetwork_topology/protein_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表", 0, "220098"],
            ["./ppinetwork_topology/*degree*pdf", "", "PPI网络节点聚类系数图", 0],
            ["./ppinetwork_topology/*centrality*pdf", "","PPI网络中心系数图", 0],
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
            self.set_error("找不到报告文件:%s", variables = (all_nodes_path), code = "12501501")
        if not os.path.isfile(interaction_path):
            self.set_error("找不到报告文件:%s", variables = (interaction_path), code = "12501502")
        if not os.path.isfile(network_stats_path):
            self.set_error("找不到报告文件:%s", variables = (network_stats_path), code = "12501503")
        # if not os.path.isfile(network_clustering_path):
        #     self.set_error("找不到报告文件:{}".format(network_clustering_path))
        if not os.path.isfile(network_centrality_path):
            self.set_error("找不到报告文件:%s", variables = (network_centrality_path), code = "12501504")
        if not os.path.isfile(degree_distribution_path):
            self.set_error("找不到报告文件:%s", variables = (degree_distribution_path), code = "12501505")
        if not os.path.isfile(network_node_degree_path):
            self.set_error("找不到报告文件:%s", variables = (network_node_degree_path), code = "12501506")

        print('start insert')
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
        print('end insert')
        self.end()

    def run(self):
        self.run_ppinetwork()
        super(ProteinsetPpiWorkflow, self).run()
