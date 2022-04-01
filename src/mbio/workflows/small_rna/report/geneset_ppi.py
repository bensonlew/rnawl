# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last update
from biocluster.workflow import Workflow
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import os,re
import json
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetPpiWorkflow(Workflow):
    """
    报告中进行ppi网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_list", "type": "infile", "format": "small_rna.common"},  # 基因列表
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 1000},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            # {"name": "seq", "type":"string"},
            {"name": "seq", "type": "infile", 'format': 'small_rna.fasta'},
            {"name": "type", "type":"string"},
            # {"name": "t2g_dir", "type":"string"},
            {"name": "t2g_dir", "type":"infile", 'format': 'small_rna.common'},
            {"name": "ppi_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/08 PPI_Analysis')
        self.inter_dirs = []
        self.ppinetwork = self.add_module('small_rna.ppinetwork_analysis')
        if self.option('type') == 'G':
            self.set_options({'seq': self.trans_t2g()})

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
        super(GenesetPpiWorkflow, self).send_log(data)

    def run_ppinetwork(self):
        # self.logger.info(self.option("diff_exp_gene").path)
        diff_exp_gene = os.path.join(self.work_dir, "geneset_list.txt")
        options = {
            'diff_exp_gene': self.option('geneset_list'),
            'species': self.option('species'),
            'seq': self.option('seq').path,
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.output_dir = self.ppinetwork.output_dir
        self.ppinetwork.run()


    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/08 PPI_Analysis", "", "靶基因蛋白互作分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白互作网络分析文件"],
            ["./map", "", "与String数据库比对结果文件"],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件"],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件"],
            # ["./ppinetwork_map/diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因列表"],
            ["./map/seq_mapped.txt", "txt", "比对到数据库的序列信息表",0],
            ["./map/unmapped_seq.txt", "txt", "未必对到数据库用的序列信息表",0],
            ["./map/unmapped_db.txt", "txt", "数据库中该物种本次没有比对到的序列信息表",0],
            ["./ppinetwork_predict/all_nodes.txt", "txt", "PPI网络节点信息列表"],
            ["./ppinetwork_predict/network_stats.txt", "txt", "PPI网络统计结果表"],
            ["./ppinetwork_predict/interaction.txt", "txt", "PPI网络边信息列表"],
            ["./ppinetwork_predict/interaction_detail.txt", "txt", "PPI网络边信息列表"],
            ["./ppinetwork_predict/gene_gene.txt", "txt", "基因id与蛋白质对应表"],
            ["./ppinetwork_topology/gene_interaction_network_centrality.txt", "txt", "PPI网络中心系数表"],
            ["./ppinetwork_topology/gene_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表"],
            ["./ppinetwork_topology/gene_interaction_network_transitivity.txt", "txt", "PPI网络传递性"],
            ["./ppinetwork_topology/gene_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表"],
            ["./ppinetwork_topology/gene_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络"],
            ["./ppinetwork_topology/gene_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetPpiWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_ppinetwork = self.api.api('small_rna.geneset_ppi')
        all_nodes_path = self.output_dir + '/ppinetwork_predict/all_nodes.txt'   # 画图节点属性文件
        interaction_path = self.output_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件

        #生成绘图所需json文件
        nodes_df = pd.read_table(all_nodes_path)
        nodes_df.sort_values(['degree'], ascending=[False], inplace=True)
        nodes_df['node_id'] = range(0, len(nodes_df))
        self.node2acc = dict(zip(nodes_df['STRING_id'], nodes_df['accession_id']))
        self.node2id = dict(zip(nodes_df['STRING_id'], nodes_df['node_id']))
        self.acc2id = dict(zip(nodes_df['accession_id'], nodes_df['node_id']))
        nodes_df = nodes_df[["method", "accession_id","node"]]
        # nodes_df["name"] = nodes_df.apply(
        #     lambda x: id2gene_name[x["accession_id"]] if id2gene_name[x["accession_id"]] != "-" else x["accession_id"],
        #     axis=1)
        # nodes_df["name"] = nodes_df.apply(lambda x:  x["node"] if x["node"] != "-" else x["accession_id"], axis=1)
        nodes_df = nodes_df.rename({"method": "group", "accession_id": "id","node":"name"}, axis=1)
        nodes_df["name"] = nodes_df.apply(lambda x: x["name"] if x["name"] != "-" else x["id"], axis=1)
        edge_df = pd.read_table(interaction_path)
        edge_df['to_accession'] = edge_df['to'].map(lambda x: self.node2acc[x])
        edge_df['from_id'] = edge_df['from'].map(lambda x: self.node2id[x])
        edge_df['to_id'] = edge_df['to'].map(lambda x: self.node2id[x])
        edge_df = edge_df.rename({"from_id": "source", "to_id": "target", "combined_score": "distance"}, axis=1)
        self.logger.info("{}".format(str(edge_df.columns)))
        edge_df = edge_df[["source", "target", 'distance']]
        with open(os.path.join(self.ppinetwork.output_dir, "ppinetwork_predict/network.json"), "w") as f:
            json.dump(dict(nodes=nodes_df.to_dict("records"), links=edge_df.to_dict("records")), f)



        network_stats_path = self.output_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
        network_centrality_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_centrality.txt'
        # network_clustering_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_clustering.txt'
        degree_distribution_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_degree_distribution.txt'
        network_node_degree_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_node_degree.txt'
        if not os.path.isfile(all_nodes_path):
            raise Exception("找不到报告文件:{}".format(all_nodes_path))
        if not os.path.isfile(interaction_path):
            raise Exception("找不到报告文件:{}".format(interaction_path))
        if not os.path.isfile(network_stats_path):
            raise Exception("找不到报告文件:{}".format(network_stats_path))
        # if not os.path.isfile(network_clustering_path):
        #     raise Exception("找不到报告文件:{}".format(network_clustering_path))
        if not os.path.isfile(network_centrality_path):
            raise Exception("找不到报告文件:{}".format(network_centrality_path))
        if not os.path.isfile(degree_distribution_path):
            raise Exception("找不到报告文件:{}".format(degree_distribution_path))
        if not os.path.isfile(network_node_degree_path):
            raise Exception("找不到报告文件:{}".format(network_node_degree_path))

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
        network = self._sheet.output + '/ppinetwork_predict/network.json'
        api_ppinetwork.update_db_record('sg_geneset_ppi', self.option('ppi_id'), json_dir=network, has_json="yes",
                                        status="end")

        print('end insert')
        self.end()

    def run(self):
        self.get_run_log()
        self.run_ppinetwork()
        super(GenesetPpiWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_geneset_ppi", main_id=self.option('ppi_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def trans_t2g(self):
        t2seq = dict()
        g2t = dict()

        # with open(self.option('t2g_dir'), 'r') as t2g_file:
        with open(self.option('t2g_dir').path, 'r') as t2g_file:
            for line in t2g_file.readlines():
                line = line.strip().split('\t')
                if line[2] == 'yes':
                    g2t[line[0]] = line[1]

        output = os.path.join(self.work_dir, 'gene.fa')
        exc_l = os.path.join(self.work_dir, 'expect.list')
        # with open(self.option('seq'), 'r') as seq_file, open(output, 'w') as gene_file, open(exc_l, 'w') as ect:
        with open(self.option('seq').path, 'r') as seq_file, open(output, 'w') as gene_file, open(exc_l, 'w') as ect:
            seq_list = seq_file.read().split('\n>')
            # seq_list=seq_list.split('\n>')
            for block in seq_list:
                block = block.split('\n')
                g = block[0].lstrip('>').split(' ')[0]
                seq = '\n'.join(block[1:]) + '\n'
                if g != '':
                    t2seq[g] = seq
            for t in t2seq:
                try:
                    gene_file.write('>{}\n{}'.format(g2t[t], t2seq[t]))
                except:
                    ect.write(t + '\n')
        return output
