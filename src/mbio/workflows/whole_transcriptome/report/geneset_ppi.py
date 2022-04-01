# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last update
from biocluster.workflow import Workflow
import os
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetPpiWorkflow(Workflow):
    """
    报告中进行ppi网络构建与分析时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_list", "type": "infile", "format": "ref_rna_v2.common"},  # 基因列表
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 1000},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "score", "type": "float", "default": 0.4},
            {"name": "seq", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "type", "type":"string"},
            {"name": "t2g_dir", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "ppi_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ppinetwork = self.add_module('ref_rna_v2.ppinetwork_analysis')
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/04 PPI_Analysis')
        self.inter_dirs = []
        '''
        if self.option('type') == 'G':
            self.set_options({'seq': self.trans_t2g()})
        '''

    def run_ppinetwork(self):
        # self.logger.info(self.option("diff_exp_gene").path)
        diff_exp_gene = os.path.join(self.work_dir, "geneset_list.txt")
        options = {
            'diff_exp_gene': self.option('geneset_list'),
            'species': self.option('species'),
            'seq': self.option('seq').prop['path'],
            'score': self.option('score'),
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.output_dir = self.ppinetwork.output_dir
        self.ppinetwork.run()

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

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/04 PPI_Analysis", "", "蛋白互作网络分析结果目录", 0]
        ]
        
        result_dir.add_relpath_rules([
            [".", "", "蛋白互作网络分析", 0],
            ["./ppinetwork_map", "", "基因id比对到string数据库文件目录", 0, "211082"], 
            ["./map", "", "与String数据库比对结果文件",0,"211560"],
            ["./map/seq_mapped.txt", "txt", "比对到数据库的序列信息表",0,"211561"],
            ["./map/unmapped_seq.txt", "txt", "未必对到数据库的序列信息",0,"211562"],
            ["./map/unmapped_db.txt", "txt", "数据库中该物种本次没有比对到的序列信息",0,"211563"],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件", 0, "211083"],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件", 0,0,"211564"],
            ["./ppinetwork_map/diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因列表", 0, "211085"],
            ["./ppinetwork_predict/all_nodes.txt", "txt", "PPI网络节点信息列表", 0, "211086"],
            ["./ppinetwork_predict/network_stats.txt", "txt", "PPI网络统计结果表", 0, "211087"],
            ["./ppinetwork_predict/interaction.txt", "txt", "PPI网络边信息列表", 0, "211088"],
            ["./ppinetwork_predict/interaction_detail.txt", "txt", "PPI网络边信息列表", 0, "211089"],
            ["./ppinetwork_predict/gene_gene.txt", "txt", "基因id与蛋白质对应表", 0, "211090"],
            ["./ppinetwork_topology/gene_interaction_network_centrality.txt", "txt", "PPI网络中心系数表", 0, "211091"],
            ["./ppinetwork_topology/gene_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表", 0, "211092"],
            ["./ppinetwork_topology/gene_interaction_network_transitivity.txt", "txt", "PPI网络传递性", 0, "211093"],
            ["./ppinetwork_topology/gene_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表", 0, "211094"],
            ["./ppinetwork_topology/gene_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络", 0, "211095"],
            ["./ppinetwork_topology/gene_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表", 0, "211096"],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetPpiWorkflow, self).end()

    def set_db(self):
        """
        报存分析结果到mongo数据库中
        """
        api_ppinetwork = self.api.api('whole_transcriptome.geneset_ppi')
        all_nodes_path = self.output_dir + '/ppinetwork_predict/all_nodes.txt'   # 画图节点属性文件
        interaction_path = self.output_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件
        network_stats_path = self.output_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
        network_centrality_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_centrality.txt'
        # network_clustering_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_clustering.txt'
        degree_distribution_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_degree_distribution.txt'
        network_node_degree_path = self.output_dir + '/ppinetwork_topology/gene_interaction_network_node_degree.txt'
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
        self.get_run_log()
        self.run_ppinetwork()
        super(GenesetPpiWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="geneset_ppi",
                                main_id=self.option('ppi_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def trans_t2g(self):
        t2seq = dict()
        g2t = dict()

        with open(self.option('t2g_dir').prop['path'], 'r') as t2g_file:
            for line in t2g_file.readlines():
                line = line.strip().split('\t')
                if line[2] == 'yes':
                    g2t[line[0]] = line[1]

        output = self.work_dir + '/gene.fa'
        exc_l = self.work_dir + '/expect.list'
        with open(self.option('seq').prop['path'], 'r') as seq_file, open(output, 'w') as gene_file, open(exc_l, 'w') as ect:
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
