# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET

class NetMerge(object):
    def __init__(self):
        self.nodes = dict()
        self.edges = dict()
        self.tf_nodes = list()
        self.mirna_nodes = list()
        self.target_nodes = list()
        self.gene_nodes = list()
        self.ce_nodes = list()
        self.ce_edges = list()
        self.mi_ce_dict = dict()  # 判断一对mi_ce影响另外的什么
        self.a_node = {
            "id": "",
            "name": "",
            "type": "",
        }
        self.a_edge = {
            "type": "",
            "score": 0,
            "energy": 0,
            "direction": 0,
            "mifam_num": 0,
            "p_value": 1,
            "corr": 0,
            "corr_pvalue": 1,
            "corr_padjust": 1,
            "ref": "unknown",
            "is_ce": "no",
            "from_type": "mRNA",
            "to_type": "mRNA",
            "ind_list": [],

        }

    def get_ce_binding(self, ce_file):
        with open(ce_file, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip("\n").split("\t")
                ce_node1 = self.a_node.copy()
                ce_node2 = self.a_node.copy()
                ce_node1["id"] = cols[0]
                ce_node2["id"] = cols[1]
                if cols[2] != "-":
                    ce_node1["name"] = cols[2]
                else:
                    ce_node1["name"] = cols[0]
                if cols[3] != "-":
                    ce_node2["name"] = cols[3]
                else:
                    ce_node2["name"] = cols[1]
                ce_node2["name"] = cols[3]
                ce_node1["type"] = cols[4]
                ce_node2["type"] = cols[5]
                if cols[0] in self.nodes:
                    pass
                else:
                    self.nodes[cols[0]] = ce_node1
                if cols[1] in self.nodes:
                    pass
                else:
                    self.nodes[cols[1]] = ce_node2

                edge = self.a_edge.copy()
                if float(cols[7]) > 0:
                    reg = "positive"
                elif float(cols[7]) < 0:
                    ref = "negative"
                else:
                    ref = "unknown"
                edge.update({
                    "type": "cerna",
                    "direction": 0,
                    "mifam_num": int(cols[6]),
                    "p_value": float(cols[7]),
                    "corr": abs(float(cols[8])),
                    "corr_pvalue": float(cols[9]),
                    "corr_padjust": float(cols[10]),
                    "reg": reg,
                    "is_ce": "yes",
                    "from_type": ce_node1["type"],
                    "to_type": ce_node2["type"]

                })
                self.edges[cols[0] + " " + cols[1]] = edge

                ce_mitarget = cols[11]
                mi_list = list()
                for fam in ce_mitarget.split(";"):
                    mi_ce_corrs = list()
                    for mirna in fam.split(":")[-1].split(","):
                        mirna_name = mirna.split("(")[0]
                        self.ce_edges.append(mirna_name + " " + cols[0])
                        self.ce_edges.append(mirna_name + " " + cols[1])
                        mi_list.append(mirna_name)
                        mi_ce = mirna_name + "||" + ce_node1["id"]
                        if mi_ce in self.mi_ce_dict:
                            self.mi_ce_dict[mi_ce].add(ce_node2["id"])
                        else:
                            self.mi_ce_dict[mi_ce] = set([ce_node2["id"]])

                        mi_ce = mirna_name + "||" + ce_node2["id"]
                        if mi_ce in self.mi_ce_dict:
                            self.mi_ce_dict[mi_ce].add(ce_node1["id"])
                        else:
                            self.mi_ce_dict[mi_ce] = set([ce_node1["id"]])

                # 存入共有的miRNA 列表方便搜索
                edge.update({"ind_list": mi_list})

    def get_mi_target(self, mi_mrna):
        '''
        获取 mirna 靶基因结果， 只保留ce中预测的靶基因
        2019.10.29 改为保留所有的靶基因
        '''
        target_ce = {
            "mi_mrna": mi_mrna,
        }
        for ce_type, ce_file in target_ce.items():
            with open(ce_file, 'r') as f:
                header = f.readline()
                header_list = header.strip("\n").split("\t")

                try:
                    score_index = header_list.index("score")
                except:
                    score_index = None
                try:
                    energy_index = header_list.index("energy")
                except:
                    energy_index = None
                try:
                    corr_index = header_list.index("corr")
                except:
                    energy_index = None
                try:
                    pvalue_index = header_list.index("pvalue")
                except:
                    pvalue_index = None
                try:
                    qvalue_index = header_list.index("qvalue")
                except:
                    qvalue_index = None
                for line in f.readlines():
                    cols = line.strip("\n").split("\t")
                    mirna = cols[0]
                    target = cols[1]
                    gene = cols[3]

                    # 确定点的属性
                    mirna_node = self.a_node.copy()
                    mirna_node["name"] = mirna
                    mirna_node["id"] = mirna
                    mirna_node["type"] = "mirna"

                    ce_node = self.a_node.copy()
                    ce_node["name"] = cols[4]
                    ce_node["id"] = target
                    ce_node["type"] = cols[2]
                    # if gene in self.nodes or target in self.nodes:
                    if mirna not in self.nodes:
                        self.nodes[mirna] = mirna_node
                    else:
                        pass

                    if target not in self.nodes:
                        self.nodes[target] = ce_node
                    else:
                        pass
                    edge = self.a_edge.copy()
                    edge["type"] = "mi_" + cols[2].lower()
                    edge["from_type"] = "miRNA"
                    edge["to_type"]= ce_node["type"]
                    mi_ce = mirna + "||" + target
                    if mi_ce in self.mi_ce_dict:
                        edge.update({"ind_list": list(self.mi_ce_dict[mi_ce])})

                    try:
                        edge["corr"] = cols[corr_index]
                        edge["corr_pvalue"] = cols[pvalue_index]
                        edge["corr_padjust"] = cols[qvalue_index]
                    except:
                        edge["corr"] = 0
                        edge["corr_pvalue"] = 1
                        edge["corr_padjust"] = 1

                    if len(cols) < 11:
                        # edge["corr"] = "0"
                        edge["reg"] = "unknown"
                        edge["mifam_num"] = 0
                        edge["direction"] = 1
                        # edge["score"] = float(cols[6])
                        try:
                            edge["energy"] = float(cols[7])
                        except:
                            pass

                    elif cols[8] != "":
                        # edge["corr"] = "0"
                        edge["reg"] = "unknown"
                        edge["mifam_num"] = 0
                        edge["direction"] = 1
                        # edge["corr_pvalue"] = 1
                        # edge["corr_padjust"] = 1
                        try:
                            edge["score"] = float(cols[score_index])
                        except:
                            edge["score"] = float(999)

                        try:
                            edge["energy"] = float(cols[energy_index])
                        except:
                            edge["energy"] = float(-99)

                        if float(cols[8]) > 0:
                            # edge["corr"] = abs(float(cols[8]))
                            edge["reg"] = "postive"
                            # edge["corr_pvalue"] = float(cols[9])
                            # edge["corr_padjust"] = float(cols[10])
                        else:
                            # edge["corr"] = abs(float(cols[8]))
                            edge["reg"] = "negative"
                            # edge["corr_pvalue"] = float(cols[9])
                            # edge["corr_padjust"] = float(cols[10])

                    if mirna + " " + gene in self.ce_edges:
                        edge["is_ce"] = "yes"
                    elif mirna + " " + target in self.ce_edges:
                        edge["is_ce"] = "yes"

                    if gene in self.nodes:
                        if mirna + " " + gene not in self.edges:
                            self.edges[mirna + " " + gene] = edge
                        else:
                            pass
                    elif target in self.nodes:
                        if mirna + " " + gene not in self.edges:
                            self.edges[mirna + " " + target] = edge
                        else:
                            pass

    def parse_node(self):
        i = 0
        # print self.nodes
        for node in sorted(self.nodes.keys(), key=lambda x:self.nodes[x]['type']):
            self.nodes[node].update({"node_num": i})
            i += 1

    def write_nodes_file(self, node_out):
        with open(node_out, 'w') as node_w:
            node_w.write("node_num\tnode_id\tnode_name\tnode_type\n")
            for node in sorted(self.nodes.keys(), key=lambda x:self.nodes[x]['type']):
                node_w.write("\t".join([
                    str(self.nodes[node]["node_num"]),
                    self.nodes[node]["id"],
                    self.nodes[node]["name"],
                    self.nodes[node]["type"],
                ]) + "\n")

    def write_edges_file(self, edge_out):
        with open(edge_out, 'w') as edge_w:
            edge_w.write("from_id\tto_id\tfrom\tto\tfrom_name\tto_name\tfrom_type\tto_type\ttype\tminum\tind_list\tpvalue\treg\tcorr\tcorr_pvalue\tcorr_padjust\tis_ce\tscore\tenergy\n")
            for edge, edge_dict in self.edges.items():
                nodes_from, nodes_to = edge.split(" ")
                edge_w.write("\t".join([
                    str(self.nodes[nodes_from]["node_num"]),
                    str(self.nodes[nodes_to]["node_num"]),
                    self.nodes[nodes_from]["id"],
                    self.nodes[nodes_to]["id"],
                    self.nodes[nodes_from]["name"],
                    self.nodes[nodes_to]["name"],
                    self.nodes[nodes_from]["type"],
                    self.nodes[nodes_to]["type"],
                    edge_dict["type"],
                    str(edge_dict["mifam_num"]),
                    str(",".join(edge_dict["ind_list"])),
                    str(edge_dict["p_value"]),
                    edge_dict["reg"],
                    str(edge_dict["corr"]),
                    str(edge_dict["corr_pvalue"]),
                    str(edge_dict["corr_padjust"]),
                    str(edge_dict["is_ce"]),
                    str(edge_dict["score"]),
                    str(edge_dict["energy"]),
                ]) + "\n")

if __name__ == "__main__":
    Net = NetMerge()
    if os.path.exists(sys.argv[1]):
        Net.get_ce_binding(sys.argv[1])
    Net.get_mi_target(sys.argv[2])
    Net.parse_node()
    Net.write_nodes_file(sys.argv[3])
    Net.write_edges_file(sys.argv[4])
