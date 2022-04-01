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
        self.a_node = {
            "id": "",
            "name": "",
            "type": "",
        }
        self.a_edge = {
            "type": "",
            "score": 0,
            "db": "",
            "reg": "unknown",
            "corr": 0.3,
            "pvalue": 0,
            "padjust": 0
        }

    def get_tf_binding(self, tf_binding):
        '''
        获取转录因子 mirna结合结果
        '''
        with open(tf_binding, 'r') as f:
            f.readline()
            for line in f.readlines():
                cols = line.strip("\n").split("\t")
                tf = cols[2]
                mirnas = cols[0].split(";")
                tf_node = self.a_node.copy()
                tf_node["name"] = tf
                tf_node["id"] = tf
                tf_node["type"] = "tf"
                if tf not in self.tf_nodes:
                    self.tf_nodes.append(tf)
                if tf not in self.nodes:
                    self.nodes[tf] = tf_node

                for mirna in mirnas:
                    mirna_node = self.a_node.copy()
                    mirna_node["name"] = mirna
                    mirna_node["id"] = mirna
                    mirna_node["type"] = "mirna"
                    if mirna not in self.mirna_nodes:
                        self.mirna_nodes.append(mirna)
                    if mirna not in self.nodes:
                        self.nodes[mirna] = mirna_node

                    edge = self.a_edge.copy()
                    edge["type"] = "tf_mirna"
                    if cols[5] != "":
                        edge["score"] = cols[4]
                    if len(cols) >= 11 and cols[10] != "":
                        edge["db"] = "transmir"
                        if cols[10].startswith("Activation"):
                            edge["reg"] = "Activation"
                        elif cols[10].startswith("Repression"):
                            edge["reg"] = "Repression"
                        else:
                            pass

                    if tf + " " + mirna not in self.edges:
                        self.edges[tf + " " + mirna] = edge

    def get_mi_target(self, mi_target):
        '''
        获取转录因子 mirna结合结果
        '''
        with open(mi_target, 'r') as f:
            f.readline()
            for line in f.readlines():
                cols = line.strip("\n").split("\t")
                mirna = cols[0]
                gene = cols[2]
                mirna_node = self.a_node.copy()
                mirna_node["name"] = mirna
                mirna_node["id"] = mirna
                mirna_node["type"] = "mirna"
                if mirna not in self.mirna_nodes:
                    self.mirna_nodes.append(mirna)
                if mirna not in self.nodes:
                    self.nodes[mirna] = mirna_node

                gene_node = self.a_node.copy()
                if cols[3] == "":
                    gene_node["name"] = gene
                else:
                    gene_node["name"] = cols[3]
                gene_node["id"] = gene
                gene_node["type"] = "gene"
                if gene not in self.gene_nodes:
                    self.gene_nodes.append(gene)
                if gene not in self.nodes:
                    self.nodes[gene] = gene_node

                edge = self.a_edge.copy()
                edge["type"] = "mirna_gene"
                if len(cols) < 7:
                    edge["corr"] = "0"
                    edge["reg"] = "unknown"
                elif cols[6] != "":
                    if float(cols[6]) > 0:
                        edge["corr"] = cols[6]
                        edge["reg"] = "postive"
                        edge["pvalue"] = cols[7]
                        edge["pajust"] = cols[8]
                    else:
                        edge["corr"] = str(abs(float(cols[6])))
                        edge["reg"] = "negtive"
                        edge["pvalue"] = cols[7]
                        edge["padjust"] = cols[8]
                if cols[5] == "yes":
                    edge["db"] = "mitarbase"

                if mirna + " " + gene not in self.edges:
                    self.edges[mirna + " " + gene] = edge

    def parse_node(self):
        i = 0
        # print self.nodes
        for node in self.tf_nodes + self.mirna_nodes + self.gene_nodes:
            self.nodes[node].update({"node_num": i})
            i += 1

    def write_nodes_file(self, node_out):
        with open(node_out, 'w') as node_w:
            node_w.write("node_num\tnode_id\tnode_name\tnode_type\n")
            for node in self.tf_nodes + self.mirna_nodes + self.gene_nodes:
                print node
                node_w.write("\t".join([
                    str(self.nodes[node]["node_num"]),
                    self.nodes[node]["id"],
                    self.nodes[node]["name"],
                    self.nodes[node]["type"],
                ]) + "\n")

    def write_edges_file(self, edge_out):
        with open(edge_out, 'w') as edge_w:
            edge_w.write("from_id\tto_id\tfrom\tto\tfrom_name\tto_name\ttype\tscore\tdb\treg\tcorr\tpvalue\tpadjust\n")
            for edge, edge_dict in self.edges.items():
                nodes_from, nodes_to = edge.split(" ")
                edge_w.write("\t".join([
                    str(self.nodes[nodes_from]["node_num"]),
                    str(self.nodes[nodes_to]["node_num"]),
                    self.nodes[nodes_from]["id"],
                    self.nodes[nodes_to]["id"],
                    self.nodes[nodes_from]["name"],
                    self.nodes[nodes_to]["name"],
                    edge_dict["type"],
                    str(edge_dict["score"]),
                    edge_dict["db"],
                    edge_dict["reg"],
                    str(edge_dict["corr"]),
                    str(edge_dict["pvalue"]),
                    str(edge_dict["padjust"]),
                ]) + "\n")

if __name__ == "__main__":
    Net = NetMerge()
    if os.path.exists(sys.argv[1]):
        Net.get_tf_binding(sys.argv[1])
    Net.get_mi_target(sys.argv[2])
    Net.parse_node()
    Net.write_nodes_file(sys.argv[3])
    Net.write_edges_file(sys.argv[4])
