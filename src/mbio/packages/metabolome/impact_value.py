# -*- coding: UTF-8 -*-
import os
import re
import glob
import pandas as pd
import xml.etree.ElementTree as ET
from collections import defaultdict
import networkx as nx
import json
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="prepare xml get Topology data")
parser.add_argument("-org", type=str, required=True, help="2019-04.hsa")
args = parser.parse_args()

db_dir = '/mnt/ilustre/users/ting.kuang/ITRAQ/db/KEGG/'
org = args.org

def cal_impactvalue(org):
    pathways = []
    xml_2_network = {}
    xml_2_imp = {}
    xml_2_nodes = {}
    xml_2_nodes_num = {}
    for xml_file in glob.glob(db_dir + org + '/*.xml'):
        _, xml_list = os.path.split(xml_file)
        xml_list = xml_list.split('.')[0]
        pathways.append(xml_list)
        tree = ET.parse(xml_file)
        root = tree.getroot()
        networks_list = []
        DG = nx.DiGraph()
        nodes = []
        node2outdegree = {}
        sub2pro = defaultdict(list)
        for member in root.findall('reaction'):
            rank = member.get('id')
            name = member.get('name')
            type = member.get('type')
            substrates_tmp = member.findall('substrate')
            substrates =[]
            for substrate in substrates_tmp:
                substrate = substrate.get('name')
                if 'cpd:' not in substrate:
                    continue
                substrates.append(substrate)
            products_tmp = member.findall('product')
            products = []
            for product in products_tmp:
                product = product.get('name')
                if 'cpd:' not in product:
                    continue
                products.append(product)
            for substrate in substrates:
                for product in products:
                    sub2pro[substrate].append(product)
                    DG.add_edge(substrate, product)
                    nodes.append(substrate)
                    nodes.append(product)
                    network_info = (substrate, product, type)
                    networks_list.append(network_info)
        nodes_num = len(set(nodes))
        xml_2_nodes[xml_list] = nodes
        xml_2_network[xml_list] = networks_list
        all_outdegree = 0
        for n in list(set(nodes)):
            node2outdegree[n] = DG.out_degree(n)
            all_outdegree += DG.out_degree(n)
        G = nx.Graph(DG)
        betweenness = nx.betweenness_centrality(G, normalized=True)
        print(betweenness)
        all_betweenness = sum(value for value in betweenness.values() if betweenness.values())
        if all_betweenness == 0:
            all_betweenness = 0.1
        nodes_info = []
        for n in list(set(nodes)):
            if len(list(set(nodes))):
                rbc = betweenness[n] / all_betweenness
                rod = node2outdegree[n] / all_outdegree
                node_info = (n, rbc, rod)
                nodes_info.append(node_info)
        xml_2_imp[xml_list] = nodes_info
    return xml_2_network, xml_2_imp, xml_2_nodes

def main():
    org = args.org
    cpd_network, cpd_imp, cpd_node = cal_impactvalue(org)
    with open(db_dir + org + "/imp.json", "w") as imp_w, open(db_dir + org + "/network.json", "w") as net_w, open(db_dir + org + "/nodes.json", "w") as nod_w:
        json.dump(cpd_node, nod_w)
        json.dump(cpd_imp, imp_w)
        json.dump(cpd_network, net_w)
        print("加载网络和重要性得分文件完成")
main()
