# -*- coding: utf-8 -*-
# __author__ = 'JieYao'

import argparse
import os

import networkx
import pandas as pd

name_list = ["0"]


def search(node_name):
    global name_list
    for i in range(len(name_list)):
        if node_name == name_list[i]:
            return i
    name_list += [node_name]
    return len(name_list) - 1


parser = argparse.ArgumentParser(description='输入蛋白质相互作用网络，输出网络信息')
parser.add_argument('-i', "--PPI_network", help="输入的PPI网络", required=True)
parser.add_argument('-c', "--cut", help='蛋白相互作用阈值', required=False)
parser.add_argument('-o', "--output", help="输出文件输出路径", required=True)
parser.add_argument('-n', "--nodes", help="节点边文件输入路径", required=True)

args = vars(parser.parse_args())

inFile = args["PPI_network"]
outFile = args["output"]
nodeFile = args["nodes"]

if not os.path.exists(outFile):
    os.mkdir(outFile)
if not args["cut"]:
    cut = -1
else:
    cut = eval(args["cut"])

nodes_table = pd.read_table(nodeFile, header=0, usecols=[0,1,2,3,4,5])
nodes_table.sort_values(['degree'],  ascending = [False], inplace=True)
nodes_table = nodes_table.drop_duplicates("accession_id")
nodes_table['node_id'] = range(0, len(nodes_table))
id2acc =dict(zip(nodes_table['node_id'],nodes_table['node'].map(str)))
acc2id =dict(zip(nodes_table['node'].map(str),nodes_table['node_id']))
id2acc2 =dict(zip(nodes_table['node_id'],nodes_table['accession_id'].map(str)))

G = networkx.Graph()
name_list = ["0"]
with open(inFile, "r") as tmp_file:
    data = tmp_file.readlines()
for i in range(1, len(data)):
    s = data[i].rstrip().split("\t")
    if eval(s[15]) >= cut:
        if s[16] in acc2id and s[17] in acc2id:
            G.add_edge(acc2id[s[16]], acc2id[s[17]], weight = eval(s[15]))

Transitivity = networkx.transitivity(G)
Clustering = networkx.clustering(G)
try:
    Degree_distribution = networkx.degree_histogram(G)
    Degree_Centrality = networkx.degree_centrality(G)
    Closeness_Centrality = networkx.closeness_centrality(G)
    Betweenness_Centrality = networkx.betweenness_centrality(G)
except:
    Degree_distribution = []
    Degree_Centrality = []
    Closeness_Centrality = []
    Betweenness_Centrality = []

with open(os.path.join(args["output"], "protein_interaction_network_degree_distribution.txt"), "w") as tmp_file:
    tmp_file.write("Degree\tNode_Num\n")
    for i in range(len(Degree_distribution)):
        tmp_file.write(str(i) + "\t" + str(Degree_distribution[i]) + "\n")

with open(os.path.join(args["output"], "protein_interaction_network_node_degree.txt"), "w") as tmp_file:
    tmp_file.write("Node_ID\tNode_Name\tDegree\tAccession_id\n")
    for i in id2acc.keys():
        try:
            print G.degree(i)
            tmp_file.write(str(i)+"\t"+str(id2acc[i])+"\t")
            tmp_file.write(str(G.degree(i)) + "\t" + str(id2acc2[i]) + "\n")
        except:
            pass
with open(os.path.join(args["output"], "protein_interaction_network_centrality.txt"), "w") as tmp_file:
    tmp_file.write("Node_ID\tNode_Name\tDegree_Centrality\t")
    tmp_file.write("Closeness_Centrality\tBetweenness_Centrality\tAccession_id\n")
    for i in id2acc.keys():
        try:
            print G.degree(i)
            tmp_file.write(str(i)+"\t"+str(id2acc[i])+"\t")
            tmp_file.write(str(Degree_Centrality[i])+"\t")
            tmp_file.write(str(Closeness_Centrality[i])+"\t")
            tmp_file.write(str(Betweenness_Centrality[i])+"\t")
            tmp_file.write(str(id2acc2[i])+"\n")
        except:
            pass

# with open(os.path.join(args["output"], "protein_interaction_network_clustering.txt"), "w") as tmp_file:
#     tmp_file.write("Node_ID\tProtein_Name\tClustering\n")
#     for i in id2acc.keys():
#         tmp_file.write(str(i)+"\t"+id2acc[i]+"\t"+str(Clustering[i])+"\n")

# with open(os.path.join(args["output"], "protein_interaction_network_transitivity.txt"), "w") as tmp_file:
#     tmp_file.write("Transitivity\t")
#     tmp_file.write(str(Transitivity)+"\n")
