# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.api.database.base import Base
from biocluster.config import Config
import os
import argparse
import re
import json
from collections import defaultdict

class PathwayNetwork(object):
    def __init__(self, json_file, output_dir):
        json_dict = json.load(open(json_file,"r"))
        self.prepare_file = json_dict["prepare_file"]
        self.anno_type = json_dict["anno_type"]
        self.sup_connect = json_dict["sup_connect"]
        self.sup_all = json_dict["sup_all"]
        self.anno_dir = json_dict["anno_dir"]
        self.output_dir = output_dir




    def run_pathway_network(self):
        if self.anno_type.lower() == "kegg":
            nodes_id_infos, link_infos = self.run_kegg_pathway_network()
        else:
            nodes_id_infos, link_infos = self.run_go_pathway_network()
        final_dict = {
            "nodes_id_infos" : nodes_id_infos,
            "link_infos" : link_infos
        }
        with open(os.path.join(self.output_dir , "network_result.json"), 'w') as json_f:
            json.dump(final_dict, json_f, sort_keys=True, indent=4)


    def run_kegg_pathway_network(self):
        raw_df = pd.read_table(self.prepare_file, index_col=0)
        raw_dict = raw_df.to_dict("index")
        raw_ids = set(raw_dict.keys())
        pathway_des_file = os.path.join(os.path.join(self.anno_dir, "kegg_pathway_releted_infos.json"))
        # path_related_file = os.path.join(os.path.join(self.anno_dir, "final_id_relate.json"))
        path_related_file = os.path.join(os.path.join(self.anno_dir, "id_related.json"))
        path_related_dict = json.load(open(path_related_file))
        pathway_des_dict = json.load(open(pathway_des_file))

        if not self.sup_connect and not self.sup_all :
            nodes_ids = set()
            connect_links = set()
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],"pvalue" : raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
            for id in raw_ids:
                if path_related_dict.get(id):
                    if  len(set(path_related_dict.get(id)) & raw_ids) >0:
                            target_ids = set(path_related_dict.get(id)) & raw_ids
                            for t_id in target_ids :
                                cmp = (id,t_id)
                                if (id,t_id) not in connect_links and (t_id,id) not in connect_links:
                                    connect_links.add(cmp)
                                    link_info = {
                                        "source" : nodes_id_list.index(id) ,
                                        "target" : nodes_id_list.index(t_id),
                                        "type" : "actual" ,
                                    }
                                    link_infos.append(link_info)

        if  self.sup_connect and not self.sup_all :
            nodes_ids = set()
            connect_links = set()
            reversed_dict = defaultdict(list)
            inter_id_stat = defaultdict(int)
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],"pvalue" : raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
                related_ids = path_related_dict.get(id)
                if related_ids:
                    for r_id in related_ids:
                        reversed_dict[r_id].append(id)
            for id in raw_ids:
                if path_related_dict.get(id):
                    related_ids = path_related_dict.get(id)
                    for t_id in related_ids:
                        if t_id in raw_ids:
                            cmp = (id, t_id)
                            if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                connect_links.add(cmp)
                                link_info = {
                                    "source": nodes_id_list.index(id),
                                    "target": nodes_id_list.index(t_id),
                                    "type": "actual",
                                }
                                link_infos.append(link_info)
                        else:
                            if len(reversed_dict[t_id]) > 1:
                                if t_id not in nodes_id_list:
                                    node_info = {"group": "node", "id": t_id, "name": pathway_des_dict[t_id]["name"],
                                                 "pvalue": 1,
                                                 "size": 0, "type": "supplementary"}
                                    nodes_id_infos.append(node_info)
                                    nodes_id_list.append(t_id)
                                    nodes_ids.add(t_id)
                                for i_id in reversed_dict[t_id]:
                                    cmp = (t_id, i_id)
                                    if (t_id, i_id) not in connect_links and (t_id, i_id) not in connect_links:
                                        connect_links.add(cmp)
                                        link_info = {
                                            "source": nodes_id_list.index(t_id),
                                            "target": nodes_id_list.index(i_id),
                                            "type": "supplementary",
                                        }
                                        link_infos.append(link_info)

        if self.sup_connect and self.sup_all:
            nodes_ids = set()
            inter_nodes = set()
            reversed_dict = defaultdict(list)
            connect_links = set()
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            all_related_ids = set()
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],
                             "pvalue": raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
                related_ids = path_related_dict.get(id)
                if related_ids:
                    for related_id in related_ids:
                        all_related_ids.add(related_id)
            if all_related_ids:
                for r_id in all_related_ids:
                    if r_id not in nodes_ids:
                        node_info = {"group": "node", "id": r_id, "name": pathway_des_dict[r_id]["name"],
                                         "pvalue": 1,
                                         "size": 0, "type": "supplementary"}
                        nodes_id_infos.append(node_info)
                        nodes_id_list.append(r_id)
                        nodes_ids.add(r_id)

            for id in raw_ids:
                if path_related_dict.get(id):
                    related_ids = path_related_dict.get(id)
                    for t_id in related_ids:
                        if t_id in raw_ids:
                            cmp = (id, t_id)
                            if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                connect_links.add(cmp)
                                link_info = {
                                    "source": nodes_id_list.index(id),
                                    "target": nodes_id_list.index(t_id),
                                    "type": "actual",
                                }
                                link_infos.append(link_info)
                        else:
                            if len(reversed_dict[t_id]) > 1:
                                for i_id in reversed_dict[t_id]:
                                    cmp = (t_id, i_id)
                                    if (t_id, i_id) not in connect_links and (t_id, i_id) not in connect_links:
                                        connect_links.add(cmp)
                                        link_info = {
                                            "source": nodes_id_list.index(t_id),
                                            "target": nodes_id_list.index(i_id),
                                            "type": "supplementary",
                                        }
                                        link_infos.append(link_info)

                            else:
                                cmp = (id, t_id)
                                if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                    connect_links.add(cmp)
                                    link_info = {
                                        "source": nodes_id_list.index(id),
                                        "target": nodes_id_list.index(t_id),
                                        "type": "supplementary",
                                    }
                                    link_infos.append(link_info)
        return nodes_id_infos, link_infos


    def run_go_pathway_network(self):
        raw_df = pd.read_table(self.prepare_file, index_col=0)
        raw_dict = raw_df.to_dict("index")
        raw_ids = set(raw_dict.keys())
        pathway_des_file = os.path.join(os.path.join(self.anno_dir, "go_id_releted_infos.json"))
        # path_related_file = os.path.join(os.path.join(self.anno_dir, "final_id_relate.json"))
        path_related_file = os.path.join(os.path.join(self.anno_dir, "id_related.json"))
        path_related_dict = json.load(open(path_related_file))
        pathway_des_dict = json.load(open(pathway_des_file))
        if not self.sup_connect and not self.sup_all:
            nodes_ids = set()
            connect_links = set()
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],
                             "pvalue": raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
            for id in raw_ids:
                if path_related_dict.get(id):
                    if len(set(path_related_dict.get(id)) & raw_ids) > 0:
                        target_ids = set(path_related_dict.get(id)) & raw_ids
                        for t_id in target_ids:
                            cmp = (id, t_id)
                            if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                connect_links.add(cmp)
                                link_info = {
                                    "source": nodes_id_list.index(id),
                                    "target": nodes_id_list.index(t_id),
                                    "type": "actual",
                                }
                                link_infos.append(link_info)
        if self.sup_connect and not self.sup_all:
            nodes_ids = set()
            connect_links = set()
            reversed_dict = defaultdict(list)
            inter_id_stat = defaultdict(int)
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],
                             "pvalue": raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
                related_ids = path_related_dict.get(id)
                if related_ids:
                    for r_id in related_ids:
                        reversed_dict[r_id].append(id)
            for id in raw_ids:
                if path_related_dict.get(id):
                    related_ids = path_related_dict.get(id)
                    for t_id in related_ids:
                        if t_id in raw_ids:
                            cmp = (id, t_id)
                            if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                connect_links.add(cmp)
                                link_info = {
                                    "source": nodes_id_list.index(id),
                                    "target": nodes_id_list.index(t_id),
                                    "type": "actual",
                                }
                                link_infos.append(link_info)
                        else:
                            if len(reversed_dict[t_id]) > 1:
                                node_info = {"group": "node", "id": t_id, "name": pathway_des_dict[t_id]["name"],
                                             "pvalue": 1,
                                             "size": 0, "type": "supplementary"}
                                nodes_id_infos.append(node_info)
                                nodes_id_list.append(t_id)
                                nodes_ids.add(t_id)
                                for i_id in reversed_dict[t_id]:
                                    cmp = (t_id, i_id)
                                    if (t_id, i_id) not in connect_links and (t_id, i_id) not in connect_links:
                                        connect_links.add(cmp)
                                        link_info = {
                                            "source": nodes_id_list.index(t_id),
                                            "target": nodes_id_list.index(i_id),
                                            "type": "supplementary",
                                        }
                                        link_infos.append(link_info)

        if self.sup_connect and self.sup_all:
            nodes_ids = set()
            inter_nodes = set()
            connect_links = set()
            reversed_dict = defaultdict(list)
            nodes_id_infos = []
            nodes_id_list = []
            link_infos = []
            for id in raw_ids:
                node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],
                             "pvalue": raw_dict[id]["padjust"],
                             "size": raw_dict[id]['gene_num'], "type": "actual"}
                nodes_id_infos.append(node_info)
                nodes_id_list.append(id)
                nodes_ids.add(id)
                related_ids = path_related_dict.get(id)
                if related_ids:
                    for r_id in related_ids:
                        if r_id not in nodes_ids:
                            node_info = {"group": "node", "id": id, "name": pathway_des_dict[id]["name"],
                                         "pvalue": 1,
                                         "size": 0, "type": "supplementary"}
                            nodes_id_infos.append(node_info)
                            nodes_id_list.append(id)
                            nodes_ids.add(id)

            for id in raw_ids:
                if path_related_dict.get(id):
                    related_ids = path_related_dict.get(id)
                    for t_id in related_ids:
                        if t_id in raw_ids:
                            cmp = (id, t_id)
                            if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                connect_links.add(cmp)
                                link_info = {
                                    "source": nodes_id_list.index(id),
                                    "target": nodes_id_list.index(t_id),
                                    "type": "actual",
                                }
                                link_infos.append(link_info)
                        else:
                            node_info = {"group": "node", "id": t_id, "name": pathway_des_dict[t_id]["name"],
                                         "pvalue": 1,
                                         "size": 0, "type": "supplementary"}
                            nodes_id_infos.append(node_info)
                            nodes_id_list.append(t_id)
                            nodes_ids.add(t_id)
                            if len(reversed_dict[t_id]) > 1:
                                for i_id in reversed_dict[t_id]:
                                    cmp = (t_id, i_id)
                                    if (t_id, i_id) not in connect_links and (t_id, i_id) not in connect_links:
                                        connect_links.add(cmp)
                                        link_info = {
                                            "source": nodes_id_list.index(t_id),
                                            "target": nodes_id_list.index(i_id),
                                            "type": "supplementary",
                                        }
                                        link_infos.append(link_info)

                            else:
                                cmp = (id, t_id)
                                if (id, t_id) not in connect_links and (t_id, id) not in connect_links:
                                    connect_links.add(cmp)
                                    link_info = {
                                        "source": nodes_id_list.index(id),
                                        "target": nodes_id_list.index(t_id),
                                        "type": "supplementary",
                                    }
                                    link_infos.append(link_info)
        return nodes_id_infos,link_infos





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_file', type=str, metavar="json_file", required=True, help="main analysis info")
    parser.add_argument('-output_dir', type=str, metavar="Output directory", required=True, help="Output directory name")

    args = parser.parse_args()

    json_file = args.json_file
    outdir = args.output_dir

    run = PathwayNetwork(json_file,outdir)
    run.run_pathway_network()
