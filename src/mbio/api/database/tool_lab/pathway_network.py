# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class PathwayNetwork(ApiBase):
    def __init__(self, bind_object):
        super(PathwayNetwork, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_pathway_network(self, main_id, network_json ,legend_title):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        node_dict = {"name": "name", 'id': 'id', 'group': 'group', "condition": {"type": "node"},"size":"size","color_value":"color_value","node_type":"node_type"}
        line_dict = {"name": "name", "value": 'value', "condition": {"type": "link"},"link_type":"link_type"}
        node_data =  json.dumps(node_dict, sort_keys=True, separators=(',', ':'))
        line_data = json.dumps(line_dict, sort_keys=True, separators=(',', ':'))
        a = json.load(open(network_json))

        node_infos = a["nodes_id_infos"]
        link_infos = a["link_infos"]
        node_df = pd.DataFrame(node_infos)
        node_df.rename(columns ={"pvalue":"color_value","type":"node_type"},inplace=True)
        node_df["type"] = "node"
        # node_df["map_id"] = node_df["id"]
        node_df["size"] = node_df["size"].apply(lambda x: int(x))
        node_df["pathway_network_id"] = main_id
        node_df_t = node_df[node_df["node_type"] != "supplementary"]
        # node_df_t = node_df[node_df["color_value"] != 1]
        pvalue_range = [0,max(node_df_t["color_value"])]



        link_df = pd.DataFrame(link_infos)
        link_df.rename(columns ={"type":"link_type"},inplace=True)
        link_df["value"] = 0.5
        link_df["name"] = link_df.index
        link_df["pathway_network_id"] = main_id
        link_df["type"] = "link"

        max_node_num = max(node_df_t["size"])
        min_node_num = min(node_df_t["size"])
        # node_domain_dict = {"data":[0,1]}
        node_domain_dict = {"data":pvalue_range}
        node_domain_data = json.dumps(node_domain_dict, sort_keys=True, separators=(',', ':'))
        # network_size_data_dict = {"data":[0,max_node_num]}
        network_size_data_dict = {"data": [min_node_num, max_node_num]}
        network_size_data = json.dumps(network_size_data_dict, sort_keys=True, separators=(',', ':'))


        legend_title = {"data": [legend_title]}
        legend_tile_data = json.dumps(legend_title, sort_keys=True, separators=(',', ':'))

        inter = pd.qcut([min_node_num,max_node_num],3)
        inter_list1 = inter.categories
        inter_list = [int(inter_list1[i].right) for i in range(3)]
        inter_list = [min_node_num] + inter_list
        network_size_customizeSize_dict = {"data": inter_list}
        network_size_customizeSize_data = json.dumps(network_size_customizeSize_dict, sort_keys=True, separators=(',', ':'))


        node_detail = node_df.to_dict("r")
        link_detail = link_df.to_dict("r")
        self.create_db_table('pathway_network_detail', node_detail)
        self.create_db_table('pathway_network_detail', link_detail)
        self.update_db_record('pathway_network', main_id, main_id=main_id, node_data=node_data,link_data = line_data,
                              network_color_domain=node_domain_data,
                              network_size_domain = network_size_data,
                              network_size_customizeSize = network_size_customizeSize_data,
                              network_color_title= legend_tile_data ,status = "end")




