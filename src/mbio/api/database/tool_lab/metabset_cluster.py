# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from collections import OrderedDict
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class MetabsetCluster(ApiBase):
    def __init__(self, bind_object):
        super(MetabsetCluster, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_metabolome_heatmap_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "MetabsetCluster"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='MetabsetCluster',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('phigaro', [main_info])
        else:
            main_id = ObjectId(main_id)
            heatmap_data = {"heatmap_data": [{"name": "content", "data": "name"}],
                            "condition": {"type": "heatmap"}}
            tree_data = {"tree_data": [{"name": "name"}], "condition": {"type": "tree"}}
            heatmap_bar = {"heatmap_bar_data": [{"name": "name"}], "condition": {"type": "heatmap_bar"}}
            heatmap_data_info = json.dumps(heatmap_data, sort_keys=False, separators=(',', ':'))
            tree_data_info = json.dumps(tree_data, sort_keys=False, separators=(',', ':'))
            heatmap_bar_info =json.dumps(heatmap_bar, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('metabolome_heatmap', main_id, status="end", main_id=main_id,
                                      heatmap_data=heatmap_data, tree_data=tree_data, heatmap_bar_data=heatmap_bar_info)
            except Exception as e:
                self.bind_object.logger.error("导入metabolome_heatmap数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入metabolome_heatmap数据成功")
        return main_id

    def add_metabolome_heatmap_detail(self, main_id, sam_tree=None, matab_tree=None, list_file=None, sample_class=None,metab_class=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "MetabsetCluster"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='MetabsetCluster',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('metabolome_heatmap', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        insert_data3 = []
        insert_data4 = []
        insert_data5 = []
        sample_list = []

        with open(list_file) as f11:
            data = f11.readlines()
            name_dict = data[0].strip().split("\t")[2:]
            for x in data[1:]:
                insert_data_dict = {
                    "heatmap_id": main_id,
                    "type": "heatmap",
                    "name": x.strip().split("\t")[0].decode("unicode_escape"),
                    "metab_id": x.strip().split("\t")[1].decode("unicode_escape"),
                }
                for i in range(len(name_dict)):
                    insert_data_dict[name_dict[i]] = x.strip().split("\t")[i+2]
                insert_data1.append(insert_data_dict)

        if sam_tree:
            with open(sam_tree,"r") as f:
                data1 = f.readlines()
                insert_data2.append({
                    "heatmap_id": main_id,
                        "tree_name": "sample_tree",
                    "data": data1[0].decode("unicode_escape"),
                    "type": "tree",
                    "direction": "h",
                })
        if matab_tree:
            with open(matab_tree, "r") as v:
                data2 = v.readlines()
                insert_data3.append({
                    "heatmap_id": main_id,
                    "tree_name": "matab_tree",
                    "data": data2[0].decode("unicode_escape"),
                    "type": "tree",
                    "direction": "v",
                })
        if sample_class:
            with open(sample_class,"r") as vv:
                data3 = vv.readlines()
                sample_dict = {}
                for i in data3[1:]:
                    if i.strip().split("\t")[1] in sample_dict:
                        sample_dict[i.strip().split("\t")[1]].append(i.strip().split("\t")[0])
                    else:
                        sample_dict[i.strip().split("\t")[1]] = [i.strip().split("\t")[0]]
                for xx in sample_dict:
                    for xxx in sample_dict[xx]:
                        insert_data4.append({
                            "heatmap_id": main_id,
                            "level_color": xx,
                            "name": xxx.decode("unicode_escape"),
                            "direction": "h",
                            "type": "heatmap_bar",
                        })
        if metab_class:
            with open(metab_class,"r") as vvv:
                data4 = vvv.readlines()
                for i in data4:
                    for ii in i.strip().split("\t")[1].split(";"):
                        insert_data5.append({
                            "heatmap_id": main_id,
                            "level_color": str(int(i.strip().split("\t")[0])+1),
                            "name": ii.decode("unicode_escape"),
                            "direction": "v",
                            "type": "heatmap_bar",
                        })
        with open(list_file) as f:
            data5 = f.readlines()
            for i in data5[0].strip().split("\t")[1:]:
                sample_list.append(i)
        table1_dict = {}
        table1_dict["column"] = [{"field":"name","filter":"false","sort":"false","title":"Metabolite","type":"string"}]
        for x in sample_list:
            table1_dict["column"].append({"field":x,"filter":"false","sort":"false","title":x,"type":"float"})
        table1_dict["condition"] = {}
        sample_list = sample_list[1:]
        heatmap_data = {"name": "name", "data": sample_list, "condition": {"type": "heatmap"}}
        tree_data = {"name": "tree_name", "condition": {"type": "tree"}}
        heatmap_bar = {"name": "name", "condition": {"type": "heatmap_bar"}}
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        heatmap_data_info = json.dumps(heatmap_data, sort_keys=False, separators=(',', ':'))
        tree_data_info = json.dumps(tree_data, sort_keys=False, separators=(',', ':'))
        heatmap_bar_info = json.dumps(heatmap_bar, sort_keys=False, separators=(',', ':'))
        try:
            self.create_db_table('metabolome_heatmap_detail', insert_data1)
            if insert_data2:
                self.create_db_table('metabolome_heatmap_detail', insert_data2)
            if insert_data3:
                self.create_db_table('metabolome_heatmap_detail', insert_data3)
            if insert_data4:
                self.create_db_table('metabolome_heatmap_classify', insert_data4)
            if insert_data5:
                self.create_db_table('metabolome_heatmap_classify', insert_data5)
            self.update_db_record('metabolome_heatmap', main_id, status="end", main_id=main_id,stat_table=table1_info,
                                  heatmap_data=heatmap_data_info, tree_data=tree_data_info, heatmap_bar_data=heatmap_bar_info)
        except Exception as e:
            self.bind_object.logger.error("导入metabolome_heatmap数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入metabolome_heatmap数据成功")
        return main_id
