# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import os
import json
import pickle
import unittest
import pandas as pd
from collections import OrderedDict
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Antismash(ApiBase):
    def __init__(self, bind_object):
        super(Antismash, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_Antismash(self, file_path, remote_path, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "Antismash"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Antismash',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_antismash', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        table1_dict = {"column": [{"field":"location","filter":"false","sort":"false","title":"location","type":"string"},
                                        {"field":"cluster_id","filter":"false","sort":"false","title":"cluster_id","type":"string"},
                                        {"field":"type","filter":"false","sort":"false","title":"type","type":"string"},
                                        {"field":"start","filter":"false","sort":"false","title":"start","type":"int"},
                                        {"field":"end","filter":"false","sort":"false","title":"end","type":"int"},
                                        {"field":"similar_cluster","filter":"false","sort":"false","title":"similar_cluster","type":"string"},
                                        {"field":"similarity","filter":"false","sort":"false","title":"similarity","type":"int"},
                                        {"field": "gene_num", "filter": "false", "sort": "false", "title": "gene_num","type": "int"},
                                        {"field":"path","filter":"false","sort":"false","title":"path","type":"string"},
                                        {"field":"predicted_structure","filter":"false","sort":"false","title":"predicted_structure","type":"string"}],
                                        "condition":{}}
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        """
        table2_dict = {"gene_data": [{"field":"gene_id","filter":"false","sort":"false","title":"gene_id","type":"string"},
                                        {"field":"gene_des","filter":"false","sort":"false","title":"gene_des","type":"string"},
                                        {"field":"location","filter":"false","sort":"false","title":"location","type":"string"},
                                        {"field":"start","filter":"false","sort":"false","title":"start","type":"int"},
                                        {"field":"end","filter":"false","sort":"false","title":"end","type":"int"},
                                        {"field":"cluster_id","filter":"false","sort":"false","title":"cluster_id","type":"string"},
                                        {"field":"gene_strand","filter":"false","sort":"false","title":"gene_strand","type":"string"},
                                        {"field":"smcog","filter":"false","sort":"false","title":"smcog","type":"string"},
                                        {"field":"type_detail","filter":"false","sort":"false","title":"type_detail","type":"string"}],
                        "condition":{'type': "table"}}
        table2_info = json.dumps(table2_dict, sort_keys=True, separators=(',', ':'))
        """
        #table3_dict = {"name": "gene_id", "id": "gene_id", "start": "start", "end": "end", "group": "type", "condition": {}}
        #table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
        table3_dict = OrderedDict()
        table3_dict["name"] = "gene_id"
        table3_dict["id"] = "vog_id"
        table3_dict["start"] = "start"
        table3_dict["end"] = "end"
        table3_dict["group"] = "type"
        table3_dict["condition"] = {}
        table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))

        # antismash注释主表注释统计表
        antismash_anno = os.path.join(file_path, 'antismash_anno.xls')
        gene_antismash = os.path.join(file_path, 'gene_antismash.xls')
        gene_list = os.path.join(file_path, 'gene_list.xls')
        insert_data1 = []; insert_data2 = []
        with open(antismash_anno, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)-1):
                Antismash = lines[i+1].strip().split("\t")
                cluster = "Cluster" + Antismash[0].split("_Cluster")[1]
                Scaffold = Antismash[0].split("_Cluster")[0]
                insert_data1.append({
                    "antismash_id": main_id,
                    "location": Scaffold,
                    "cluster_id": cluster,
                    "type": Antismash[1],
                    "start": Antismash[2],
                    "end": Antismash[3],
                    "similar_cluster": Antismash[7],
                    "similarity": Antismash[8],
                    "genes": Antismash[5],
                    "gene_num": Antismash[4],
                    "predicted_structure":Antismash[6],
                    "path":remote_path
                })
        self.create_db_table('sg_antismash_stat', insert_data1)
        # 基因注释表
        with open(gene_antismash, "r") as a, open(gene_list, "r") as b, open(antismash_anno, "r") as f:
            lines1 = a.readlines()[1:]
            lines2 = b.readlines()[1:]
            lines = f.readlines()[1:]
            for i in range(len(lines1)):
                for x in range(len(lines2)):
                    if lines1[i].strip().split("\t")[0] == lines2[x].strip().split("\t")[0]:
                        c = lines1[i].strip().split("\t")
                        d = lines2[x].strip().split("\t")
                        for k in range(len(lines)):
                            if c[0] in lines[k].strip().split("\t")[9]:
                                d[5] = "core"
                        insert_data2.append({
                            "antismash_id": main_id,
                            "gene_id": c[0],
                            "gene_des": c[1],
                            "location": c[3],
                            "start": d[1],
                            "end": d[2],
                            "cluster_id": c[4],
                            "gene_strand": d[3],
                            "smcog": d[4],
                            "type": d[5],
                        })
        self.create_db_table('sg_antismash_detail', insert_data2)
        self.update_db_record('sg_antismash', main_id, status="end", main_id=main_id, antismash_stat_table=table1_info, arrows_axis_data=table3_info)
        return main_id

    def add_Antismash2(self, main_id):
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        table1_dict = {}
        table1_info = json.dumps(table1_dict)
        table3_dict = OrderedDict()
        table3_dict["name"] = "gene_id"
        table3_dict["id"] = "vog_id"
        table3_dict["start"] = "start"
        table3_dict["end"] = "end"
        table3_dict["group"] = "type"
        table3_dict["condition"] = {}
        table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
        self.update_db_record('sg_antismash', main_id, status="end", main_id=main_id, antismash_stat_table=table1_info, arrows_axis_data=table3_info)