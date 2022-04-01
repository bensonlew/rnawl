# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class Ani(ApiBase):
    def __init__(self, bind_object):
        super(Ani, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_ani_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "ani"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ani',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_ani', [main_info])
        else:
            main_id = ObjectId(main_id)

    def add_ani_detail(self, main_id, detail_file, cluster_file, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "ani"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ani',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('ani', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        all_sample = [];insert_data1 = [];insert_data2 = []
        table_dict = {"column": [{"field":"sample_name","filter":"false","sort":"false","title":"sample name","type":"string"}], "condition": {}}
        with open(detail_file) as d,open(cluster_file) as c:
            detail = d.readlines()
            row_tree = c.readline().strip()
            for i in detail[0].split("\t")[1:]:
                all_sample.append(i.strip())
                table_dict["column"].append({"field":i.strip(),"filter":"false","sort":"false","title":i.strip(),"type":"float"})
            for k in detail[1:]:
                tmp_dict = {}
                tmp_dict["sample_name"] = k.split("\t")[0]
                tmp_dict["ani_id"] = main_id
                for num in range(len(detail[0].split("\t"))-1):
                    tmp_dict[detail[0].split("\t")[num+1].strip()] = round(float(k.split("\t")[num+1].strip()),2)
                insert_data1.append(tmp_dict)

            insert_data2.append({
                "name": "cluster_tree",
                "direction" : "h",
                "data": row_tree,
                "type": "tree",
                "ani_id": main_id
            })
            insert_data2.append({
                "name": "cluster_tree",
                "direction": "v",
                "data": row_tree,
                "type": "tree",
                "ani_id": main_id
            })
        heatmap_dict = {"name":"sample_name","data":all_sample, "condition": {}}
        table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
        heatmap_info = json.dumps(heatmap_dict, sort_keys=False, separators=(',', ':'))
        cluster_dict = {"name": "cluster_tree","condition":{"type":"tree"}}
        cluster_info = json.dumps(cluster_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.create_db_table('sg_ani_table', insert_data1)
            self.create_db_table('sg_ani_cluster', insert_data2)
            self.update_db_record('sg_ani', main_id, status="end", main_id=main_id,
                                  stat_table=table_info, heatmap_data=heatmap_info, cluster_data = cluster_info)
        except Exception as e:
            self.bind_object.logger.error("导入ani数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入ani数据成功")

