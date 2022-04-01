# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from collections import OrderedDict

class GcContent(ApiBase):
    def __init__(self, bind_object):
        super(GcContent, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_gc_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "gc_content"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='gc_content',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('gc_content', [main_info])
        else:
            main_id = ObjectId(main_id)

    def add_gc_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "gc_content"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='gc_content',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('gc_content', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        table1_dict = {
            "column": [{"field":"sample_name","filter":"false","sort":"false","title":"sample name","type":"string"},
                       {"field":"seq_count","filter":"false","sort":"false","title":"seq count","type":"int"},
                       {"field":"gc_content","filter":"false","sort":"false","title":"GC content (%)","type":"float"},
                       {"field":"n_ratio","filter":"false","sort":"false","title":"N ratio(%)","type":"float"},
                       {"field":"atgc_ratio","filter":"false","sort":"false","title":"AT/GC ratio","type":"float"}],
                        "condition": {}}
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        table2_dict = {
            "column": [{"field":"sample_name","filter":"false","sort":"false","title":"sample name","type":"string"},
                       {"field":"scaffold_name","filter":"false","sort":"false","title":"scaffold name","type":"string"},
                       {"field":"length","filter":"false","sort":"false","title":"length (bp)","type":"int"},
                       {"field":"gc_content","filter":"false","sort":"false","title":"GC content (%)","type":"float"},
                       {"field":"n_ratio","filter":"false","sort":"false","title":"N ratio (%)","type":"float"},
                       {"field":"atgc_ratio","filter":"false","sort":"false","title":"AT/GC ratio","type":"float"}],
                        "condition": {}}
        table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
        table3_dict = OrderedDict()
        table3_dict["name"] = "scaffold_name"
        table3_dict["data"] = "gc_content"
        table3_dict["condition"] = {}
        table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
        table4_dict = OrderedDict()
        table4_dict["name"] = "scaffold_name"
        table4_dict["y"] = "length"
        table4_dict["x"] = "scaffold_name"
        table4_dict["condition"] = {}
        table4_info = json.dumps(table4_dict, sort_keys=False, separators=(',', ':'))

        insert_data1 = [];insert_data2 = []
        for file in os.listdir(file_path):
            if file.endswith("all.result.xls"):
                result_file = file
                with open(file_path + '/' + result_file, "r") as f:
                    lines = f.readlines()
                    data1 = lines[1:]
                    for j in data1:
                        tmp_data = j.strip("\n").split("\t")
                        insert_data1.append({
                            "gc_id": main_id,
                            "sample_name": tmp_data[0],
                            "seq_count": tmp_data[1],
                            "gc_content": tmp_data[2],
                            "n_ratio": tmp_data[3],
                            "atgc_ratio": tmp_data[4],
                        })
            elif file.endswith("detail.result.xls"):
                detail_file = file
                with open(file_path + '/' + detail_file, "r") as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        data = line.strip("\n").split("\t")
                        insert_data2.append({
                            "gc_id": main_id,
                            "sample_name": data[0],
                            "scaffold_name": data[1],
                            "length": int(data[2]),
                            "gc_content": float(data[3]),
                            "n_ratio": float(data[4]),
                            "atgc_ratio": float(data[5]),
                        })
                insert_data2 = sorted(insert_data2, key=lambda x: x["gc_content"],reverse = True)
                for i in range(len(insert_data2)):
                    insert_data2[i]["scaffold_number"] = i + 1
                try:
                    self.create_db_table('gc_content_detail', insert_data2)
                except Exception as e:
                    self.bind_object.logger.error("导入gc_content_detail数据出错:%s" % e)
                else:
                    self.bind_object.logger.info("导入gc_content_detail数据成功")
                insert_data2 = []
        try:
            self.create_db_table('gc_content_stat', insert_data1)
            self.update_db_record('gc_content', main_id, status="end", main_id=main_id,
                                  stat_table=table1_info, detail_table=table2_info, column_data = table3_info, line_data = table4_info)
        except Exception as e:
            self.bind_object.logger.error("导入geneclust_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入geneclust_detail数据成功")

