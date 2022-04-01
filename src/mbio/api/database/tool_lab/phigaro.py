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

class Phigaro(ApiBase):
    def __init__(self, bind_object):
        super(Phigaro, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_phigaro_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "phigaro"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='phigaro',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('phigaro', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "prophage_id", "filter": "false", "sort": "false", "title": "Prophage ID", "type": "string"},
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location","type": "string"},
                    {"field": "start", "filter": "false", "sort": "false", "title": "Start","type": "int"},
                    {"field": "end", "filter": "false", "sort": "false", "title": "End","type": "int"},
                    {"field": "len", "filter": "false", "sort": "false", "title": "Len (bp)","type": "int"},
                    {"field": "vog_no", "filter": "false", "sort": "false", "title": "VOG No.","type": "int"},
                    {"field": "taxonomy", "filter": "false", "sort": "false", "title": "Taxonomy","type": "string"},
                    {"field": "gc_ontent", "filter": "false", "sort": "false", "title": "GC Content (%)","type": "float"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "prophage_id", "filter": "false", "sort": "false", "title": "Prophage ID", "type": "string"},
                    {"field": "vog_id", "filter": "false", "sort": "false", "title": "VOG ID", "type": "string"},
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                    {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
                    {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                    {"field": "len", "filter": "false", "sort": "false", "title": "Len (bp)", "type": "int"},
                    {"field": "gc_ontent", "filter": "false", "sort": "false", "title": "GC Content (%)", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            table3_dict = OrderedDict()
            table3_dict["name"] = "vog_id"
            table3_dict["id"] = "vog_id"
            table3_dict["start"] = "start"
            table3_dict["end"] = "end"
            table3_dict["group"] = "prophage_id"
            table3_dict["condition"] = {}
            table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('phigaro', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info,arrows_axis_data=table3_info)
            except Exception as e:
                self.bind_object.logger.error("导入phigaro数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入phigaro数据成功")
        return main_id

    def add_phigaro_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "phigaro"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='phigaro',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('phigaro', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        detail_data = {}

        dir = os.listdir(file_path)[0]
        with open(file_path+"/"+dir+"/phigaro_stat.xls","r") as f, open(file_path+"/"+dir+"/phigaro_detail.xls","r") as v:
            data1 = f.readlines()
            for i in data1[1:]:
                insert_data1.append({
                    "phigaro_id": main_id,
                    "prophage_id": i.strip().split("\t")[0],
                    "sample_name": i.strip().split("\t")[1],
                    "location": i.strip().split("\t")[2],
                    "start": int(i.strip().split("\t")[3]),
                    "end": int(i.strip().split("\t")[4]),
                    "len": int(i.strip().split("\t")[5]),
                    "vog_no": int(i.strip().split("\t")[6]),
                    "taxonomy": i.strip().split("\t")[7],
                    "gc_ontent": float(i.strip().split("\t")[8]),
                })
            data2 = v.readlines()
            for x in data2[1:]:
                insert_data2.append({
                    "phigaro_id": main_id,
                    "prophage_id": x.strip().split("\t")[0],
                    "vog_id": x.strip().split("\t")[1],
                    "sample_name": x.strip().split("\t")[2],
                    "location": x.strip().split("\t")[3],
                    "start": int(x.strip().split("\t")[4]),
                    "end": int(x.strip().split("\t")[5]),
                    "len": int(x.strip().split("\t")[6]),
                    "gc_ontent": float(x.strip().split("\t")[7]),
                })
            insert_data2 = sorted(insert_data2, key=lambda k: k['start'])
        try:
            self.create_db_table('phigaro_detail', insert_data2)
            self.create_db_table('phigaro_stat', insert_data1)
        except Exception as e:
            self.bind_object.logger.error("导入phigaro数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入phigaro数据成功")
        return main_id
