# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class Housegene(ApiBase):
    def __init__(self, bind_object):
        super(Housegene, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_hgene_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "hgene"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='hgene',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('ani', [main_info])
        else:
            main_id = ObjectId(main_id)

    def add_hgene_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "hgene"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='hgene',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('housegene', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = [];insert_data2 = [];table3_dict = []
        stat_file = file_path + "/all.hgene.xls"
        with open(stat_file) as f:
            data = f.readlines()[1:]
        for i in data:
            if i:
                insert_data1.append({
                    "sample_name": i.split("\t")[0],
                    "hgene_num": int(i.strip().split("\t")[1]),
                    "housegene_id": main_id
                })
        files = os.listdir(file_path)
        for file in files:
            print file_path + "/" + file + "/" + file + "_hgene_blast.xls"
            if os.path.exists(file_path + "/" + file + "/" + file + "_hgene_blast.xls"):
                with open(file_path + "/" + file + "/" + file + "_hgene_blast.xls") as f:
                    for i in f.readlines()[1:]:
                        insert_data2.append({
                            "sample_name": file,
                            "hgene": i.strip().split("\t")[5],
                            "location": i.strip().split("\t")[2],
                            "start": int(i.strip().split("\t")[3]),
                            "end": int(i.strip().split("\t")[4]),
                            "identity": float(i.strip().split("\t")[6]),
                            "coverage": round(float(i.strip().split("\t")[7]) * 100,2),
                            "evalue": float(i.strip().split("\t")[8]),
                            "score": float(i.strip().split("\t")[9]),
                            "housegene_id": main_id
                        })

        table1_dict = {
            "column": [
                {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                {"field": "hgene_num", "filter": "false", "sort": "false", "title": "House-keeping gene", "type": "int"}],
            "condition": {}}
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        table2_dict = {
            "column": [
                {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                {"field": "hgene", "filter": "false", "sort": "false", "title": "House-keeping gene", "type": "string"},
                {"field": "location", "filter": "false", "sort": "false", "title": " Location", "type": "string"},
                {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
                {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                {"field": "identity", "filter": "false", "sort": "false", "title": "Identity(%)", "type": "float"},
                {"field": "coverage", "filter": "false", "sort": "false", "title": "Coverage(%)", "type": "float"},
                {"field": "evalue", "filter": "false", "sort": "false", "title": "Evalue", "type": "float"},
                {"field": "score", "filter": "false", "sort": "false", "title": "Score", "type": "float"}],
            "condition": {}}
        table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
        table3_dict.append({
            "field": str(main_id),
            "title": str(main_id)
        })
        table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.create_db_table('housegene_stat', insert_data1)
            self.create_db_table('housegene_detail', insert_data2)
            self.update_db_record('housegene', main_id, status="end", main_id=main_id,
                                  stat_table=table1_info, detail_table=table2_info,housegene_id=table3_info)
        except Exception as e:
            self.bind_object.logger.error("导入housegene数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入housegene数据成功")

