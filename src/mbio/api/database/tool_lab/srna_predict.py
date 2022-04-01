# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class SrnaPredict(ApiBase):
    def __init__(self, bind_object):
        super(SrnaPredict, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_srna_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "srna_predict"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='srna_predict',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('srna_predict', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "srna_num", "filter": "false", "sort": "false", "title": "sRNA Num","type": "int"},
                    {"field": "average_length", "filter": "false", "sort": "false", "title": "Average Length（bp）","type": "int"},
                    {"field": "percent", "filter": "false", "sort": "false", "title": "Percent（%）","type": "float"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "id", "filter": "false", "sort": "false", "title": "ID", "type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                    {"field": "acc_num", "filter": "false", "sort": "false", "title": "Acc Num", "type": "string"},
                    {"field": "family", "filter": "false", "sort": "false", "title": "Family", "type": "string"},
                    {"field": "description", "filter": "false", "sort": "false", "title": "Description", "type": "string"},
                    {"field": "srna_start", "filter": "false", "sort": "false", "title": "sRNA Start", "type": "int"},
                    {"field": "srna_end", "filter": "false", "sort": "false", "title": "sRNA End", "type": "int"},
                    {"field": "srna_length", "filter": "false", "sort": "false", "title": "sRNA Length（bp）", "type": "int"},
                    {"field": "strand", "filter": "false", "sort": "false", "title": "Strand", "type": "string"},
                    {"field": "evalue", "filter": "false", "sort": "false", "title": "E-value", "type": "float"},
                    {"field": "score", "filter": "false", "sort": "false", "title": "Score", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('srna_predict', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入housegene数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入housegene数据成功")
        return main_id

    def add_srna_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "srna_predict"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='srna_predict',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('srna_predict', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []; insert_data2 = []; num = 0; srna_length = 0
        for file in os.listdir(file_path):
            if file != "length_stat.txt":
                sample_name = "_".join(file.split("_")[:-1])
                break
        with open(file_path + '/' + sample_name + '_assembly.tblout') as f,open(file_path + "/length_stat.txt") as v:
            data2 = v.readlines()
            all_length = int(data2[0])
            data = f.readlines()
            for i in data:
                tmp = re.split(r'\s+',i)
                if not tmp[0].startswith("#"):
                    num += 1
                    if len(str(num)) == 1:
                        srna_id = "sRNA00" + str(num)
                    elif len(str(num)) == 2:
                        srna_id = "sRNA0" + str(num)
                    elif len(str(num)) == 3:
                        srna_id = "sRNA" + str(num)
                    insert_data1.append({
                        "sample_name": sample_name,
                        "id": srna_id,
                        "location": tmp[3],
                        "acc_num": tmp[2],
                        "family": tmp[1],
                        "description": " ".join(tmp[26:]).rstrip(),
                        "srna_start": int(tmp[9]),
                        "srna_end": int(tmp[10]),
                        "srna_length":  abs(int(tmp[10]) - int(tmp[9])) + 1,
                        "strand": tmp[11],
                        "evalue": float(tmp[17]),
                        "score": round(float(tmp[16]),1),
                        "spre_id": main_id
                    })
                    srna_length += (abs(int(tmp[10]) - int(tmp[9])) + 1)
            if num ==0:
                insert_data2.append({
                    "sample_name": sample_name,
                    "srna_num": 0,
                    "average_length": "-",
                    "Percent": "-",
                    "spre_id": main_id
                })
            else:
                insert_data2.append({
                    "sample_name": sample_name,
                    "srna_num": num,
                    "average_length": srna_length/num,
                    "percent": round(float(srna_length) / all_length * 100,2),
                    "spre_id": main_id
                })



        try:
            self.create_db_table('srna_predict_detail', insert_data1)
            self.create_db_table('srna_predict_stat', insert_data2)
        except Exception as e:
            self.bind_object.logger.error("导入housegene数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入housegene数据成功")

