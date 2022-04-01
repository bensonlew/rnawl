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

class Tagseq(ApiBase):
    def __init__(self, bind_object):
        super(Tagseq, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_tagseq_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "Tagseq"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Tagseq',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('tagseq', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "tagseq_no", "filter": "false", "sort": "false", "title": "TagSeq Num","type": "int"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                    {"field": "identity", "filter": "false", "sort": "false", "title": "Identity（%）", "type": "float"},
                    {"field": "alignment_len", "filter": "false", "sort": "false", "title": "Alignment Len（bp）", "type": "int"},
                    {"field": "mismatch", "filter": "false", "sort": "false", "title": "Mismatch（bp）", "type": "int"},
                    {"field": "tag_start", "filter": "false", "sort": "false", "title": "TagSeq Start", "type": "int"},
                    {"field": "tag_end", "filter": "false", "sort": "false", "title": "TagSeq End", "type": "int"},
                    {"field": "s_start", "filter": "false", "sort": "false", "title": "Subject Start", "type": "int"},
                    {"field": "s_end", "filter": "false", "sort": "false", "title": "Subject End", "type": "int"},
                    {"field": "evalue", "filter": "false", "sort": "false", "title": "Evalue", "type": "float"},
                    {"field": "bit_score", "filter": "false", "sort": "false", "title": "Bit Score", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('tagseq', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入tagseq数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入tagseq数据成功")
        return main_id

    def add_tagseq_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "tagseq"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='tagseq',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('tagseq', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []

        file = os.listdir(file_path)
        if len(file) > 0:
            sample_name = file[0].split("/")[-1].split("_vs_")[-1].rstrip(".xls")
            with open(file_path + "/" + file[0], "r") as f:
                data1 = f.readlines()
                number = 0
                for i in data1[1:]:
                    number += 1
                    insert_data1.append({
                        "tagseq_id": main_id,
                        "sample_name": sample_name,
                        "location": i.strip().split("\t")[10],
                        "identity": float(i.strip().split("\t")[3]),
                        "alignment_len": int(i.strip().split("\t")[2]),
                        "mismatch": int(i.strip().split("\t")[9]),
                        "tag_start": int(i.strip().split("\t")[7]),
                        "tag_end": int(i.strip().split("\t")[8]),
                        "s_start": int(i.strip().split("\t")[12]),
                        "s_end": int(i.strip().split("\t")[13]),
                        "evalue": float(i.strip().split("\t")[1]),
                        "bit_score": float(i.strip().split("\t")[0]),
                    })
                    print insert_data1
                insert_data2.append({
                    "tagseq_id": main_id,
                    "sample_name": sample_name,
                    "tagseq_no": number,
                })
        try:
            if insert_data2:
                self.create_db_table('tagseq_stat', insert_data2)
            if insert_data1:
                self.create_db_table('tagseq_detail', insert_data1)
        except Exception as e:
            self.bind_object.logger.error("导入tagseq数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入tagseq数据成功")
        return main_id
