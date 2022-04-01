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

class Vfdb(ApiBase):
    def __init__(self, bind_object):
        super(Vfdb, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_vfdb_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "vfdb"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='vfdb',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('vfdb', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "gene_no", "filter": "false", "sort": "false", "title": "Gene Num","type": "int"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "geneid", "filter": "false", "sort": "false", "title": "ID", "type": "string"},
                    {"field": "vfdbid", "filter": "false", "sort": "false", "title": "VFDB ID", "type": "string"},
                    {"field": "vfs", "filter": "false", "sort": "false", "title": "VFs", "type": "string"},
                    {"field": "species", "filter": "false", "sort": "false", "title": "Species", "type": "string"},
                    {"field": "identity", "filter": "false", "sort": "false", "title": "Identity（%）", "type": "float"},
                    {"field": "evalue", "filter": "false", "sort": "false", "title": "E-value", "type": "float"},
                    {"field": "score", "filter": "false", "sort": "false", "title": "Score", "type": "float"},
                ],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('vfdb', main_id, status="end", main_id=main_id, stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入vfdb数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入vfdb数据成功")
        return main_id

    def add_vfdb_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "vfdb"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='vfdb',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('vfdb', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []

        all_dir = os.listdir(file_path)
        anno_file = ""
        for file in os.listdir(file_path):
            if file.endswith("vfdb_anno.xls"):
                anno_file = file_path + "/" + file
        with open(anno_file,"r") as f:
            data1 = f.readlines()
            num = 0
            for i in data1[1:]:
                name = i.strip().split("\t")[2]
                insert_data1.append({
                    "vfdb_id": main_id,
                    "sample_name": i.strip().split("\t")[2],
                    "geneid": i.strip().split("\t")[0],
                    "vfdbid": i.strip().split("\t")[3],
                    "vfs": i.strip().split("\t")[4],
                    "species": i.strip().split("\t")[5],
                    "identity": float(i.strip().split("\t")[9]),
                    "evalue": float(i.strip().split("\t")[10]),
                    "score": float(i.strip().split("\t")[11]),
                })
                num += 1
            insert_data2.append({
                "vfdb_id": main_id,
                "sample_name": name,
                "gene_no": str(num),
            })

        try:
            if insert_data1:
                self.create_db_table('vfdb_detail', insert_data1)
            self.create_db_table('vfdb_stat', insert_data2)
        except Exception as e:
            self.bind_object.logger.error("导入vfdb数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入vfdb数据成功")
        return main_id
