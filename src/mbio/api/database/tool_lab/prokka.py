# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class Prokka(ApiBase):
    def __init__(self, bind_object):
        super(Prokka, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_prokka_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "prokka"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='prokka',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prokka', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "scaffold_number", "filter": "false", "sort": "false", "title": "scaffold number","type": "int"},
                    {"field": "base", "filter": "false", "sort": "false", "title": "bases (bp)","type": "int"},
                    {"field": "n_ratio", "filter": "false", "sort": "false", "title": "N ratio (%)","type": "float"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "predict_type", "filter": "false", "sort": "false", "title": "predict type", "type": "string"},
                    {"field": "predict_number", "filter": "false", "sort": "false", "title": "predict number", "type": "int"},
                    {"field": "bases_percent", "filter": "false", "sort": "false", "title": "bases percent (%)", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('prokka', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入prokka数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入prokka数据成功")
        return main_id

    def add_prokka_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "prokka"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='prokka',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prokka', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1=[];insert_data2=[];detail_data = {}
        with open(file_path) as f:
            data = f.readlines()
            for i in data:
                if i.startswith("contigs"):
                    contig = int(i.split(":")[1].strip().split("\t")[0])
                elif i.startswith("bases"):
                    base = int(i.split(":")[1].strip().split("\t")[0])
                elif i.startswith("N_base"):
                    n_ratio = float(i.strip().split("\t")[1])
                elif i.startswith("sample_name"):
                    sample_name = i.strip().split("\t")[1]
                elif i.startswith("gene"):
                    detail_data["gene"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("CDS"):
                    detail_data["CDS"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("misc_RNA"):
                    detail_data["misc_RNA"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("rRNA"):
                    detail_data["rRNA"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("repeat_region"):
                    detail_data["repeat_region"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("tRNA"):
                    detail_data["tRNA"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
                elif i.startswith("tmRNA"):
                    detail_data["tmRNA"] = [int(i.split(":")[1].strip().split("\t")[0]),int(i.split(":")[1].strip().split("\t")[1])]
            insert_data1.append({
                "sample_name": sample_name,
                "scaffold_number": contig,
                "base": base,
                "n_ratio": n_ratio,
                "prokka_id": main_id
            })
            if "gene" not in detail_data.keys():
                detail_data["gene"] = [0,0]
            if "CDS" not in detail_data.keys():
                detail_data["CDS"] = [0,0]
            if "misc_RNA" not in detail_data.keys():
                detail_data["misc_RNA"] = [0,0]
            if "rRNA" not in detail_data.keys():
                detail_data["rRNA"] = [0,0]
            if "repeat_region" not in detail_data.keys():
                detail_data["repeat_region"] = [0,0]
            if "tRNA" not in detail_data.keys():
                detail_data["tRNA"] = [0,0]
            if "tmRNA" not in detail_data.keys():
                detail_data["tmRNA"] = [0,0]
            for i in detail_data.keys():
                insert_data2.append({
                    "sample_name": sample_name,
                    "predict_type": i,
                    "predict_number": detail_data[i][0],
                    "bases_percent": round((float(detail_data[i][1]) / base) * 100,2),
                    "prokka_id": main_id
                })
        try:
            self.create_db_table('prokka_detail', insert_data2)
            self.create_db_table('prokka_stat', insert_data1)
        except Exception as e:
            self.bind_object.logger.error("导入prokka数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入prokka数据成功")

