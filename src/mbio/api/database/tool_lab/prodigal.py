# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class Prodigal(ApiBase):
    def __init__(self, bind_object):
        super(Prodigal, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_prodigal_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "prodigal"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='prodigal',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prodigal', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "gene_no", "filter": "false", "sort": "false", "title": "Gene No.","type": "int"},
                    {"field": "gene_total_len", "filter": "false", "sort": "false", "title": "Gene Total Len（bp）","type": "int"},
                    {"field": "gene_average_len", "filter": "false", "sort": "false", "title": "Gene Average Len（bp）","type": "float"}   ,
                    {"field": "gc_ontent", "filter": "false", "sort": "false", "title": "GC Content  (%)","type": "float"},
                    {"field": "gene_genome", "filter": "false", "sort": "false", "title": "Gene/Genome (%)","type": "float"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "gene_id", "filter": "false", "sort": "false", "title": "Gene ID", "type": "string"},
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                    {"field": "strand", "filter": "false", "sort": "false", "title": "Strand", "type": "string"},
                    {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
                    {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                    {"field": "len", "filter": "false", "sort": "false", "title": "Len (bp)", "type": "int"},
                    {"field": "in_condon", "filter": "false", "sort": "false", "title": "Initiator Codon", "type": "string"},
                    {"field": "ter_condon", "filter": "false", "sort": "false", "title": "Terminator Codon", "type": "string"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('prodigal', main_id, status="end", main_id=main_id, stat_table=table1_info, detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入phigaro数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入phigaro数据成功")
        return main_id

    def add_prodigal_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "prodigal"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='prodigal',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prodigal', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        detail_data = {}

        dir = os.listdir(file_path)[0]
        with open(file_path+"/"+dir+"/Prodigal_Gene prediction_stat.xls","r") as f, open(file_path+"/"+dir+"/Prodigal_Gene prediction_detail.xls","r") as v:
            data1 = f.readlines()
            for i in data1[1:]:
                insert_data1.append({
                    "prodigal_id": main_id,
                    "sample_name": i.strip().split("\t")[0],
                    "gene_no": int(i.strip().split("\t")[1]),
                    "gene_total_len": int(i.strip().split("\t")[2]),
                    "gene_average_len": float(i.strip().split("\t")[3]),
                    "gc_ontent": float(i.strip().split("\t")[4]),
                    "gene_genome": float(i.strip().split("\t")[5]),
                })
            data2 = v.readlines()
            for x in data2[1:]:
                insert_data2.append({
                    "prodigal_id": main_id,
                    "gene_id": x.strip().split("\t")[0],
                    "sample_name": x.strip().split("\t")[1],
                    "location": x.strip().split("\t")[2],
                    "strand": x.strip().split("\t")[3],
                    "start": int(x.strip().split("\t")[4]),
                    "end": int(x.strip().split("\t")[5]),
                    "len": int(x.strip().split("\t")[6]),
                    "in_condon": x.strip().split("\t")[7],
                    "ter_condon": x.strip().split("\t")[8],
                })
        try:
            if insert_data2:
                self.create_db_table('prodigal_detail', insert_data2)
            self.create_db_table('prodigal_stat', insert_data1)
        except Exception as e:
            self.bind_object.logger.error("导入prodigal数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入prodigal数据成功")
        return main_id
