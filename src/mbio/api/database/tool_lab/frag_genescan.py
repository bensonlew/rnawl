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

class FragGenescan(ApiBase):
    def __init__(self, bind_object):
        super(FragGenescan, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_frag_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "FragGenescan"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='FragGenescan',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('frag_genescan', [main_info])
        else:
            main_id = ObjectId(main_id)

            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name", "type": "string"},
                    {"field": "gene_no", "filter": "false", "sort": "false", "title": "Gene No.", "type": "int"},
                    {"field": "gene_len", "filter": "false", "sort": "false", "title": "Gene Total Len（bp）", "type": "int"},
                    {"field": "gene_average_len", "filter": "false", "sort": "false", "title": "Gene Average Len（bp）", "type": "float"},
                    {"field": "gc_content", "filter": "false", "sort": "false", "title": "GC Content  (%)","type": "float"},
                    {"field": "gene_ratio", "filter": "false", "sort": "false", "title": "Gene/Genome (%)","type": "float"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "gene_id", "filter": "false", "sort": "false", "title": "Gene ID", "type": "string"},
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name", "type": "string"},
                    {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                    {"field": "strand", "filter": "false", "sort": "false", "title": "Strand", "type": "string"},
                    {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
                    {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                    {"field": "gene_len", "filter": "false", "sort": "false", "title": "Gene Len（bp）", "type": "int"},
                    {"field": "protein_len", "filter": "false", "sort": "false", "title": "Protein Len", "type": "int"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('frag_genescan', main_id, status="end", main_id=main_id,stat_table=table1_info,
                                      detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入frag_genescan数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入frag_genescan数据成功")
        return main_id

    def add_frag_detail(self, main_id, detail_path, stat_file ,project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "FragGenescan"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='FragGenescan',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('frag_genescan', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        with open(stat_file, "r") as f,open(detail_path,"r") as v:
            data1 = f.readlines()
            for i in data1[1:]:
                tmp = i.strip().split("\t")
                insert_data1.append({
                    "frag_genescan_id": main_id,
                    "sample_name": tmp[0],
                    "gene_no": "-" if tmp[1] == "-" else int(tmp[1]),
                    "gene_len": "-" if tmp[2] == "-" else int(tmp[2]),
                    "gene_average_len": "-" if tmp[3] == "-" else float(tmp[3]),
                    "gc_content": "-" if tmp[4] == "-" else float(tmp[4]),
                    "gene_ratio": "-" if tmp[5] == "-" else float(tmp[5]),
                })
            data2 = v.readlines()
            for x in data2[1:]:
                tmp = x.strip().split("\t")
                insert_data2.append({
                    "frag_genescan_id": main_id,
                    "gene_id": tmp[0],
                    "sample_name": tmp[1],
                    "location": tmp[2],
                    "strand": tmp[3],
                    "start": int(tmp[4]) ,
                    "end": int(tmp[5]) ,
                    "gene_len": int(tmp[6]) ,
                    "protein_len": int(tmp[7]) ,
                })

        try:
            self.create_db_table('frag_genescan_stat', insert_data1)
            self.create_db_table('frag_genescan_detail', insert_data2)
        except Exception as e:
            self.bind_object.logger.error("导入frag_genescan数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入frag_genescan数据成功")
        return main_id
