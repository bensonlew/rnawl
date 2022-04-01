# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class Quast(ApiBase):
    def __init__(self, bind_object):
        super(Quast, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_quast_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "quast"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='quast',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('quast', [main_info])
        return main_id

    def add_quast_detail(self, main_id, file_path, picture_path, is_ref, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "quast"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='quast',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('quast', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data = []

        with open(file_path,"r") as f:
            data = f.readlines()
            if is_ref:
                for i in data[1:]:
                    insert_data.append({
                        "quast_id": main_id,
                        "sample_name": i.strip().split("\t")[0],
                        "contigs": int(i.strip().split("\t")[1]),
                        "largest_contig": int(i.strip().split("\t")[3]),
                        "len": float(i.strip().split("\t")[2]),
                        "gc_content": float(i.strip().split("\t")[5]),
                        "n50": float(i.strip().split("\t")[6]),
                        "l50": float(i.strip().split("\t")[7]),
                        "n75": float(i.strip().split("\t")[8]),
                        "l75": float(i.strip().split("\t")[9]),
                        "n100": float(i.strip().split("\t")[10]),
                        "fraction": float(i.strip().split("\t")[11]),
                        "duplication": float(i.strip().split("\t")[12]),
                        "max_align": float(i.strip().split("\t")[13]),
                        "all_algin": float(i.strip().split("\t")[14]),
                        "mismatch": float(i.strip().split("\t")[15]),
                        "indels": float(i.strip().split("\t")[16]),
                    })
            else:
                for i in data[1:]:
                    insert_data.append({
                        "quast_id": main_id,
                        "sample_name": i.strip().split("\t")[0],
                        "contigs": int(i.strip().split("\t")[1]),
                        "largest_contig": int(i.strip().split("\t")[3]),
                        "len": float(i.strip().split("\t")[2]),
                        "gc_content": float(i.strip().split("\t")[4]),
                        "n50": float(i.strip().split("\t")[5]),
                        "l50": float(i.strip().split("\t")[6]),
                        "n75": float(i.strip().split("\t")[7]),
                        "l75": float(i.strip().split("\t")[8]),
                        "n100": float(i.strip().split("\t")[9]),
                    })
        if is_ref:
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "contigs", "filter": "false", "sort": "false", "title": "Contigs", "type": "int"},
                    {"field": "largest_contig", "filter": "false", "sort": "false", "title": "Largest contig（bp）", "type": "int"},
                    {"field": "len", "filter": "false", "sort": "false", "title": "Total length（bp）","type": "int"},
                    {"field": "gc_content", "filter": "false", "sort": "false", "title": "GC Content（%）", "type": "float"},
                    {"field": "n50", "filter": "false", "sort": "false", "title": "N50 (bp)","type": "int"},
                    {"field": "l50", "filter": "false", "sort": "false", "title": "L50", "type": "int"},
                    {"field": "n75", "filter": "false", "sort": "false", "title": "N75 (bp)", "type": "int"},
                    {"field": "l75", "filter": "false", "sort": "false", "title": "L75", "type": "int"},
                    {"field": "n100", "filter": "false", "sort": "false", "title": "N/100 kbp", "type": "float"},
                    {"field": "fraction", "filter": "false", "sort": "false", "title": "Genome fraction（%）", "type": "float"},
                    {"field": "duplication", "filter": "false", "sort": "false", "title": "Duplication ratio", "type": "float"},
                    {"field": "max_align", "filter": "false", "sort": "false", "title": "Largest alignment（bp）", "type": "int"},
                    {"field": "all_algin", "filter": "false", "sort": "false", "title": "Total aligned length（bp）", "type": "int"},
                    {"field": "mismatch", "filter": "false", "sort": "false", "title": "Mismatches /100 kbp", "type": "float"},
                    {"field": "indels", "filter": "false", "sort": "false", "title": "Indels/100 kbp", "type": "float"}],
                "condition": {}}
        else:
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "contigs", "filter": "false", "sort": "false", "title": "Contigs", "type": "int"},
                    {"field": "largest_contig", "filter": "false", "sort": "false", "title": "Largest contig（bp）","type": "int"},
                    {"field": "len", "filter": "false", "sort": "false", "title": "Total length（bp）", "type": "int"},
                    {"field": "gc_content", "filter": "false", "sort": "false", "title": "GC Content（%）", "type": "float"},
                    {"field": "n50", "filter": "false", "sort": "false", "title": "N50 (bp)", "type": "int"},
                    {"field": "l50", "filter": "false", "sort": "false", "title": "L50", "type": "int"},
                    {"field": "n75", "filter": "false", "sort": "false", "title": "N75 (bp)", "type": "int"},
                    {"field": "l75", "filter": "false", "sort": "false", "title": "L75", "type": "int"},
                    {"field": "n100", "filter": "false", "sort": "false", "title": "N/100 kbp", "type": "float"},],
                "condition": {}}

        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.create_db_table('quast_stat', insert_data)
            if is_ref:
                self.update_db_record('quast', main_id, status="end", main_id=main_id, stat_table=table1_info,picture=picture_path)
            else:
                self.update_db_record('quast', main_id, status="end", main_id=main_id, stat_table=table1_info,picture=picture_path)
        except Exception as e:
            self.bind_object.logger.error("导入quast数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入quast数据成功")
        return main_id