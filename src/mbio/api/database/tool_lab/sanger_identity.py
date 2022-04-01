# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import math


class SangerIdentity(ApiBase):
    def __init__(self, bind_object):
        super(SangerIdentity, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_main(self, file1, file2=None,project_sn='tool_lab', method = "false",main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "sanger_identity"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='sanger_identity',
                params= params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('volcano', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)

        if method == "false":
            table1_dict = {
                "column": [{"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name",
                            "type": "string"},
                           {"field": "f_quality", "filter": "false", "sort": "false", "title": "F Quality",
                            "type": "string"},
                           {"field": "f_raw_len", "filter": "false", "sort": "false", "title": "F Raw Len (bp)",
                            "type": "int"},
                           {"field": "f_trim_len", "filter": "false", "sort": "false", "title": "F Trim Len (bp)",
                            "type": "int"},
                           {"field": "f_identity", "filter": "false", "sort": "false", "title": "F Identity",
                            "type": "float"},
                           {"field": "r_quality", "filter": "false", "sort": "false", "title": "R Quality",
                            "type": "string"},
                           {"field": "r_raw_len", "filter": "false", "sort": "false", "title": "R Raw Len (bp)",
                            "type": "int"},
                           {"field": "r_trim_len", "filter": "false", "sort": "false", "title": "R Trim Len (bp)",
                            "type": "int"},
                           {"field": "r_identity", "filter": "false", "sort": "false", "title": "R Identity",
                            "type": "float"},],
                "condition": {'type': 'table'}}
        else:
            table1_dict = {
                "column": [{"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name",
                            "type": "string"},
                           {"field": "f_quality", "filter": "false", "sort": "false", "title": "F Quality",
                            "type": "string"},
                           {"field": "f_raw_len", "filter": "false", "sort": "false", "title": "F Raw Len (bp)",
                            "type": "int"},
                           {"field": "f_trim_len", "filter": "false", "sort": "false", "title": "F Trim Len (bp)",
                            "type": "int"},
                           {"field": "f_identity", "filter": "false", "sort": "false", "title": "F Identity",
                            "type": "float"},
                           {"field": "r_quality", "filter": "false", "sort": "false", "title": "R Quality",
                            "type": "string"},
                           {"field": "r_raw_len", "filter": "false", "sort": "false", "title": "R Raw Len (bp)",
                            "type": "int"},
                           {"field": "r_trim_len", "filter": "false", "sort": "false", "title": "R Trim Len (bp)",
                            "type": "int"},
                           {"field": "r_identity", "filter": "false", "sort": "false", "title": "R Identity(%)",
                            "type": "float"},
                           {"field": "nt_blast", "filter": "false", "sort": "false", "title": "NT Blast",
                            "type": "string"}],
                "condition": {'type': 'table'}}
        table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
        if method == "false":
            pass
        else:
            table2_dict = {"column": [
                {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name", "type": "string"},
                {"field": "subject", "filter": "false", "sort": "false", "title": "Subject", "type": "string"},
                {"field": "taxonomy", "filter": "false", "sort": "false", "title": "Taxonomy", "type": "string"},
                {"field": "identity", "filter": "false", "sort": "false", "title": "Identity(%)", "type": "float"},
                {"field": "align_len", "filter": "false", "sort": "false", "title": "Align Len (bp)", "type": "int"},
                {"field": "mismatches", "filter": "false", "sort": "false", "title": "Mismatches", "type": "int"},
                {"field": "gap_open", "filter": "false", "sort": "false", "title": "Gap_Open", "type": "int"},
                {"field": "q_start", "filter": "false", "sort": "false", "title": "Q_Start", "type": "int"},
                {"field": "q_end", "filter": "false", "sort": "false", "title": "Q_End", "type": "int"},
                {"field": "s_start", "filter": "false", "sort": "false", "title": "S_Start", "type": "int"},
                {"field": "s_end", "filter": "false", "sort": "false", "title": "S_End", "type": "int"},
                {"field": "evalue", "filter": "false", "sort": "false", "title": "Evalue", "type": "float"},
                {"field": "bit_score", "filter": "false", "sort": "false", "title": "Bit Score", "type": "int"}],
                           "condition": {'type': 'table'}}
            table2_info = json.dumps(table2_dict, sort_keys=True, separators=(',', ':'))

        insert_data1 = [];insert_data2 =[]
        with open(file1) as f:
            data = f.readlines()
            for i in data[1:]:
                tmp = i.strip().split("\t")
                insert_data1.append({
                    "si_id": main_id,
                    "sample_name": tmp[0],
                    "f_quality": tmp[1],
                    "f_raw_len": tmp[2],
                    "f_trim_len": tmp[3],
                    "f_identity": tmp[4],
                    "r_quality": tmp[5],
                    "r_raw_len": tmp[6],
                    "r_trim_len": tmp[7],
                    "r_identity": tmp[8],
                    "nt_blast": tmp[9] if len(tmp) > 9 else "",
                    "type":"table"
                })
        if file2 is None:
            pass
        else:
            with open(file2) as v:
                data = v.readlines()
                for i in data[1:]:
                    tmp = i.strip().split("\t")
                    insert_data2.append({
                        "si_id": main_id,
                        "sample_name": tmp[0],
                        "subject": tmp[1],
                        "taxonomy": tmp[2],
                        "identity": tmp[3],
                        "align_len": tmp[4],
                        "mismatches": tmp[5],
                        "gap_open": tmp[6],
                        "q_start": tmp[7],
                        "q_end": tmp[8],
                        "s_start": tmp[9],
                        "s_end": tmp[10],
                        "evalue": tmp[11],
                        "bit_score": tmp[12],
                        "type": "table"
                    })
        try:
            self.create_db_table('sanger_identity_stat', insert_data1)
            if insert_data2:
                self.create_db_table('sanger_identity_detail', insert_data2)
            if method == "false":
                self.update_db_record('sanger_identity', main_id, status="end", main_id=main_id, stat_table=table1_info)
            else:
                self.update_db_record('sanger_identity', main_id, status="end", main_id=main_id, stat_table=table1_info,detail_table=table2_info)
        except Exception as e:
            self.bind_object.logger.error("导入sanger_identity数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入sanger_identity数据成功")
