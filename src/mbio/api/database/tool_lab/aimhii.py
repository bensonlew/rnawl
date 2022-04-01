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

class Aimhii(ApiBase):
    def __init__(self, bind_object):
        super(Aimhii, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_aimhii_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "aimhii"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='aimhii',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('tagseq', [main_info])
        else:
            main_id = ObjectId(main_id)

            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "sample", "type": "string"},
                    {"field": "aim_num", "filter": "false", "sort": "false", "title": "AIM Num", "type": "string"},
                    {"field": "type", "filter": "false", "sort": "false", "title": "type", "type": "string"},
                    {"field": "ref_chrom", "filter": "false", "sort": "false", "title": "ref_chrom", "type": "string"},
                    {"field": "left_start_b1", "filter": "false", "sort": "false", "title": "left_start_b1", "type": "int"},
                    {"field": "left_junction_b1", "filter": "false", "sort": "false", "title": "left_junction_b1", "type": "int"},
                    {"field": "left_numreads", "filter": "false", "sort": "false", "title": "left_numreads", "type": "int"},
                    {"field": "gap_length", "filter": "false", "sort": "false", "title": "gap_length", "type": "int"},
                    {"field": "insert_length", "filter": "false", "sort": "false", "title": "insert_length", "type": "int"},
                    {"field": "insert_chrom", "filter": "false", "sort": "false", "title": "insert_chrom", "type": "string"},
                    {"field": "insert_start", "filter": "false", "sort": "false", "title": "insert_start", "type": "int"},
                    {"field": "insert_end", "filter": "false", "sort": "false", "title": "insert_end", "type": "int"},
                    {"field": "insert_strand", "filter": "false", "sort": "false", "title": "insert_strand", "type": "int"},
                    {"field": "right_junction_b1", "filter": "false", "sort": "false", "title": "right_junction_b1", "type": "int"},
                    {"field": "right_end_b1", "filter": "false", "sort": "false", "title": "right_end_b1", "type": "int"},
                    {"field": "right_numreads", "filter": "false", "sort": "false", "title": "right_numreads", "type": "float"},
                    {"field": "picture", "filter": "false", "sort": "false", "title": "picture", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('aimhii', main_id, status="end", main_id=main_id,
                                      detail_table=table2_info)
            except Exception as e:
                self.bind_object.logger.error("导入aimhii数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入aimhii数据成功")
        return main_id

    def add_aimhii_detail(self, main_id, file_path, picture_path ,project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "aimhii"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='aimhii',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('tagseq', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []

        with open(file_path, "r") as f:
            data1 = f.readlines()
            for i in data1[1:]:
                tmp = i.strip().split(",")
                insert_data1.append({
                    "tagseq_id": main_id,
                    "sample_name": tmp[0],
                    "aim_num": tmp[1],
                    "type": tmp[2],
                    "ref_chrom": tmp[3],
                    "left_start_b1": tmp[4] if tmp[4] else "NA",
                    "left_junction_b1": tmp[5] if tmp[5] else "NA",
                    "left_numreads": tmp[6] if tmp[6] else "NA",
                    "gap_length": tmp[7] if tmp[7] else "NA",
                    "insert_length": tmp[8] if tmp[8] else "NA",
                    "insert_chrom": tmp[9] if tmp[9] else "NA",
                    "insert_start": tmp[10] if tmp[10] else "NA",
                    "insert_end": tmp[11] if tmp[11] else "NA",
                    "insert_strand": tmp[12] if tmp[12] else "NA",
                    "right_junction_b1": tmp[13] if tmp[13] else "NA",
                    "right_end_b1": tmp[14] if tmp[14] else "NA",
                    "right_numreads": tmp[15] if tmp[15] else "NA",
                    "picture": picture_path + "/" + "{}.{}.pdf".format(tmp[0],tmp[1])
                })

        try:
            self.create_db_table('aimhii_detail', insert_data1)
        except Exception as e:
            self.bind_object.logger.error("导入aimhii数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入aimhii数据成功")
        return main_id
