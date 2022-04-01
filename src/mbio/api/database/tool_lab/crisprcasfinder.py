# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210430

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Crisprcasfinder(ApiBase):
    def __init__(self, bind_object):
        super(Crisprcasfinder, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_crisprcasfinder(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "Crisprcasfinder" + "_"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Crisprcasfinder',
                params=params,
                status="start")
            main_id = self.create_db_table('crisfinder', [main_info])
        else:
            main_id = ObjectId(main_id)

            try:
                self.update_db_record('crisfinder', main_id, )
            except Exception as e:
                self.bind_object.logger.error("导入crisfinder数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入crisfinder数据成功")
        return main_id

    def add_crisprcasfinder_stat(self, main_id, file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "cris_id": ObjectId(main_id),
                    "sample_name": lin[0],
                    "num": lin[1],
                }
                data_list.append(insert_data)
        self.create_db_table('crisfinder_stat', data_list)
        table_dict = {
            "column": [
            {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name", "type": "string"},
            {"field": "num", "filter": "false", "sort": "false", "title": " CRISPR-Cas No.", "type": "string"}],
            "condition":{}}
        table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('crisfinder', main_id, main_id=main_id, stat_table=table_info, status="end")
        except Exception as e:
            self.bind_object.logger.error("导入crisfinder_stat数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入crisfinder_stat数据成功")


    def add_crisprcasfinder_detail(self, main_id, file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "cris_id": ObjectId(main_id),
                    "cong_id": lin[0],
                    "sample_name": lin[1],
                    "location": lin[2],
                    "start": lin[3],
                    "end": lin[4],
                    "cris_len": lin[5],
                    "dr_num": lin[6],
                    "dr_len": lin[7],
                    "spa_num": lin[8],
                    "spa_len": lin[9],
                }
                data_list.append(insert_data)
        self.create_db_table('crisfinder_detail', data_list)
        table_dict = {
            "column": [
                {"field": "cong_id", "filter": "false", "sort": "false", "title": "CRISPR ID", "type": "string"},
                {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample Name", "type": "string"},
                {"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
                {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "string"},
                {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                {"field": "cris_len", "filter": "false", "sort": "false", "title": "CRISPR Length（bp）", "type": "int"},
                {"field": "dr_num", "filter": "false", "sort": "false", "title": "DR Num", "type": "string"},
                {"field": "dr_len", "filter": "false", "sort": "false", "title": "DR Average Len (bp)", "type": "float"},
                {"field": "spa_num", "filter": "false", "sort": "false", "title": "SPA Num", "type": "int"},
                {"field": "spa_len", "filter": "false", "sort": "false", "title": "SPA Average Len (bp)", "type": "int"}],
            "condition": {}}
        table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('crisfinder', main_id, main_id=main_id, detail_table=table_info, status="end")
        except Exception as e:
            self.bind_object.logger.error("导入crisfinder_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入crisfinder_detail数据成功")

    def add_crisprcasfinder_type(self, main_id, file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "cris_id": ObjectId(main_id),
                    "cong_id": lin[0],
                    "element": lin[1],
                    "name": lin[1],
                    "start": int(lin[2]),
                    "strand": "+",
                    "icon": "rect",
                    "end": int(lin[3]),
                    "seq": lin[4],
                    "type": "element_structure",
                }
                data_list.append(insert_data)
        self.create_db_table('crisfinder_type', data_list)
        table_dict = {
            "column": [
                {"field": "cong_id", "filter": "false", "sort": "false", "title": "CRISPR ID", "type": "string"},
                {"field": "name", "filter": "false", "sort": "false", "title": "Type", "type": "string"},
                {"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
                {"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
                {"field": "seq", "filter": "false", "sort": "false", "title": "Sequence", "type": "string"}],
            "condition": {}}
        table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
        element_structure_data = {"name":"name", "element":"element", "icon":"icon", "strand":"strand", "start":"start", "end":"end", "condition": {"type":"element_structure"}}
        element_structure_info = json.dumps(element_structure_data, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('crisfinder', main_id, main_id=main_id, type_table=table_info, element_structure_data = element_structure_info,status="end")
        except Exception as e:
            self.bind_object.logger.error("导入crisfinder_type数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入crisfinder_type数据成功")