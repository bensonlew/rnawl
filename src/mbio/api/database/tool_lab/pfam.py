# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
import os, re
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from collections import OrderedDict

class Pfam(ApiBase):
    def __init__(self, bind_object):
        super(Pfam, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_pfam_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "pfam"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='pfam',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('prokka', [main_info])
        else:
            main_id = ObjectId(main_id)
            table1_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name","type": "string"},
                    {"field": "seq_count", "filter": "false", "sort": "false", "title": "seq count","type": "int"},
                    {"field": "pfam_num", "filter": "false", "sort": "false", "title": "Pfam No.","type": "int"}],
                "condition": {}}
            table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
            table2_dict = {
                "column": [
                    {"field": "sample_name", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                    {"field": "seq_id", "filter": "false", "sort": "false", "title": "seq id", "type": "string"},
                    {"field": "pfam", "filter": "false", "sort": "false", "title": "Pfam ID", "type": "string"},
                    {"field": "domain", "filter": "false", "sort": "false", "title": "Domain", "type": "string"},
                    {"field": "type", "filter": "false", "sort": "false", "title": "type", "type": "string"},
                    {"field": "clan_id", "filter": "false", "sort": "false", "title": "Clan ID", "type": "string"},
                    {"field": "domain_description", "filter": "false", "sort": "false", "title": "Domain description", "type": "string"},
                    {"field": "evalue", "filter": "false", "sort": "false", "title": "E-value", "type": "float"}],
                "condition": {}}
            table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
            table3_dict = OrderedDict()
            table3_dict["data"] = "num"
            table3_dict["name"] = "pfam_type"
            table3_dict["condition"] = {"type":"pie"}
            table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
            try:
                self.update_db_record('pfam_predict', main_id, status="end", main_id=main_id,
                                      stat_table=table1_info, detail_table=table2_info, pie_plot=table3_info)
            except Exception as e:
                self.bind_object.logger.error("导入pfam数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入pfam数据成功")
        return main_id

    def add_pfam_detail(self, main_id, file_path, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "pfam"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='pfam',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('pfam_predict', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data1=[];insert_data2=[];insert_data3=[]
        pie_dict = {"Domain":0,"Family":0,"Motif":0,"Repeat":0,"Disordered":0,"Coiled-coil":0,}
        for file in os.listdir(file_path):
            if not file.startswith("stat"):
                with open(file_path + '/' + file) as g,open(file_path + '/stat') as f:
                    tmp1 = f.readlines()
                    data1 = tmp1[0].strip().split("\t")
                    sample_name = data1[0]
                    insert_data1.append({
                        "sample_name": sample_name,
                        "seq_count": int(data1[1]),
                        "pfam_num": int(data1[2]),
                        "pfam_id": main_id
                    })
                    if int(data1[1]) == 0:
                        pass
                    else:
                        tmp2 = g.readlines()
                        for i in tmp2[1:]:
                            data2 = i.strip().split("\t")
                            if data2[7] in pie_dict.keys():
                                pie_dict[data2[7]] += 1
                            insert_data2.append({
                                "sample_name": sample_name,
                                "seq_id": data2[0],
                                "pfam": data2[5],
                                "domain": data2[6],
                                "type": data2[7],
                                "clan_id": data2[14],
                                "domain_description": data2[15],
                                "evalue": float(data2[12]),
                                "pfam_id":main_id
                            })
        if int(data1[1]) == 0:
            pass
        else:
            for i in pie_dict.keys():
                if pie_dict[i] != 0:
                    insert_data3.append({
                        "sample_name": data1[0],
                        "pfam_type": i,
                        "num": pie_dict[i],
                        "type": "pie",
                        "pfam_id": main_id
                    })
        try:
            self.create_db_table('pfam_predict_detail', insert_data2)
            self.create_db_table('pfam_predict_stat', insert_data1)
            self.create_db_table('pfam_predict_plot', insert_data3)
        except Exception as e:
            self.bind_object.logger.error("导入pfam数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入pfam数据成功")

