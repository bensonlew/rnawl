# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd

class Bugbase(Base):
    def __init__(self, bind_object):
        super(Bugbase, self).__init__(bind_object)
        self._project_type = "metaasv"

    #@report_check
    def add_bugbase(self,level_id,otu_id, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "otu_id" : otu_id,
            "status": "end",
            "desc": "正在计算",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_bugbase"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    #@report_check
    def add_detail(self, prediction_file,threshold_file,normalized_file,table_id=None,project_sn='bugbase',task_id='bugbase', params=None):
        collection = self.db["bugbase"]
        if (table_id is None) or (table_id == ""):
            name = "Bugbase" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Bugbase',
                params=params if params else "null",
                status="start",
            )
            #main_id = self.create_db_table('sg_bugbase', [main_info])
            main_id = collection.insert_one(main_info).inserted_id
        else:
            main_id = ObjectId(table_id)

        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        insert_data3 = []
        insert_data1_dict = {}
        insert_data2_dict = {}
        insert_data3_dict = {}
        data = pd.read_table(prediction_file, sep='\t')
        with open(prediction_file) as v:
            a = v.readlines()
            sample_list = a[0].strip().split("\t")[1:]

            for i in range(len(data)):
                insert_data1_dict = {
                    "bugbase_id": main_id,
                    "phenotypes": data["Phenotypes"][i],
                    }
                for x in sample_list:
                    insert_data1_dict[x] = float(data[x][i]) if data[x][i] != "" else 0.0
                insert_data1.append(insert_data1_dict)
        data2 = pd.read_table(threshold_file, sep='\t')
        phenotype_name = data2.columns
        for ii in range(len(phenotype_name)-2):
            insert_data2_dict = {
                "bugbase_id": main_id,
                "phenotypes": phenotype_name[ii+2],
                }
            for xx in range(len(data2)):
                if data2[phenotype_name[ii+2]][xx]:
                    num = 1
                else:
                    num = 0
                insert_data2_dict[str(data2[phenotype_name[0]][xx])] = num
            insert_data2.append(insert_data2_dict)
        sample_list.insert(0, "phenotypes")
        table_dict = {"table_data": sample_list}
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        with open(normalized_file) as v1:
            data3 = v1.readlines()
            for i in data3[1:]:
                insert_data3_dict = {
                    "bugbase_id": main_id,
                    "otu_id": i.strip().split("\t")[0],
                }
                for x in range(len(data3[0].strip().split("\t")[1:])):
                    insert_data3_dict[data3[0].strip().split("\t")[x + 1]] = i.strip().split("\t")[x + 1]
                insert_data3.append(insert_data3_dict)
        try:
            collection1 = self.db["bugbase_detail"]
            collection1.insert_many(insert_data1)
            collection2 = self.db["bugbase_detail_contribution"]
            collection2.insert_many(insert_data2)
            collection3 = self.db["bugbase_detail_normalized"]
            collection3.insert_many(insert_data3)
            main = self.db['bugbase']
            main.update_one({'_id': main_id},{"$set": {"main_id": main_id,"specimen_list": ",".join(sample_list),"table_data":table_info}})
        except Exception as  e:
            self.bind_object.logger.info("导入bugbase 详情数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入bugbased 详情数据成功")

