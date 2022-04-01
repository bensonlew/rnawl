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

class Faprotax(Base):
    def __init__(self, bind_object):
        super(Faprotax, self).__init__(bind_object)
        self._project_type = "metaasv"

    #@report_check
    def add_Faprotax(self, prediction_file,group_file,main_id=None,project_sn='faprotax',task_id='faprotax', params=None):
        collection = self.db["faprotax"]
        if (main_id is None) or (main_id == ""):
            name = "Faprotax" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Faprotax',
                params=params if params else "null",
                status="start",
            )
            main_id = collection.insert_one(main_info).inserted_id
        else:
            main_id = ObjectId(main_id)

        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        insert_data1 = []
        insert_data2 = []
        insert_data1_dict = {}
        insert_data2_dict = {}
        data = pd.read_table(prediction_file, sep='\t')
        with open(prediction_file) as v:
            a = v.readlines()
            sample_list = a[0].strip().split("\t")[1:]

            for i in range(len(data)):
                insert_data1_dict = {
                    "faprotax_id": main_id,
                    "functional_groups": data["#group"][i],
                    "type":"sample"
                    }
                all_num1 = 0
                for x in sample_list:
                    insert_data1_dict[x] = float(data[x][i])
                    all_num1 += float(data[x][i])
                insert_data1.append(insert_data1_dict)
                insert_data1_dict["all_num"] = all_num1
        data1 = pd.read_table(prediction_file, sep='\t')   #index_col=0
        data1.rename(columns={'#group': '#sample'}, inplace=True)
        data1 = data1.set_index('#sample')
        #data2 = pd.read_table(group_file, sep='\t')
        group = pd.DataFrame(pd.read_table(group_file, sep='\t', dtype={"#sample": "object"}))
        table_name = [col for col in data1.columns if col not in['ASV ID']]
        table = pd.DataFrame(data1[table_name])
        df = group.join(table.T, on="#sample",how="right")

        table_name1 = [col for col in df.columns]
        for i in table_name1[2:]:
            df[i] = df[i].astype("float")
        group_sample = df.groupby(group.columns[1]).mean().T
        for i in range(len(group_sample)):
            insert_data2_dict = {
                "faprotax_id": main_id,
                "functional_groups": group_sample.index[i],
                "type":"group"
            }
            all_num2 = 0
            for x in group_sample.columns:
                insert_data2_dict[x] = float(group_sample[x][i])
                all_num2 += float(group_sample[x][i])
            insert_data2_dict["all_num"] = all_num2
            insert_data2.append(insert_data2_dict)

        heatmap_data = {"heatmap_data":{"name":"functional_groups","data":sample_list}}
        heatmap_info = json.dumps(heatmap_data, sort_keys=True, separators=(',', ':'))
        #sample_list.insert(0,"functional_groups")
        table_dict = {"table_data": sample_list}
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        try:
            collection1 = self.db["faprotax_detail"]
            collection1.insert_many(insert_data1)
            collection2 = self.db["faprotax_detail"]
            collection2.insert_many(insert_data2)
            main = self.db['faprotax']
            main.update_one({'_id': main_id},{"$set": {"main_id": main_id,"specimen_list" : ",".join(sample_list),"table_data" : table_info, "heatmap_data" : heatmap_info}})
        except Exception as  e:
            self.bind_object.logger.info("导入faprotax详情数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入faprotax详情数据成功")

