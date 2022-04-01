# -*- coding: utf-8 -*-


from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
import pandas as pd
import re
from types import StringTypes


class SampleVenn(Base):
    def __init__(self, bind_object):
        super(SampleVenn, self).__init__(bind_object)
        self._project_type = "metabolome"

    #@report_check
    def add_venn(self, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "完成",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "version" : "3.0"

        }
        collection = self.db["sample_venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': inserted_id}, {'$set': {'main_id': inserted_id}})
        return inserted_id

    #@report_check
    def add_venn_detail(self, venn_path, venn_id,group_num, metab_desc):
        data_list = []
        data = pd.read_table(metab_desc, sep='\t', index_col=0)
        sub_data = data[data['Metabolite'].map(lambda x : False if re.match("^metab_\d+$",x)  or x=='-' else True)]
        has_name_metab_set = set(sub_data.index.tolist())

        with open(venn_path, 'rb') as r:
            for line in r:
                line = line.strip().split("\t")
                label = line[0]
                if '&' not in label:
                    if 'only'in label:
                        label = label.split('only')[0].strip()
                else:
                    if group_num > 6 and len(label.split("&"))!=group_num :
                        continue
                    label = '&'.join((sorted([i.strip() for i in label.split("&")])))
                insert_data = {
                    'venn_id': ObjectId(venn_id),
                    'label': label,
                    'metab_num': int(line[1]),
                    }
                if int(line[1]):
                    insert_data['metab_list'] = line[2]
                    c_metabs = set(line[2].split(',')) & has_name_metab_set
                    insert_data['metab_num_f'] = len(c_metabs)
                    insert_data['metab_list_f'] = ','.join(sorted(c_metabs, key=lambda x: int(x.split("_")[1])))
                else:
                    insert_data['metab_list'] = ""
                    insert_data['metab_num_f'] = 0
                    insert_data['metab_list_f'] = ""

                data_list.append(insert_data)

        try:
            collection = self.db["sample_venn_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入venn数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入venn数据成功")

