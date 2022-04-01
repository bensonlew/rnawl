# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd


class Funguild(Base):
    def __init__(self, bind_object):
        super(Funguild, self).__init__(bind_object)
        self._project_type = "meta"

    @report_check
    def add_funguild(self,level_id,otu_id, name=None, params=None):
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
        collection = self.db["sg_funguild"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_detail(self, file_path,table_id):
        data_list = []
        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        data = pd.read_table(file_path, sep='\t')
        data.fillna('-',inplace=True)
        taxon_list = ['taxon_' + i for i in 'd,k,p,c,o,f,g,s'.split(',')]
        head_list = data.columns.tolist()
        sample_list = head_list[1:head_list.index('taxonomy')]
        for i in range(len(data)):
            insert_data = {
                "guild_id": table_id,
                "otu_name": data['OTU ID'][i],
                "guild": data['Guild'][i],
                "confidence" : data['Confidence Ranking'][i],
                "growth" : data['Growth Morphology'][i],
                "trait" : data['Trait'][i] ,
                "note" : data['Notes'][i],
                "source" : data['Citation/Source'][i],
                'trophic':data['Trophic Mode'][i]
            }
            for sample in sample_list:
                insert_data[sample] = data[sample][i]

            split_taxon = data['taxonomy'][i].split(';')
            for id,taxon in enumerate(taxon_list,0):
                try:
                    insert_data[taxon]= split_taxon[id].strip()
                except Exception as e:
                    insert_data[taxon]= '-'

            data_list.append(insert_data)
        try:
            collection = self.db["sg_funguild_detail"]
            collection.insert_many(data_list)

            main = self.db['sg_funguild']
            #main.update_one({'_id': table_id},{"$set": {"main_id": table_id}})

        except Exception as  e:
            self.bind_object.logger.info("导入Funguild 详情数据出错:%s" % e)

        else:
            self.bind_object.logger.info("导入Funguild 详情数据成功")

    @report_check
    def add_sum(self,file_path, table_id):

        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        data = pd.read_table(file_path,sep='\t',header=0)
        data.fillna('-',inplace=True)
        samples = data.columns[1:]
        insert_data = []
        for index in range(len(data)):
            tmp_data = {
                "guild_id": table_id,
                "guild" : data['Guild'][index]
            }

            for sample in samples:
                tmp_data[sample] = data[sample][index]

            insert_data.append(tmp_data)

        try:
            collection = self.db["sg_funguild_stat"]
            collection.insert_many(insert_data)
            main=self.db['sg_funguild']
            main.update_one({'_id':table_id},{"$set":{'specimen_list': ','.join(samples)}})

        except Exception as e:
            self.bind_object.logger.error("导入Funguild统计数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入Funguild统计数据成功")

    def add_path(self, main_id, path):
        """
        添加结果路径到主表
        :param task_id:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            table_id = ObjectId(main_id)
        collection = self.db["sg_funguild"]
        collection.update_one({'_id': ObjectId(main_id)}, {'$set': {'path': path}})