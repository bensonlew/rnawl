# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20191104

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class PanCategory(Base):
    def __init__(self, bind_object):
        super(PanCategory, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        # self.id = 'tsg_123'
        # self.project_sn = '188_5b5acb3018'

    @report_check
    def add_category(self, pan_id, pan_category, params=None, name=None, category_scheme=None):
        if not isinstance(pan_id, ObjectId):
            new_inserted_id = ObjectId(pan_id)
        else:
            new_inserted_id = pan_id
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        # task_id = self.id
        # project_sn = self.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "分组方案合并主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pan_category_Origin",
            "category_scheme": category_scheme if category_scheme else "number1",
        }
        collection = self.db["pan_genome"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id':inserted_id, 'pan_id': new_inserted_id, "pan_category":pan_category}})
        return inserted_id

    def add_category_detail(self, inserted_id, category_path, pan_id=None, type=None):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        self.bind_object.logger.info("category_path:path %s" % (category_path))
        with open(category_path, "r") as f:
            lines = f.readlines()
            header = lines[0].lower()
            head_list = header.strip().split("\t")
            # categ_list = head_list[3:]
            for line in lines[1:]:
                line = line.strip().split("\t")
                sample_name = line[0]
                if re.search(r'.faa', sample_name):
                    sample_name.replace(".faa", "")
                total_list = line[1]
                total_num = float(line[2])
                data = {
                    "pan_genome_id": new_inserted_id,
                    "spe_name": sample_name,
                    "total_list": total_list,
                    "total_num": total_num
                    }
                if type:
                    data["type"] = type
                for i in range(3, len(head_list)):
                    class_name = head_list[i]
                    if i % 2 == 1:
                        data[class_name] = line[i]
                    else:
                        data[class_name] = float(line[i])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["pan_genome_specimen"]
            collection.insert_many(data_list)
            if pan_id:
                if not isinstance(pan_id, ObjectId):
                    new_pan_id = ObjectId(pan_id)
                else:
                    new_pan_id = pan_id
                main_collection = self.db["pan_genome"]
                main_collection.update_one({'_id': new_inserted_id}, {'$set': {'pan_id': new_pan_id,
                                                                               "main_id": new_inserted_id}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (category_path, e))

    def add_category_graph(self, inserted_id, distribution_path):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        self.bind_object.logger.info("distribution_path:path %s" % (distribution_path))
        with open(distribution_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                sample_num = float(line[0])
                cluster_num = float(line[1])
                data = {
                    "pan_genome_id": new_inserted_id,
                    "spe_num": sample_num,
                    "cluster_num": cluster_num,
                }
                data_son = SON(data)
                data_list.append(data_son)
            try:
                collection = self.db["pan_genome_bar"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (distribution_path, e))

    def add_pan_group(self, category, category_name, value_type, number, params=None, name=None):
        """
        工作流插入方案一和方案二使用
        :param category:
        :param category_name:
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "分组方案主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pan_group_Origin",
            "category_scheme": number,
        }
        collection = self.db["pan_group"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id':inserted_id, 'category_names': category_name, "category":category, "value_type": value_type}})
        return inserted_id

    def update_sg_status(self, table_id, data):
        self.update_table('sg_status', table_id, data, search_id='table_id')

    def update_table(self, table, table_id, data, search_id='_id'):
        tb = self.db[table]
        tb.update({search_id: ObjectId(table_id)},{'$set': data})