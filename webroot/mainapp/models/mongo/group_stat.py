# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

from bson.objectid import ObjectId
import datetime
from types import StringTypes
from mainapp.models.mongo.core.base import Base
# from mainapp.config.db import get_mongo_client
# from biocluster.config import Config


class GroupStat(Base):
    def __init__(self, bind_object=None):
        super(GroupStat, self).__init__(bind_object)
        self._project_type = 'meta'
        # self.client = get_mongo_client()
        # self.db = self.client[Config().MONGODB]

    def create_species_difference_check(self, level, check_type, params, category_name, group_id=0, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception("group_id必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "正在进行组间差异性检验"
        if check_type == 'tow_sample':
            insert_data = {
                "type": check_type,
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": from_otu_table,
                "name": name if name else "difference_stat_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
                "level_id": int(level),
                "params": params,
                "desc": desc,
                "status": "start",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        else:
            insert_data = {
                "type": check_type,
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": from_otu_table,
                "group_id": group_id,
                "name": name if name else "difference_stat_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
                "level_id": int(level),
                "params": params,
                "desc": desc,
                "status": "start",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "category_name": category_name
            }
        collection = self.db["sg_species_difference_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def create_species_difference_lefse(self, params, group_id=0, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception("group_detail必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "正在进行lefse分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name if name else "Lefse_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "lda_cladogram_id": "",
            "lda_png_id": "",
            "desc": desc,
            "status": "start",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_species_difference_lefse"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def get_group_name(self, group_id, lefse=False, second_group=''):
        """
        根据分组方案id获取分组方案名字
        :param group_id: 分组方案id
        :return: 分组方案名字
        """
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception("group_detail必须为ObjectId对象或其对应的字符串!")
        collection = self.db['sg_specimen_group']
        result = collection.find_one({'_id': group_id})
        gname = result['group_name']
        if lefse and second_group:
            gname = gname + ',' + 'second_group'
        return gname

    def get_otu_sample_name(self, otu_id):
        """
        获取otu表样本名字
        """
        if isinstance(otu_id, StringTypes):
            otu_id = ObjectId(otu_id)
        elif isinstance(otu_id, ObjectId):
            otu_id = otu_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        # collection_1 = self.db['sg_otu']
        # result = collection_1.find_one({"_id": otu_id})
        # _id = result["otu_id"]
        collection_2 = self.db['sg_otu_specimen']
        result_2 = collection_2.find({"otu_id": otu_id})
        collection_3 = self.db['sg_specimen']
        sample_name = []
        for i in result_2:
            specimen_id = ObjectId(i['specimen_id'])
            result_3 = collection_3.find_one({"_id": specimen_id})
            sample_name.append(result_3["specimen_name"])
        return sample_name
