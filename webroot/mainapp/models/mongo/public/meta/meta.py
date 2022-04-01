# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
import re
import json
from mainapp.models.mongo.core.base import Base
# from biocluster.config import Config


class Meta(Base):
    def __init__(self, bind_object=None):
        super(Meta, self).__init__(bind_object)
        self._project_type = "meta"
        #self.client = get_mongo_client()
        #self.db = self.client[Config().MONGODB]

    def get_otu_table_info(self, otu_id):

        if isinstance(otu_id, types.StringTypes):
            otu_id = ObjectId(otu_id)
        elif isinstance(otu_id, ObjectId):
            otu_id = otu_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_otu']
        result = collection.find_one({"_id": otu_id})
        return result

    def get_task_info(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def sampleIdToName(self, sampleIds):
        """
        将一个用逗号隔开的样本ID的集合转换成样本名，返回一个用逗号隔开的样本名的集合
        """
        myIds = re.split("\s*,\s*", sampleIds)
        collection = self.db["sg_specimen"]
        mySampleNames = list()
        for id_ in myIds:
            if id_ == "":
                raise Exception("存在空的sample_id")
            if not isinstance(id_, ObjectId):
                if isinstance(id_, types.StringTypes):
                    id_ = ObjectId(id_)
                else:
                    raise Exception("样本id必须为ObjectId对象或者其对应的字符串！")
            result = collection.find_one({"_id": id_})
            if not result:
                raise Exception("无法根据传入的_id:{}在sg_speciem表里找到相应的记录".format(str(id_)))
            mySampleNames.append(result["specimen_name"])
        mySamples = ",".join(mySampleNames)
        return mySamples

    def group_detail_to_table(self, group_detail, group_path):
        """
        输入group_detail，返回一个group表
        """
        print group_detail
        with open(group_path, "wb") as f:
            f.write("#sample\t" + "group_name" + "\n")
        if not isinstance(group_detail, dict):
            try:
                table_dict = json.loads(group_detail)
            except Exception:
                raise Exception("生成group表失败，传入的group_datail不是一个字典或者是字典对应的字符串")
        if not isinstance(table_dict, dict):
            raise Exception("生成group表失败，传入的group_datail不是一个字典或者是字典对应的字符串")
        sample_table = self.db["sg_specimen"]
        with open(group_path, "ab") as f:
            for k in table_dict:
                for sp_id in table_dict[k]:
                    sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                    if not sp:
                        raise Exception("group_detal中的样本_id:{}在样本表sg_specimen中未找到".format(sp_id))
                    else:
                        sp_name = sp["specimen_name"]
                    f.write("{}\t{}\n".format(sp_name, k))
        return group_path
