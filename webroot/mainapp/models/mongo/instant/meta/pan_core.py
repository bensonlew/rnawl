# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import datetime
import re
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.models.mongo.core.base import Base
# from biocluster.config import Config


class PanCoreMongo(Base):
    def __init__(self, bind_object):
        super(PanCoreMongo, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self._params = self.PackParams()

    def PackParams(self):
        data = self.bind_object.data
        params = dict()
        params['otu_id'] = str(data.otu_id)
        params['level_id'] = str(data.level_id)
        params['group_id'] = str(data.group_id)
        c_name = re.split(',', data.category_name)
        c_name.sort()
        new_cname = ','.join(c_name)
        params['category_name'] = new_cname
        params["submit_location"] = data.submit_location
        params["task_type"] = data.task_type
        params = self.SortDict(params)
        return params

    def CreatePanCoreTable(self, pan_core_type, name=None):
        data = self.bind_object.data
        if data.otu_id != 0 and not isinstance(data.otu_id, ObjectId):
            if isinstance(data.otu_id, StringTypes):
                data.otu_id = ObjectId(data.otu_id)
            else:
                raise Exception("传入的otu_id必须为ObjectId对象或其对应的字符串")
        if data.group_id != "all":
            if not isinstance(data.group_id, ObjectId):
                if isinstance(data.group_id, StringTypes):
                    data.group_id = ObjectId(data.group_id)
                else:
                    raise Exception("传入的group_id必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": data.otu_id})
        if not result:
            raise Exception("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(data.otu_id)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        if pan_core_type == 1:
            desc = "pan otu表格"
        else:
            desc = "core otu表格"
        insertData = {
            "type": pan_core_type,
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": data.level_id,
            "otu_id": data.otu_id,
            "group_id": data.group_id,
            "status": "end",
            "desc": desc,
            "name": name if name else "pan_core表格",
            "params": self._params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_otu_pan_core"]
        print "正在往Mongo的sg_pan_core表里插入一条记录，名称为{}".format(name)
        self.bind_object.logger.info("正在往Mongo的sg_pan_core表里插入一条记录，名称为{}".format(name))
        inserted_id = collection.insert_one(insertData).inserted_id
        return inserted_id

    def AddPanCoreDetail(self, filePath, panCoreId):
        if not isinstance(panCoreId, ObjectId):
            if isinstance(panCoreId, StringTypes):
                panCoreId = ObjectId(panCoreId)
            else:
                raise Exception("panCoreId必须为ObjectId对象或其对应的字符串!")
        with open(filePath, 'rb') as r:
            header = r.next().rstrip("\n")
            header = re.split('\t', header)
            panCoreDetail = list()
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                value = list()
                length = len(line)
                for i in range(1, length):
                    myDict = dict()
                    myDict["title"] = header[i]
                    myDict["value"] = line[i]
                    value.append(myDict)
                insertData = {
                    "pan_core_id": panCoreId,
                    "category_name": line[0],
                    "value": value
                }
                panCoreDetail.append(insertData)
            try:
                collection = self.db['sg_otu_pan_core_detail']
                collection.insert_many(panCoreDetail)
            except Exception as e:
                self.bind_object.logger.error("导入pan_core_detail表格{}信息出错:{}".format(filePath, e))
            else:
                print "导入pan_core_detail表格{}成功".format(filePath)
                self.bind_object.logger.info("导入pan_core_detail表格{}成功".format(filePath))
