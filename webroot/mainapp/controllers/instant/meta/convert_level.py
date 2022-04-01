# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# last_modify by qiuping 20170113

import web
import json
import datetime
import os
from mainapp.models.mongo.export_file import ExportFile
from biocluster.config import Config
from mainapp.models.mongo.core.base import Base
from bson.objectid import ObjectId
from types import StringTypes
import re
from mainapp.libs.signature import check_sig
from mainapp.libs.input_check import meta_check


class ConvertLevelAction(Base):
    def __init__(self, bind_object=None):
        super(ConvertLevelAction, self).__init__(bind_object)
        self._project_type = "meta"
        # self._client = Config().mongo_client
        # self._db_name = Config().MONGODB
        # self.db = self._client[self._db_name]

    @check_sig
    @meta_check
    def POST(self):
        data = web.input()
        postArgs = ['level_id', 'submit_location', "otu_id", "task_type"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)
        otu_path = os.path.join(Config().WORK_DIR, "tmp", datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + ".otu.xls")
        ExportFile().export_otu_table_by_level(data.otu_id, otu_path, data.level_id)
        self.add_otu_detail(otu_path, data.otu_id, data.level_id)
        info = dict()
        info["success"] = True
        info["info"] = "已成功完成计算"
        return json.dumps(info)

    def add_otu_detail(self, otu_path, from_otu_table, level):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            raise Exception("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        covered_level = list()
        if "level_id" in result:
            covered_level = json.loads(result["level_id"])
            covered_level.append(int(level))
        else:
            covered_level.append(int(level))
        covered_level = list(set(covered_level))
        covered_level.sort()
        result["level_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['otu_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                count = 0
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                    count += int(sample_num[i])
                otu_detail["total_"] = count
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_otu_detail_level']
            collection.insert_many(insert_data)
        except Exception as e:
            raise Exception("导入sg_otu_detail_level表格失败：{}".format(e))
        else:
            print "导入sg_otu_detail_level表格成功"
