# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180522
from biocluster.api.database.base import Base, report_check
import os, json
import datetime
import types
from bson.son import SON
from bson.objectid import ObjectId


class Repeat(Base):
    def __init__(self, bind_object):
        super(Repeat, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def rm_main(self, collect_name, name=None):
        """
        根据task_id, project_sn, name中Origin名字删除orgin表
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        collection = self.db[collect_name]
        if name:
            result = collection.find({"task_id": task_id, "project_sn": project_sn, "name": name})
        else:
            result = collection.find({"task_id": task_id, "project_sn": project_sn, "name": {'$regex': '_Origin'}})
        if result:
            for each in result:
                mongo_id = each["_id"]
                try:
                    collection.remove({"_id": mongo_id})
                except Exception as e:
                    self.bind_object.set_error("删除主表%s信息出错:%s" , variables=(collect_name, e), code="54701501")
                else:
                    self.bind_object.logger.info("删除主表%s信息成功!" % collect_name)

    @report_check
    def rm_each(self, collect_name, name=None, name_detail=None):
        collection = self.db[collect_name]
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if name:
            if not name_detail:
                self.bind_object.set_error('if name exists, name_detial must exists!', code="54701502")
            result = collection.find({"task_id": task_id, "project_sn": project_sn, name: name_detail})
        else:
            result = collection.find({"task_id": task_id, "project_sn": project_sn})
        if result:
            for each in result:
                mongo_id = each["_id"]
                try:
                    collection.remove({"_id": mongo_id})
                except Exception as e:
                    self.bind_object.set_error("删除主表%s信息出错:%s" , variables=(collect_name, e), code="54701503")
                else:
                    self.bind_object.logger.info("删除主表%s信息成功!" % collect_name)