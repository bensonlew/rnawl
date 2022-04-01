# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.3.12

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re,os
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import shutil


class Cgview(Base):
    def __init__(self, bind_object):
        super(Cgview, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_cgview(self,file,genome_type=None,specimen_id=None,name = None,params=None, project_sn=None, task_id=None,main_id=None,main=False,update=True,link_path =None):
        task_id = task_id if task_id else self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        collection = self.db["cgview"]
        if main:
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "status": "faile",
                "name": name,
                "version":"3.0",
                "desc": 'cgview分析画图',
                "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "genome_type":genome_type,
                "specimen_id":specimen_id,
            }
            self.bind_object.logger.info(insert_data)
            main_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        else:
            if not main_id:
                #raise Exception('不写入主表时，需要提供主表ID')
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="51402101")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        collection.update_one({"_id":main_id},{'$set': {'img_path': link_path + '','main_id':main_id}})
        if update:
            collection = self.db["cgview"]
            xml_path = link_path + 'xml/cgview.xml'
            self.bind_object.logger.info(main_id)
            collection.update({"_id": ObjectId(main_id)},
                                   {"$set": {"status": 'end', "xml_path": xml_path,'main_id': main_id}})
        else:
            collection = self.db["cgview"]
            self.bind_object.logger.info(main_id)
            collection.update({"_id": ObjectId(main_id)},
                              {"$set": {"status": 'end'}})


