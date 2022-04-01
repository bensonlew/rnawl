# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2018.04.16
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId


class Permanova(Base):
    def __init__(self, bind_object):
        super(Permanova, self).__init__(bind_object)
        self._project_type = "meta"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_permanova(self, file_path, main=False, task_id=None, main_id=None,
                      anno_type=None, params=None, name=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'PERMANOVA_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'anno_type': anno_type
            }
            collection = self.db['sg_permanova']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                self.bind_object.set_error("main为False时需提供main_id!", code="51004901")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        data_list = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            collection_detail = self.db["sg_permanova_detail"]
            for line in lines[1:]:
                data = [("permanova_id", main_id)]
                line = line.strip().split('\t')
                data.extend([("cha", line[0]), ("sumsqs", line[1]), ("meanseq", line[2]), ("f", line[3]),
                             ("r2", line[4]), ("p", line[5]), ("padjust", line[6])])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            main_collection = self.db["sg_permanova"]
            #main_collection.update({"_id":main_id},{"$set":{"main_id": main_id}})

        except Exception as e:
            self.bind_object.logger.error("导入permanova%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入permanova信息出错", code="51004902")
        else:
            self.bind_object.logger.info("导入permanova%s信息成功!" % file_path)
        return main_id
