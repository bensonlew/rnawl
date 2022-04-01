# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last modified guhaidong 20171116
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import gridfs
import datetime
import os
# from biocluster.config import Config


class PipeCheckAll(Base):
    def __init__(self, bind_object):
        super(PipeCheckAll, self).__init__(bind_object)
        self._project_type = 'meta'
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add main_task_id by guhaidong 20171116
        # self._db_name = Config().MONGODB

    #@report_check
    def check_all(self, all_results, main_table_id):
        self.bind_object.logger.info('测试check_all函数')
        collection_pipe = self.db["sg_pipe_batch"]
        anaysis_num = []
        for id in all_results:
            self.bind_object.logger.info(id['id'])
            collection = self.db["sg_network"]
            result = collection.find_one({"_id": ObjectId(id['id']), "task_id": self.main_task_id})  # add main_task_id by guhaidong 201711116
            self.bind_object.logger.info(result['status'])
            if result['status'] != 'start':
                anaysis_num.append(result['status'])
        if len(all_results) == len(anaysis_num):
            self.bind_object.logger.info(len(anaysis_num))
            self.bind_object.logger.info("++++++++++")
            self.bind_object.logger.info(len(all_results))
            data = {
                "status": "end",
                "sub_analysis_id": all_results
            }
            try:
                collection_pipe.update({"_id": ObjectId(main_table_id)}, {'$set': data}, upsert=False)
            except Exception as e:
                self.bind_object.logger.error("任务信息导入sg_pipe_batch失败：%s" % e)
                self.bind_object.set_error("任务信息导入sg_pipe_batch失败", code="51005101")
            self.bind_object.logger.info('任务信息导入sg_pipe_batch成功。')
            m = True
        else:
            pass

        # insert_data = {
        #     "sub_analysis_id": all_results
        # }
        # collection = self.db["sg_pipe_batch"]
        # collection.insert_one(insert_data)
        # self.bind_object.logger.info('任务信息导入sg_pipe_batch成功。')
        return m

