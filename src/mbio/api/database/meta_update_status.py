# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.api.database.base import Base, report_check
# import re
from bson.objectid import ObjectId
import datetime
import json
from bson.son import SON
from types import DictType, StringTypes
# from biocluster.config import Config


class MetaUpdateStatus(Base):
    def __init__(self, bind_object):
        super(MetaUpdateStatus, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_meta_status(self, table_id=None, table_name=None, type_name=None, status='end', task_id=None, isnew='new', desc=None, submit_location=None, params=None):
        collection = self.db["sg_status"]
        if not task_id:
            task_id = self.bind_object.sheet.id
        else:
            if params is not None:
                if not isinstance(params, dict):
                    self.bind_object.set_error('提供的参数params必须为字典', code="51004001")
        if isinstance(table_id, StringTypes):
            try:
                table_id = ObjectId(table_id)
            except Exception as e:
                self.bind_object.logger.error('参数table_id不是mongo表ID类型:{}'.format(e))
                self.bind_object.set_error("table_id不是ObjectId类型", code="51004002")
        elif isinstance(table_id, ObjectId):
            pass
        else:
            self.bind_object.set_error('提供的table_id不是字符串或者ObjectId类型', code="51004003")
        if not type_name:
            self.bind_object.set_error('必须提供type_name!!!', code="51004004")
        if not table_name or not params:
            temp_collection = self.db[type_name]
            tempfind = temp_collection.find_one({'_id': table_id})
            if not tempfind:
                self.bind_object.logger.error('提供的ID:%s无法在表:%s中找到' % (table_id, type_name))
                self.bind_object.set_error("找不到相应记录", code="51004005")
            if 'name' in tempfind:
                table_name = tempfind['name']
            else:
                self.bind_object.logger.error('表：%s中没有name字段' % type_name)
                self.bind_object.set_error("表中没有name字段", code="51004006")
            if 'params' in tempfind:
                params = tempfind['params']
        if not submit_location:
            if isinstance(params, StringTypes):
                if 'submit_location' in params:
                    temp_params = json.loads(params)
                    submit_location = temp_params['submit_location']
            elif isinstance(params, DictType):
                submit_location = params['submit_location'] if 'submit_location' in params else None
        insert_data = {
            "table_id": table_id,
            "table_name": table_name,
            "task_id": task_id,
            "type_name": type_name,
            "status": status,
            "is_new": isnew,
            "desc": desc,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, DictType) else params,
            "submit_location": submit_location,
            "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        insert_data = SON(insert_data)
        return_id = collection.insert_one(insert_data).inserted_id
        if return_id:
            return return_id
        else:
            self.bind_object.set_error('插入sg_status表出错', code="51004007")

    # by houshuang 20190925 更新用户选择的数据库至sg_task表
    @report_check
    def add_database(self, task_id, database):
        try:
            collection = self.db["sg_task"]
            tempfind = collection.find_one({'task_id': task_id})
            collection.update_one({'_id': tempfind["_id"]}, {'$set': {'database': database}})
        except Exception as e:
            self.bind_object.logger.error("database导入sg_task表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("database导入sg_task表格成功")

    @report_check
    def add_sample_numbers(self, task_id, sample_numbers):
        try:
            collection = self.db["sg_task"]
            tempfind = collection.find_one({'task_id': task_id})
            collection.update_one({'_id': tempfind["_id"]}, {'$set': {'sample_numbers': sample_numbers}})
        except Exception as e:
            self.bind_object.logger.error("sample_numbers导入sg_task表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("sample_numbers导入sg_task表格成功")

    @report_check
    def add_save_pdf(self, task_id, save_pdf):
        try:
            collection = self.db["sg_task"]
            tempfind = collection.find_one({'task_id': task_id})
            collection.update_one({'_id': tempfind["_id"]}, {'$set': {'save_pdf': save_pdf}})
        except Exception as e:
            self.bind_object.logger.error("sample_numbers导入sg_task表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("sample_numbers导入sg_task表格成功")