# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20160920
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
import json
from bson.son import SON
from types import DictType, StringTypes
from biocluster.config import Config


class DenovoUpdateStatus(Base):
    def __init__(self, bind_object):
        super(DenovoUpdateStatus, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'

    @report_check
    def add_denovo_status(self, table_id=None, table_name=None, type_name=None, status='end', task_id=None, isnew='new', desc=None, submit_location=None, params=None):
        collection = self.db["sg_status"]
        if not task_id:
            task_id = self.bind_object.sheet.id
        else:
            if params is not None:
                if not isinstance(params, dict):
                    raise Exception('提供的参数params必须为字典')
        if isinstance(table_id, StringTypes):
            try:
                table_id = ObjectId(table_id)
            except Exception as e:
                raise Exception('参数table_id不是mongo表ID类型:{}'.format(e))
        elif isinstance(table_id, ObjectId):
            pass
        else:
            raise Exception('提供的table_id不是字符串或者ObjectId类型')
        if not type_name:
            raise Exception('必须提供type_name!!!')
        if not table_name or not params:
            temp_collection = self.db[type_name]
            tempfind = temp_collection.find_one({'_id': table_id})
            if not tempfind:
                raise Exception('提供的ID:%s无法在表:%s中找到' % (table_id, type_name))
            if 'name' in tempfind:
                table_name = tempfind['name']
            else:
                raise Exception('表：%s中没有name字段' % type_name)
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
            raise Exception('插入sg_status表出错')
