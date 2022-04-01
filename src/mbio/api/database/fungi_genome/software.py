# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Software(Base):
    def __init__(self, bind_object):
        super(Software, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_software(self,task_id=None, project_sn=None,params=None, name=None):
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn if project_sn else self.bind_object.sheet.project_sn,
            'task_id': task_id,
            'desc': '软件列表',
            'created_ts': created_ts,
            'name': name if name else 'Soft_Origin',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
        }
        try:
            collection = self.db["software"]
            inserted_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % ("software", e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %"software")
        return inserted_id

    @report_check
    def add_software_detail(self,main_id):
        """
        将软件改成与任务挂钩，便于以后更新软件版本
        :param main_id:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        software_old = self.db['software_config']
        collection = self.db['software_detail']
        data_list = []
        results = software_old.find({})
        if results:
            for result in results:
                new_result = result
                new_result.pop("_id")
                new_result["soft_id"] = main_id
                data = SON(new_result)
                data_list.append(data)
            try:
                collection.insert_many(data_list)
                main_collection = self.db['software']
                main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
            except Exception as e:
                self.bind_object.logger.error("导入software_detail%s信息出错:%s" % (main_id, e))
                self.bind_object.set_error("导入software_detail%s信息出错" )
            else:
                self.bind_object.logger.info("导入software_detail%s信息成功!" % main_id)

