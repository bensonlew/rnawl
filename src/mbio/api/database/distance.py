# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
from biocluster.api.database.base import Base, report_check
import re
import json
import datetime
from bson.son import SON
# from types import StringTypes
from bson.objectid import ObjectId
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort


class Distance(Base):
    def __init__(self, bind_object):
        super(Distance, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_dist_table(self, file_path, major=False, level=None, dist_id=None, otu_id=None, task_id=None, name=None, params=None, spname_spid=None):
        if level and level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在允许范围内！", code="51002401")
        # if task_id is None:
        #     task_id = self.bind_object.sheet.id
        data_list = []
        if otu_id:
            if not isinstance(otu_id, ObjectId):
                otu_id = ObjectId(otu_id)
            params['otu_id'] = str(otu_id)  # otu_id在再metabase中不可用
        if spname_spid:
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        # insert major
        if major:
            if not otu_id:
                self.bind_object.set_error('写主表时必须提供otu_id', code="51002402")
            collection = self.db["sg_otu"]
            result = collection.find_one({"_id": otu_id})
            task_id = result['task_id']
            insert_data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": task_id,
                "otu_id": otu_id,
                "level_id": level,
                "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "distance_origin",
                "status": "end",
                "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
            collection = self.db["sg_beta_specimen_distance"]
            dist_id = collection.insert_one(insert_data).inserted_id
            #collection.update({"_id": dist_id}, {"$set": {"main_id": dist_id}})
        else:
            if dist_id is None:
                self.bind_object.set_error("major为False时需提供dist_id!", code="51002403")
            else:
                if not isinstance(dist_id, ObjectId):
                    dist_id = ObjectId(dist_id)
        # insert detail
        with open(file_path, 'r') as f:
            l = f.readline().strip('\n')
            if not re.match(r"^\t", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的distance表格文件" % file_path)
                self.bind_object.set_error("distance表格文件格式错误", code="51002404")
            sample_list = l.split("\t")
            sample_list.pop(0)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                sample_name = line_data.pop(0)
                data = [("specimen_distance_id", dist_id), ("specimen_name", sample_name)]
                i = 0
                for smp in sample_list:
                    data.append((smp, float(line_data[i])))
                    i += 1
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_beta_specimen_distance_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_beta_specimen_distance"]
            #main_collection.update({"_id": dist_id}, {"$set": {"main_id": dist_id}})
        except Exception as e:
            self.bind_object.logger.error("导入distance%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入distance%s信息成功!" % file_path)
        return dist_id
