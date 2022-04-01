# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import datetime
import json
import re
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.models.mongo.core.base import Base
from bson.son import SON
# from biocluster.config import Config


class EstimatorsMongo(Base):
    def __init__(self, bind_object):
        super(EstimatorsMongo, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self._params = self.PackParams()

    def PackParams(self):
        data = self.bind_object.data
        my_param = dict()
        my_param['otu_id'] = str(data.otu_id)
        my_param['level_id'] = int(data.level_id)
        # my_param['indices'] = data.index_type
        sort_index = data.index_type.split(',')
        sort_index.sort()
        sort_index = ','.join(sort_index)
        my_param['indices'] = sort_index
        my_param['submit_location'] = data.submit_location
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        return params

    def add_est_collection(self, name=None):
        data = self.bind_object.data
        if int(data.level_id) not in range(1, 10):
            raise Exception("level参数%s为不在允许范围内!" % data.level_id)
        if data.otu_id != 0 and not isinstance(data.otu_id, ObjectId):
            if isinstance(data.otu_id, StringTypes):
                data.otu_id = ObjectId(data.otu_id)
            else:
                raise Exception("传入的otu_id必须为ObjectId对象或其对应的字符串")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": data.otu_id})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "多样性指数表"
        insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": data.otu_id,
                "name": name if name else "estimators_origin",
                "level_id": data.level_id,
                "status": "end",
                "desc": desc,
                "params": self._params,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        collection = self.db["sg_alpha_diversity"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_est_detail(self, file_path, est_id):
        data_list = []
        with open(file_path, 'r') as f:
            # self.bind_object.logger.info(file_path)
            l = f.readline().strip('\n')
            if not re.match(r"^sample", l):
                raise Exception("文件%s格式不正确，请选择正确的estimator表格文件" % file_path)
            est_list = l.split("\t")
            # est_list.pop()
            est_list.pop(0)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                sample_name = line_data.pop(0)
                # line_data.pop()
                data = [("alpha_diversity_id", est_id), ("specimen_name", sample_name)]
                i = 0
                for est in est_list:
                    # self.bind_object.logger.info(str(line_data))
                    data.append((est, float(line_data[i])))
                    i += 1
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_alpha_diversity_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入estimator%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入estimator%s信息成功!" % file_path)
