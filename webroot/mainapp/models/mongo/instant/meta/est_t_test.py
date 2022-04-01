# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last modified by guhaidong 20171116

import datetime
import json
import types
from bson.son import SON
import re
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.models.mongo.core.base import Base
from mainapp.libs.param_pack import group_detail_sort
# from biocluster.config import Config


class EstTTestMongo(Base):
    def __init__(self, bind_object):
        super(EstTTestMongo, self).__init__(bind_object)
        self._project_type = 'meta'
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add main_task_id by guhaidong 20171116
        # self._db_name = Config().MONGODB
        self._params = self.PackParams()

    def PackParams(self):
        data = self.bind_object.data
        my_param = dict()
        my_param['alpha_diversity_id'] = data.est_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        my_param['submit_location'] = data.submit_location
        est_params = self.get_est_params(data.est_id)
        my_param['otu_id'] = str(est_params[0])
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        return params

    def get_est_table_info(self, est_id):
        if isinstance(est_id, types.StringTypes):
            est_id = ObjectId(est_id)
        elif isinstance(est_id, ObjectId):
            est_id = est_id
        else:
            raise Exception("输入est_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_alpha_diversity']
        result = collection.find_one({"_id": est_id, "task_id": self.main_task_id})  # add task_id by guhaidong 20171116
        return result

    def get_est_params(self, est_id):
        est_info = self.get_est_table_info(est_id)
        otu_id = est_info['otu_id']
        return [otu_id]

    def add_est_t_test_collection(self, name=None):
        data = self.bind_object.data
        if data.est_id != 0 and not isinstance(data.est_id, ObjectId):
            if isinstance(data.est_id, StringTypes):
                data.est_id = ObjectId(data.est_id)
            else:
                raise Exception("传入的est_id必须为ObjectId对象或其对应的字符串")
        collection = self.db["sg_alpha_diversity"]
        result = collection.find_one({"_id": data.est_id, "task_id": self.main_task_id})  # add task_id by guhaidong 20171116
        project_sn = result['project_sn']
        task_id = result['task_id']
        otu_id = result['otu_id']
        level_id = result['level_id']
        desc = "多样性指数T检验表"
        insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "otu_id": otu_id,
                "alpha_diversity_id": data.est_id,
                "name": name if name else "多样性指数T检验结果表",
                "level_id": int(level_id),
                "group_id": data.group_id,
                "status": "end",
                "desc": desc,
                "params": self._params,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        collection = self.db["sg_alpha_ttest"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_est_t_test_detail(self, file_path, table_id):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip('\n')
            group_list = re.findall(r'mean\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                length = len(line_data)
                i = 1
                for name in group_list:
                    data = [("alpha_ttest_id", table_id), ("index_type", line_data[0]), ("qvalue", line_data[length-1]), ("pvalue", line_data[length-2])]
                    data.append(("category_name", name))
                    data.append(("compare_name", self.get_another_name(name, group_list)))
                    data.append(("mean", str('%0.5g' % float(line_data[i]))))
                    data.append(("sd", str('%0.5g' % float(line_data[i+1]))))
                    i += 2
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_alpha_ttest_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)

    def get_another_name(self, name, group_list):
        another = ''
        for n in group_list:
            if n == name:
                pass
            else:
                another = n
        return another
