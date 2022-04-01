# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
import datetime
import json
from bson.son import SON
from types import StringTypes
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort


class Estimator(Base):
    def __init__(self, bind_object):
        super(Estimator, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_est_table(self, file_path, level, major=False, otu_id=None, est_id=None, task_id=None, name=None,
                      params=None, spname_spid=None,index_type=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level不在允许范围内", code="51003101")
        # if task_id is None:
        #     task_id = self.bind_object.sheet.id
        data_list = []
        if not otu_id:
            self.bind_oject.set_error("需要otu_id", code="51003102")
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": otu_id})
        if task_id is None:
            task_id = result['task_id']
        # if otu_id:
        #     if not isinstance(otu_id, ObjectId):
        #         otu_id = ObjectId(otu_id)
        # params['otu_id'] = str(otu_id)  # otu_id在再metabase中不可用
        # if spname_spid:
        #     group_detail = {'All': [str(i) for i in spname_spid.values()]}
        #     params['group_detail'] = group_detail_sort(group_detail)
        # insert major
        if major:
            params['otu_id'] = str(otu_id)  # otu_id在再metabase中不可用

            if spname_spid:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
                params['group_detail'] = group_detail_sort(group_detail)
            insert_data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": task_id,
                "otu_id": otu_id,
                "index_type": index_type,
                "name": name if name else "Estimators_Origin",
                "level_id": level,
                "status": "end",
                "desc": "",
                "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            # if params is not None:
            #     insert_data['params'] = params
            collection = self.db["sg_alpha_diversity"]
            est_id = collection.insert_one(insert_data).inserted_id
            #collection.update({"_id": est_id}, {"$set": {"main_id": est_id}})
        else:
            if est_id is None:
                self.bind_object.set_error("major为False时需提供est_id!", code="51003103")
            if not isinstance(est_id, ObjectId):
                if isinstance(est_id, StringTypes):
                    est_id = ObjectId(est_id)
                else:
                    self.bind_object.set_error("est_id必须为ObjectId对象或其对应的字符串!", code="51003104")
        # insert detail
        with open(file_path, 'r') as f:
            # self.bind_object.logger.info(file_path)
            l = f.readline().strip('\n')
            if not re.match(r"^sample", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的estimator表格文件" % file_path)
                self.bind_object.set_error("estimator表格格式不正确", code="51003105")
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
            main_collection = self.db["sg_alpha_diversity"]
            #main_collection.update({"_id": est_id}, {"$set": {"main_id": est_id}})
        except Exception, e:
            self.bind_object.logger.error("导入estimator%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入estimator%s信息成功!" % file_path)
        return est_id
