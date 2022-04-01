# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson.objectid import ObjectId


class Specimen(Base):
    def __init__(self, bind_object):
        super(Specimen, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.sample_table_ids = list()
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])

    @report_check
    def add_samples_info(self):
        self.bind_object.logger.info("add_samples_info start")
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            'task_id': self.main_task_id,
            'name': 'Specimen_Origin',
            'status': 'end',
            'desc': '-',
            "submit_location": "specimen",
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        try:
            collection = self.db["specimen"]
            inserted_id = collection.insert_one(insert_data).inserted_id
            return inserted_id
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")

    def add_specimen_detail(self, update_id, file_path):
        """
        导入详情表
        :param update_id：主表specimen的id
        :param file_path: 要导入的样本路径
        :return:
        """
        self.bind_object.logger.info("update_id:{}\n{}".format(update_id, file_path))
        if not isinstance(update_id, ObjectId):
            main_id = ObjectId(update_id)
        else:
            main_id = update_id
        self.bind_object.logger.info("add_samples_info start")
        data_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "specimen": line_data[0],
                    "new_specimen": line_data[0],
                    "desc": '-',
                    "specimen_id":main_id,
                    'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                }
                data_list.append(data)
        try:
            collection = self.db["specimen_detail"]
            result = collection.insert_many(data_list)
            self.sample_table_ids = result.inserted_ids[:]
            self.sample_ids = self.get_spname_spid()
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % self.sample_ids)
        try:
            main_collection = self.db["specimen"]
            main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
        except Exception, e:
            self.bind_object.logger.error("更新主表信息数据出错:%s" % e)
            self.bind_object.set_error("更新主表信息数据出错")
        else:
            self.bind_object.logger.info("更新信息数据成功！")
        return self.sample_ids

    @report_check
    def get_spname_spid(self):
        if not self.sample_table_ids:
            self.bind_object.logger.error("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
            self.bind_object.set_error("没有找到样本信息")
        collection = self.db["specimen_detail"]
        spname_spid = dict()
        for id_ in self.sample_table_ids:
            result = collection.find_one({"_id": id_})
            if not result:
                self.bind_objecct.set_error("意外错误，无法找到样本id")
            spname_spid[result["specimen"]] = id_
        return spname_spid

    def _find_specimen_id(self, results):
        specimen_id = ""
        for result in results:
            if result["_id"] in self.sample_table_ids:
                specimen_id = result["_id"]
                break
        return specimen_id
        
