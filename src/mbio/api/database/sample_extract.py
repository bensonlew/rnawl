# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# last modified by guhaidong 20171109

import os
import json
import datetime
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config
from biocluster.api.database.base import Base, report_check


class SampleExtract(Base):
    def __init__(self, bind_object):
        super(SampleExtract, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def update_sg_seq_sample(self, list_path, table_id):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51007201")
        file_sample = list()
        length_sample = list()
        #workdir_sample = list()
        with open(list_path, "rb") as r:
            for line in r:
                if "#" in line:
                    continue
                line = line.rstrip().split("\t")
                name = os.path.basename(line[0])
                path = line[1]
                work_path = line[2]
                line.pop(0)
                line.pop(0)
                line.pop(0)
                length = "\t".join(line)
                sample = name + "/" + path
                file_sample.append({name: path})
         #       if {name:work_path} not in workdir_sample:
         #           workdir_sample.append({name:work_path})
                length_sample.append({sample:length})
        collection = self.db["sg_seq_sample"]
        results = collection.find_one({"_id": table_id})
        if not results:
            self.bind_object.logger.error("table_id:{}在sg_seq_sample表里未找到".format(table_id))
            self.bind_object.set_error("sg_seq_sample没有对应信息", code="51007202")
        results["sample_info"] = json.dumps(file_sample)
        results["name"] = "SampleExtract"
        results["length"] = json.dumps(length_sample)  # 新增字段length
        results["info_path"] = list_path  # add by zhujuan 2018.04.28 新增该字段，解决当“检测样品”后工作流中重复进行检测的bug
        #results["workdir_sample"] = json.dumps(workdir_sample)  # 新增字段，对应样本检测时的工作目录
        try:
            # collection.find_one_and_update({"_id": table_id}, {'$set': results})
            results["main_id"] =  table_id
            collection.update({"_id": table_id}, {'$set': results})
            self.bind_object.logger.info("表格导入成功")
        except Exception as e:
            self.bind_object.logger.error("表格导入出错:{}".format(e))
            self.bind_object.set_error("表格导入出错", code="51007203")


    @report_check
    def add_seq_sample_detail(self, list_path, table_id, real_name=None, s3_name=None):
        """
        插入序列的信息
        :param list_path:整理形成的info_path
        :param table_id: 主表ID
        :return:
        """
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        self.bind_object.logger.info("list_path: {}".format(list_path))
        self.bind_object.logger.info("table_id: {}".format(table_id))
        with open(list_path, "r") as r:
            r.readline()
            for line in r:
                line = line.strip().split("\t")
                if "#file" in line[0]:
                    continue
                else:
                    file_name = os.path.basename(line[0])
                    sample_name = line[1]
                    new_sample = line[1]
                    min = line[6]
                    mean = line[5]
                    max = line[7]
                    insert_data = {
                        "check_id": table_id,
                        "specimen": sample_name,
                        "is_check": "true",
                        "origin_specimen": new_sample,
                        "read_num": int(line[3]),
                        "base_num": int(line[4]),
                        "min": min,
                        "mean": mean,
                        "max": max
                    }
                    if real_name:
                        insert_data[ "file_name"] = real_name
                    else:
                        insert_data[ "file_name"] = file_name
                    if s3_name:
                        true_file_name = ""
                        if file_name in s3_name:
                            true_file_name = s3_name[file_name]
                        elif real_name in s3_name:
                            true_file_name = s3_name[real_name]
                        insert_data[ "true_file_name"] = true_file_name
                    data_list.append(insert_data)

        try:
            collection = self.db["sg_sample_check_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_sample_check"]
            main_collection.update({"_id": table_id}, {"$set": {"main_id": table_id, "status": "end"}})
            self.bind_object.logger.info("表格导入成功")
        except Exception as e:
            self.bind_object.logger.error("表格导入出错:{}".format(e))
            self.bind_object.set_error("表格导入出错")

    def add_sample_check(self, task_id, file_path, params, query_id):
        """
        多样性导入样本检测主表
        :param task_id: 任务id
        :param file_path: info_path
        :param params: params
        :param query_id: 文件或者文件夹的id
        :return:
        """
        insert_data = {
            "task_id": task_id,
            "file_path": file_path,
            "params": params,
            "query_id": query_id,
            'status': 'start',
            'desc': 'sample_check主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["sg_sample_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    #@report_check
    def check_sample(self, task_id, query_id):
        """
        检查sample_check主表是否存在
        :param task_id:
        :return:
        """
        collection = self.db["sg_sample_check"]
        result = collection.find_one({"task_id": task_id, "query_id": query_id, "status": "end"})
        if result:
            return True
        else:
            return False

    def get_file_info(self, task_id, query_id):
        """
        从样本检测主表中获取file_info信息
        :param task_id:
        :return:
        """
        collection = self.db["sg_sample_check"]
        result = collection.find_one({"task_id": task_id, "query_id": query_id})
        if result:
            if 'params' in result:
                file_list = json.loads(result['params'])['file_info']['file_list']
                return file_list
        else:
            return False
