# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os,re
import json
import datetime
from bson.objectid import ObjectId
from types import StringTypes
from biocluster.api.database.base import Base, report_check


class SampleCheck(Base):
    def __init__(self, bind_object):
        super(SampleCheck, self).__init__(bind_object)
        self._project_type = 'metaasv'

    @report_check
    def add_seq_sample(self, params,task_id, name=None, query_id=None):
        project_sn = self.bind_object.sheet.project_sn
        if not name:
            name = "SeqSample_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.bind_object.logger.info("params:{}".format(params))
        insert_data = {
            "project_sn": project_sn,
            'task_id': task_id,
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            'desc': 'sample_check主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        if query_id:
            insert_data["query_id"] = str(query_id)
        else:
            insert_data["query_id"] = ""
        collection = self.db["sample_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_seq_sample_detail(self, list_path, table_id,real_name=None,s3_name=None):
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
                    # if re.search(r".s.fq", file_name):
                    #     file_name = file_name.rstrip(".s.fq")
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
                        if file_name in s3_name :
                            true_file_name = s3_name[file_name]
                        elif real_name in s3_name:
                            true_file_name = s3_name[real_name]
                        insert_data[ "true_file_name"] = true_file_name
                    data_list.append(insert_data)

        try:
            collection = self.db["sample_check_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["sample_check"]
            main_collection.update({"_id": table_id}, {"$set": {"main_id": table_id}})
            self.bind_object.logger.info("表格导入成功")
        except Exception as e:
            self.bind_object.logger.error("表格导入出错:{}".format(e))
            self.bind_object.set_error("表格导入出错")

    @report_check
    def check_sample(self, task_id, query_id):
        """
        检查sample_check主表是否存在
        :param task_id:
        :return:
        """
        collection = self.db["sample_check"]
        result = collection.find_one({"task_id": task_id, "query_id": query_id, "status": "end"})
        if result:
            return True
        else:
            return False

    def delete_sg_status(self, task_id, query_id):
        """
        删除sg_status
        :param task_id:
        :param query_id:
        :return:
        """
        collection = self.db["sg_status"]
        result = collection.find_one({"task_id": task_id, "type_name" : "sample_check"})
        if result:
            main_id = result["_id"]
            try:
                collection.remove({"_id": main_id})
            except:
                self.bind_object.set_error("删除sg_status记录失败！")
