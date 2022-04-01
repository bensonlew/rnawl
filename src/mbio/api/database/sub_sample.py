# -*- coding: utf-8 -*-
# __author__ = 'xuting'

from biocluster.api.database.base import Base, report_check
import re
import datetime
import json
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config


class SubSample(Base):
    def __init__(self, bind_object):
        super(SubSample, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.name_id = dict()  # otu表中样本名和id对照的字典
        self.otu_rep = dict()  # o

    def add_sg_otu(self, params, my_size, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007501")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51007502")
        project_sn = result['project_sn']
        task_id = result['task_id']
        if not name:
            name = "otu_subsample" + str(my_size) + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            'task_id': task_id,
            'from_id': str(from_otu_table),
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            'status': 'end',
            "level_id": json.dumps([9]),
            'desc': 'otu table after Otu Subsampe',
            "type": "otu_statistic",
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_otu"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def _get_name_id(self, from_otu_id):
        collection = self.db['sg_otu_specimen']
        results = collection.find({"otu_id": from_otu_id})
        if not results.count():
            self.bind_object.logger.error("otu_id:{}未在otu_sg_specimen表里找到相应的记录".format(from_otu_id))
            self.bind_object.set_error("otu_sg_specimen表找不到相应记录", code="51007503")
        sp_ids = list()
        for result in results:
            sp_ids.append(result['specimen_id'])
        collection = self.db['sg_specimen']
        for id_ in sp_ids:
            result = collection.find_one({"_id": id_})
            if not result:
                self.bind_object.logger.error("意外错误， id: {}在sg_otu_specimen表中找到，但未在sg_specimen表中出现")
                self.bind_object.set_error("sg_specimen表找不到相应记录", code="51007504")
            self.name_id[result["specimen_name"]] = id_

    def _get_task_info(self, otu_id):
        collection = self.db['sg_otu']
        result = collection.find_one({'_id': otu_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(otu_id))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51007502")
        self.project_sn = result['project_sn']
        self.task_id = result['task_id']

    def _prepare_otu_rep(self, from_otu_id):
        """
        实际运行的时候发现对每一行(即每一个otu)都去数据库里查询一次，并获取otu_rep的时候，效率非常低，需要很长的时间， 因此，需要对mongo做一次查询， 将属于一个otu_id的otu_rep全部去读出来， 放到内存当中， 以提高效率
        """
        self.bind_object.logger.info("开始依据otu_id: {}查询所有的代表序列".format(from_otu_id))
        collection = self.db["sg_otu_detail"]
        results = collection.find({"otu_id": from_otu_id})
        for result in results:
            self.otu_rep[result['otu']] = result["otu_rep"]
        self.bind_object.logger.info("代表序列查询完毕")

    @report_check
    def add_sg_otu_detail(self, file_path, from_otu_id, new_otu_id, new_samples=False):
        if not isinstance(from_otu_id, ObjectId):
            if isinstance(from_otu_id, StringTypes):
                from_otu_id = ObjectId(from_otu_id)
            else:
                self.bind_object.set_error("from_otu_id必须为ObjectId对象或其对应的字符串!", code="51007501")
        if not isinstance(new_otu_id, ObjectId):
            if isinstance(new_otu_id, StringTypes):
                new_otu_id = ObjectId(new_otu_id)
            else:
                self.bind_object.set_error("new_otu_id必须为ObjectId对象或其对应的字符串!", code="51007505")
        self._get_name_id(from_otu_id)
        self._get_task_info(new_otu_id)
        # 导入sg_otu_detail表
        self.bind_object.logger.info("开始导入sg_otu_detail表")
        self.bind_object.logger.info("file_path:%s"%file_path)
        self._prepare_otu_rep(from_otu_id)
        insert_data = list()
        with open(file_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                otu_detail = dict()
                # OTU表可能有taxonomy列， 也可能没有， 需要适应
                if len(re.split(";", line[0])) > 1:
                    sample_num = line[1:]
                    otu = re.split(";", line[0])[-1]
                    classify_list = re.split(r"\s*;\s*", line[0])
                else:
                    sample_num = line[1:-1]
                    otu = line[0]
                    otu_detail['otu'] = line[0]
                    classify_list = re.split(r"\s*;\s*", line[-1])

                otu_detail['otu_id'] = new_otu_id
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(sample_num[i])

                if otu not in self.otu_rep:
                    self.bind_object.logger.error("意外错误，otu_id: {}和otu: {}在sg_otu_detail表里未找到".format(from_otu_id, line[0]))
                    self.bind_object.set_error("sg_otu_detail表中找不到相应记录", code="51007506")
                otu_detail['otu_rep'] = self.otu_rep[otu]
                otu_detail['task_id'] = self.task_id
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_otu_detail']
            collection.insert_many(insert_data)

            main_collection = self.db["sg_otu"]
            #main_collection.update({"_id": new_otu_id},{"$set": { "main_id": new_otu_id}})

        except Exception as e:
            self.bind_object.logger.error("导入sg_otu_detail表格信息出错:{}".format(e))
            self.bind_object.set_error("导入sg_otu_detail表格信息出错", code="51007507")
        else:
            self.bind_object.logger.info("导入sg_otu_detail表格成功")
        # 导入sg_otu_specimen表
        self.bind_object.logger.info("开始导入sg_otu_specimen表")
        if not new_samples:
            insert_data = list()
            for sp in new_head:
                my_data = dict()
                my_data['otu_id'] = new_otu_id
                my_data["specimen_id"] = self.name_id[sp]
                insert_data.append(my_data)
            collection = self.db['sg_otu_specimen']
            collection.insert_many(insert_data)

    @report_check
    def add_sg_otu_detail_level(self, otu_path, from_otu_table, level):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007501")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51007502")
        project_sn = result['project_sn']
        task_id = result['task_id']
        covered_level = list()
        if "level_id" in result:
            covered_level = json.loads(result["level_id"])
            covered_level.append(int(level))
        else:
            covered_level.append(int(level))
        covered_level = list(set(covered_level))
        covered_level.sort()
        result["level_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                otu_detail = dict()

                if len(re.split("; ", line[0])) > 1:
                    sample_num = line[1:]
                    classify_list = re.split(r"\s*;\s*", line[0])
                elif len(re.split(";", line[0])) > 1:
                    sample_num = line[1:]
                    classify_list = re.split(r"\s*;\s*", line[0])
                else:
                    sample_num = line[1:-1]
                    otu_detail['otu'] = line[0]
                    classify_list = re.split(r"\s*;\s*", line[-1])

                otu_detail['otu_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                count = 0
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(sample_num[i])
                    count += int(sample_num[i])
                otu_detail["total_"] = count
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_otu_detail_level']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入sg_otu_detail_level表格失败：{}".format(e))
            self.bind_object.set_error("导入sg_otu_detail_level表格失败", code="51007508")
        else:
            self.bind_object.logger.info("导入sg_otu_detail_copy表格成功")

    ###guanqing.zou 20180505
    
    @report_check
    def add_sg_otu_seq_summary(self,otu_path, from_otu_table):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007501")

        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            if head[-1] == "taxonomy":
                new_head = head[1:-1]
            else:
                new_head = head[1:]
            sample_num = len(new_head)
            sample_count ={}
            for index in range(sample_num):    #guanqing.zou 20180525
                sample_count[new_head[index]] = 0
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                for index in range(sample_num):    #guanqing.zou 20180525  
                    sample_count[new_head[index]] += int(line[index+1])

        otu_sum = 0
        insert_data = []
        for k in sample_count.keys():
            insert_tmp = {"otu_id":from_otu_table}
            insert_tmp["specimen_name"] = k
            insert_tmp["read_number"] = sample_count[k]
            insert_data.append(insert_tmp)
            otu_sum += sample_count[k]
        try:
            collection = self.db['sg_otu_seq']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入sg_otu_seq表格失败：{}".format(e))
            self.bind_object.set_error("导入sg_otu_seq表格失败", code="51007509")
        else:
            self.bind_object.logger.info("导入sg_otu_seq表格成功")

        find_task = self.db['sg_otu'].find_one({"_id":from_otu_table})
        task_id = find_task["task_id"]
        sg = self.db['sg_valid_sequence_info'].find_one({"task_id":task_id})
        # add excution of empty result situation; by liulinmeng 20180611
        if not sg:
            amplified_region = "--"
        else:
            amplified_region = sg['amplified_region']
        try:
            collection = self.db['sg_otu_summary']
            collection.insert_one({"otu_id":from_otu_table,"samples":sample_num,"sequences":otu_sum,"amplified_region":amplified_region})
        except Exception as e:
            self.bind_object.logger.error("导入sg_otu_summary表格失败：{}".format(e))
            self.bind_object.set_error("导入sg_otu_summary表格失败", code="51007510")
        else:
            self.bind_object.logger.info("导入sg_otu_summary表格成功")

    @report_check
    def add_otu_detail(self, otu_path, from_otu_table, level):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            raise Exception("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        covered_level = list()
        if "level_id" in result:
            covered_level = json.loads(result["level_id"])
            covered_level.append(int(level))
        else:
            covered_level.append(int(level))
        covered_level = list(set(covered_level))
        covered_level.sort()
        result["level_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['otu_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                count = 0
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                    count += int(float(sample_num[i]))
                otu_detail["total_"] = count
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_otu_detail_level']
            collection.insert_many(insert_data)
        except Exception as e:
            raise Exception("导入sg_otu_detail_level表格失败：{}".format(e))
        else:
            print "导入sg_otu_detail_level表格成功"