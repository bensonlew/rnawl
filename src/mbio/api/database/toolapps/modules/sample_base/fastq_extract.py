# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import os
import re
from bson import SON
import datetime

class FastqExtract(Base):
    def __init__(self, bind_object):
        super(FastqExtract, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.check()

    def check(self):
        pass

    @report_check
    def run(self):
        """
        运行函数
        """
        self.export2database()
        for step in (20, 50, 100, 200):
            reads_len_info_path = self.work_dir + "/FastqExtract/LengthStat/output/" \
                                                  "reads_len_info/step_{}.reads_len_info.txt".format(str(step))
            self.bind_object.logger.info(reads_len_info_path)
            if not os.path.isfile(reads_len_info_path):
                raise Exception("找不到报告文件")
            self.add_reads_len_info(step, reads_len_info_path)

    def export2database(self):
        results_list = []
        with open(self.work_dir + "/FastqExtract/info.txt") as r:
            main_id = self.db['fastq_extract'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='fastq_extract',
                desc='样本拆分表',
                status='failed',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )).inserted_id
            r.readline()
            for line in r:
                result = {}
                tmp = line.strip().split("\t")
                result["main_id"] = main_id
                result["name"] = os.path.basename(tmp[0])
                result["sample"] = tmp[1]
                result["seq_num"] = tmp[3]
                result["base_num"] = tmp[4]
                result["average"] = tmp[5]
                result["min"] = tmp[6]
                result["max"] = tmp[7]
                results_list.append(result)
        self.db['fastq_info'].insert_many(results_list)
        self.db['fastq_extract'].update_one({'_id': main_id}, {'$set': {'status': 'end', }})

    def add_reads_len_info(self, step_length, file_path):
        data_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            if not re.match(r"^sample\t", l):
                raise Exception("文件%s格式不正确，请选择正确的碱基统计文件" % file_path)
            else:
                length_list = l.strip("\n").split("\t")
                length_list.pop(0)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                step_data = {}
                i = 0
                for step in length_list:
                    i += 1
                    step_data[step] = int(line_data[i])
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "sample": line_data[0],
                    "step": step_length,
                    "value": step_data
                }
                data_list.append(data)
        collection = self.db["fastq_length_stat"]
        collection.insert_many(data_list)
        self.bind_object.logger.info("导入步%s的步长序列长度统计成功" % step_length)
        self.bind_object.logger.info("导表成功")
