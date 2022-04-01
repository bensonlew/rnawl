# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import json
import datetime
from bson.objectid import ObjectId


class DataStat(Base):
    def __init__(self, bind_object):
        super(DataStat, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.sample_table_ids = list()
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])

    @report_check
    def add_datastat(self):
        self.bind_object.logger.info("add_datastat start")
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            'task_id': self.main_task_id,
            'name': 'Datastat_Origin',
            'status': 'end',
            'desc': 'Datastat_Origin',
            "submit_location": "data_stat",
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        try:
            collection = self.db["data_stat"]
            inserted_id = collection.insert_one(insert_data).inserted_id
            return inserted_id
        except Exception, e:
            self.bind_object.logger.error("导入样本统计表失败:%s！" % e)
            self.bind_object.set_error("导入样本统计表失败！")

    @report_check
    def add_datastat_detail(self, update_id, file_path, type, column_number=None):
        """
        导入详情表
        :param update_id：主表specimen的id
        :param file_path: 要导入的样本路径
        :param type: type类型
        :param column_number 列数，只有5和6 列两种类型
        :return:
        """
        if not isinstance(update_id, ObjectId):
            main_id = ObjectId(update_id)
        else:
            main_id = update_id
        self.bind_object.logger.info("add_datastat_detail start")
        data_list = []
        total_reads = 0
        total_base = 0
        with open(file_path, 'r') as f:
            l = f.readline()
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                if column_number in [5]:
                    data = {
                        "specimen": "-",
                        "amplified": line_data[0],
                        "insert": int(line_data[1]),
                        "length": line_data[2],
                        "reads_num": str(line_data[3]),
                        "base_num": str(line_data[4]),
                        "stat_id": main_id,
                        "type": type
                    }
                else:
                    # if line_data[0] != "Total":
                    ##注释掉的原因是产品线要求raw的情况不做加和的计算 20200617 张俊彪
                    data = {
                        "specimen": line_data[0],
                        "amplified": line_data[1],
                        "insert": int(line_data[2]),
                        "length": line_data[3],
                        "reads_num": str(line_data[4]),
                        "base_num": str(line_data[5]),
                        "stat_id": main_id,
                        "type": type
                    }
                    total_reads += int(line_data[4])
                    total_base += int(line_data[5])
                data_list.append(data)
        ##注释掉的原因是产品线要求raw的情况不做加和的计算 20200617 张俊彪
        # if column_number not in [5]:
        #     insert_data = {
        #         "specimen": "Total",
        #         "amplified": "--",
        #         "insert": "--",
        #         "length":"--",
        #         "reads_num": total_reads,
        #         "base_num": total_base,
        #         "stat_id": main_id,
        #         "type": type
        #     }
        #     data_list.append(insert_data)
        try:
            collection = self.db["data_stat_detail"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)
        try:
            main_collection = self.db["data_stat"]
            raw_table_data = {
                "table_data": ["specimen", "amplified", "insert", "length", "reads_num", "base_num"],
                            "condition": {"type": "raw"}}
            table_data_json = json.dumps(raw_table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id,
                                                                "raw_table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功")

    @report_check
    def add_datastat_clean(self, update_id, file_path, type):
        """
        导入详情表
        :param update_id：主表specimen的id
        :param file_path: 要导入的样本路径
        :param type: type类型
        :return:
        """
        if not isinstance(update_id, ObjectId):
            main_id = ObjectId(update_id)
        else:
            main_id = update_id
        self.bind_object.logger.info("add_datastat_clean start")
        data_list = []
        total_reads = 0
        total_base = 0
        total_mean = 0
        sample_num = 0
        min_list = []
        max_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "specimen": line_data[0],
                    "reads_num": int(line_data[1]),
                    "base_num": int(line_data[2]),
                    "mean": float(line_data[3]),
                    "min": int(line_data[4]),
                    "max":int(line_data[5]),
                    "stat_id": main_id,
                    "type": type
                }
                sample_num += 1
                total_reads += int(line_data[1])
                total_base += int(line_data[2])
                total_mean += float(line_data[3])
                min_list.append(int(line_data[4]))
                max_list.append(int(line_data[5]))
                data_list.append(data)
        max_list.sort(reverse=True)
        min_list.sort(reverse=False)

        insert_data = {
            "specimen": "Total",
            "mean": round(total_mean / float(sample_num),6),
            "min": min_list[0],
            "max": max_list[0],
            "reads_num": int(total_reads),
            "base_num": int(total_base),
            "stat_id": main_id,
            "type": type
        }
        data_list.append(insert_data)
        try:
            collection = self.db["data_stat_detail"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)
        try:
            main_collection = self.db["data_stat"]
            clean_table_data = {
                "table_data": ["specimen", "reads_num", "base_num", "mean", "min", "max"],
                            "condition": {"type": "clean"}}
            table_data_json = json.dumps(clean_table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id,
                                                                "clean_table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功")

    @report_check
    def add_datastat_denoise(self, update_id, file_path, type, method):
        """
        导入详情表
        :param update_id：主表specimen的id
        :param file_path: 要导入的样本路径
        :param type: type类型
        :param method: method降噪方法
        :return:
        """
        if not isinstance(update_id, ObjectId):
            main_id = ObjectId(update_id)
        else:
            main_id = update_id
        self.bind_object.logger.info("add_datastat_clean start")
        data_list = []
        total_reads = 0
        total_base = 0
        with open(file_path, 'r') as f:
            l = f.readline()
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "specimen": line_data[0],
                    "reads_num": int(float(line_data[2])),
                    "asv_num": int(float(line_data[1])),
                    "stat_id": main_id,
                    "type": type
                }
                # total_reads += int(line_data[2])
                # total_base += int(line_data[1])
                data_list.append(data)
        try:
            collection = self.db["data_stat_detail"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)

        try:
            main_collection = self.db["data_stat"]
            if method in ['DADA2', 'dada2']:
                settled_params = {'software': "fastp (v0.19.6),FLASH (v1.2.7), DADA2",}
            else:
                settled_params = {'software': "fastp (v0.19.6), FLASH(v1.2.7), Deblur"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            denoise_table_data = {
                "table_data": ["specimen", "asv_num", "reads_num"],
                            "condition": {"type": "denoise"}}
            table_data_json = json.dumps(denoise_table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id,
                                                                "settled_params": settled_params_json,
                                                                "denoise_table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错")
        else:
            self.bind_object.logger.info("导入样品信息数据成功")