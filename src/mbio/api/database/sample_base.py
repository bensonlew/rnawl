# -*- coding: utf-8 -*-
# __author__ = 'shijin'

import os
import json
from bson.objectid import ObjectId
from types import StringTypes
from biocluster.config import Config
from biocluster.api.database.base import Base, report_check


class SampleBase(Base):
    def __init__(self, bind_object):
        super(SampleBase, self).__init__(bind_object)
        self._db_name = "samplebase"

    @report_check
    def add_sg_test_specimen_rna(self, sample, stat_path, file_sample):
        collection = self.db["sg_test_specimen"]
        results = {}
        with open(stat_path, "r") as file:
            for line in file:
                if line.startswith(sample):
                    tmp = line.strip().split("\t")
                    results["specimen_name"] = tmp[0]
                    results["total_reads"] = tmp[1]
                    results["total_bases"] = tmp[2]
                    results["total_reads_with_ns"] = tmp[3]
                    results["n_reads%"] = tmp[4]
                    results["A%"] = tmp[6]
                    results["T%"] = tmp[7]
                    results["C%"] = tmp[8]
                    results["G%"] = tmp[9]
                    results["N%"] = tmp[10]
                    results["Error%"] = tmp[11]
                    results["Q20%"] = tmp[12]
                    results["Q30"] = tmp[13]
                    results["GC"] = tmp[14]
        results["file_name"] = []
        for key in file_sample.keys():
            if file_sample[key] == sample:
                results["file_name"].append(key)
        results["file_name"] = ",".join(results["file_name"])  # 考虑到双端序列，进行这样的处理
        sample_id = collection.insert_one(results).insert_id
        self.bind_object.logger.info("表格导入成功")
        return sample_id

    @report_check
    def add_sg_test_batch_specimen(self, table_id, sample_id, sample):
        collection = self.db["sg_test_batch_specimen"]
        results_list = []
        results = {}
        results["test_batch_id"] = ObjectId(table_id)
        results["test_specimen_id"] = ObjectId(sample_id)
        results["alias_name"] = sample
        results_list.append(results)
        try:
            collection.insert_many(results_list)
            self.bind_object.logger.info("表格导入成功")
        except Exception as e:
            self.bind_object.logger.error("表格导入出错:{}".format(e))
            self.bind_object.set_error("表格导入出错", code="51007101")

    @report_check
    def add_sg_test_specimen_meta(self, sample, info_txt, file_sample):
        collection = self.db["sg_test_specimen"]
        results = {}
        with open(info_txt, "r") as file:
            for line in file:
                tmp = line.strip().split("\t")
                if not tmp[1] == sample:
                    continue
                # file\tsample\tworkdir\tseqs_num\tbase_num\tmean_length\tmin_length\tmax_length
                results["specimen_name"] = tmp[1]
                results["sequence_num"] = tmp[3]
                results["base_num"] = tmp[4]
                results["mean_length"] = tmp[5]
                results["min_length"] = tmp[6]
                results["max_length"] = tmp[7]
        for key in file_sample.keys():
            if file_sample[key] == sample:
                results["file_name"] = key
        sample_id = collection.insert_one(results).insert_id
        self.bind_object.logger.info("表格导入成功")
        return sample_id