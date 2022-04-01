# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# last_modify:20171212
from __future__ import division
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from bson.objectid import ObjectId
# from cStringIO import StringIO
# import bson.binary
import datetime
# import pandas
# import numpy
from bson.son import SON
from collections import OrderedDict
import json
import glob
import re
import os
import unittest
from mbio.api.database.small_rna.api_base import ApiBase

class SmallRnaTrim(ApiBase):
    def __init__(self, bind_object):
        super(SmallRnaTrim, self).__init__(bind_object)
        self._project_type = 'small_rna'

    @report_check
    def add_samples_info(self, qc_stat, qc_adapt=None, fq_type='se', about_qc='after'):
        """
        :param qc_stat: 统计结果文件夹，即module.output_dir
        :param qc_adapt:去接头率文件，由于需求变动，可不传
        :param fq_type:测序类型
        :return:
        """
        stat_file = qc_stat + "/fastq_stat.xls"
        dup_file = qc_stat + "/dup.xls"
        dup = ""
        dup_rate = {}
        adapter = False
        if not os.path.exists(stat_file):
            raise Exception("%s文件不存在" % stat_file)
        if qc_adapt is not None:
            adapt_rate = {}
            adapter = True
            with open(qc_adapt, "r") as a:
                for line in a:
                    line = line.split()
                    adapt_rate[line[0]] = line[1]
        if os.path.exists(dup_file):
            with open(dup_file, "r") as d:
                col_num = len(d.readline().split())
                if col_num == 4:
                    dup = "pe"
                    for line in d:
                        line = line.split()
                        dup_rate[line[0]] = [float(line[1]), float(line[2]), float(line[3])]
                if col_num == 2:
                    dup = "se"
                    for line in d:
                        line = line.split()
                        dup_rate[line[0]] = line[1]
        with open(stat_file, "r") as f:
            data_list = []
            first_line = f.readline()
            num = 1
            if not re.match(r"#Sample_ID", first_line):
                raise Exception("%s文件类型不正确" % stat_file)
            for line in f:
                line = line.split()
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "created_ts" : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "old_name": line[0],
                    "new_name": line[0],
                    "desc" : "",
                    "group" : "",
                    "alias" : "S" + str(num),
                    "total_reads": int(line[1]),
                    "total_bases": int(line[2]),
                    "reads_with_ns": int(line[3]),
                    "n_reads_rate": float(line[4]),
                    "a_rate": float(line[5]),
                    "t_rate": float(line[6]),
                    "c_rate": float(line[7]),
                    "g_rate": float(line[8]),
                    "n_rate": float(line[9]),
                    "error_rate": float(line[10]),
                    "q20_rate": float(line[11]),
                    "q30_rate": float(line[12]),
                    "gc_rate": float(line[13]),
                    "about_qc": about_qc,
                    "type": fq_type
                    }
                num = num + 1
                if line[0] in dup_rate:
                    if dup == "se":
                        data["reads1_dup_rate"] = dup_rate[line[0]]
                    if dup == "pe":
                        data["reads1_dup_rate"] = dup_rate[line[0]][0]
                        data["reads2_dup_rate"] = dup_rate[line[0]][1]
                        data["paired_dup_rate"] = dup_rate[line[0]][1]
                if adapter:
                    if line[0] in adapt_rate:
                        data["adapt_rate"] = adapt_rate[line[0]]
                        data["about_qc"] = "after"
                data_list.append(data)

            try:
                collection = self.db["sg_specimen"]
                result = collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入样品信息数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入样品信息数据成功")
                self.sample_ids = result.inserted_ids
        sample_ids = [str(i) for i in result.inserted_ids]
        return sorted(sample_ids)

    @report_check
    def add_gragh_info(self, quality_stat, about_qc="before"):
        """
        :param quality_stat: ouput_dir里一个叫qualityStat的文件夹，即~/output_dir/qualityStat
        :param about_qc:质控后的或是质控前的统计，质控前统计为before，指控后传after
        :return:
        """
        stat_files = glob.glob("{}/*".format(quality_stat))
        data_list = []
        for sf in sorted(stat_files):
            sample_name = os.path.basename(sf).split(".")[0]
            self.bind_object.logger.info('%s,%s' % (sf, sample_name))
            #spname_alias = self.get_spname_alias()
            #sp_alias = spname_alias[sample_name]
            site = os.path.basename(sf).split(".")[1]
            if site == "l": site_type = "left"
            elif site == "r": site_type = "right"
            else: site_type = "single"
            with open(sf, "r") as f:
                f.readline()
                for line in f:
                    line = line.strip().split()
                    data = {
                        "project_sn": self.bind_object.sheet.project_sn,
                        "task_id": self.bind_object.sheet.id,
                        "specimen_name": sample_name,
                        "type": site_type,
                        "about_qc": about_qc,
                        "column": int(line[0]),
                        "min": int(line[10]),
                        "max": int(line[11]),
                        "q1": int(line[6]),
                        "q3": int(line[8]),
                        "median": int(line[7]),
                        "error": 10 ** (float(line[5])/(-10)) * 100,
                        "A": int(line[12])/int(line[17]) * 100,
                        "T": int(line[15])/int(line[17]) * 100,
                        "C": int(line[13])/int(line[17]) * 100,
                        "G": int(line[14])/int(line[17]) * 100,
                        "N": int(line[16])/int(line[17]) * 100
                    }
                    data_list.append(data)
        try:
            collection = self.db["sg_specimen_graphic"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入样品画图数据信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品画图数据信息成功")

    @report_check
    def add_specimen_group(self, file):
        #spname_alias = self.get_spname_alias()
        category_names = list()
        specimen_names = list()
        group_dict = OrderedDict()
        with open(file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                group_dict.setdefault(tmp[1], list())
                if tmp[0] not in group_dict[tmp[1]]:
                    group_dict[tmp[1]].append(tmp[0])
        col = self.db["sg_specimen"]
        for key in group_dict:
            category_names.append(key)
            specimen_names.append(group_dict[key])
            for sample in group_dict[key]:
                col.update({"task_id" : self.bind_object.sheet.id, "old_name": sample, "about_qc":"before"}, {"$set": {"group": key}})
                col.update({"task_id" : self.bind_object.sheet.id, "old_name": sample, "about_qc":"after"}, {"$set": {"group": key}})

        data = {
            "task_id" : self.bind_object.sheet.id,
            "category_names": category_names,
            "specimen_names": specimen_names,
            "group_name": os.path.basename(file),
            "project_sn": self.bind_object.sheet.project_sn,
            "is_use" : 1,
        }
        col = self.db["sg_specimen_group"]
        group_id = col.insert_one(data).inserted_id
        col.update({"_id": group_id, "task_id" : self.bind_object.sheet.id}, {"$set": {"main_id": group_id}}, upsert=True)
        self.bind_object.logger.info("导入样本分组信息成功")
        return group_id, specimen_names, category_names

    @report_check
    def add_control_group(self,file, group_id):
        con_list = list()
        with open(file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                string = tmp[0] + "|" + tmp[1]
                con_list.append(string)
        col = self.db["sg_specimen_group"]
        result = col.find_one({"_id": group_id})
        group_name = result["group_name"]
        category_names = str(result["category_names"])
        data = {
            "task_id": self.bind_object.sheet.id,
            "compare_group_name": group_name,
            "compare_names": json.dumps(con_list),
            "compare_category_name": "all",
            "specimen_group_id": str(group_id),
            "is_use" : 1,
        }
        col = self.db["sg_specimen_group_compare"]
        try:
            com_id = col.insert_one(SON(data)).inserted_id
            col.update({"_id": com_id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": com_id}}, upsert=True)
        except Exception,e:
            self.bind_object.set_error("导入样本对照组信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样本对照组信息成功")
            return com_id, con_list

    def get_spname_id(self):##获取name-id键值对
        if not self.sample_ids:
            raise Exception("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
        collection = self.db["sg_specimen"]
        spname_id = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_id[results['new_name']] = id_
        return spname_id

    def get_spname_alias(self):##获取name-alias键值对
        if not self.sample_ids:
            raise Exception("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
        collection = self.db["sg_specimen"]
        spname_alias = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_alias[results['new_name']] = results['alias']
        return spname_alias

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    from mbio.workflows.denovo_rna_v2.denovo_test_api import SmallTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random

    data = {
        # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
        "id": "test1",
        "project_sn": "test1",
        "type": "workflow",
        "name": "denovo_rna_v2.denovo_test_api",
        "options": {
        },
    }
    wsheet = Sheet(data=data)
    wf = SmallTestApiWorkflow(wsheet)
    wf.IMPORT_REPORT_DATA = True
    wf.IMPORT_REPORT_AFTER_END = False
    wf.test_api = wf.api.api("small_rna.denovo_rna_qc")
    wf.test_api.run1()

if __name__ == '__main__':
    unittest.main()
