# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# last_modify:20181015
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
import types
import unittest
from mbio.api.database.small_rna.api_base import ApiBase

class SmallRnaQc(ApiBase):
    def __init__(self, bind_object):
        super(SmallRnaQc, self).__init__(bind_object)
        self._project_type = 'small_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    def add_before_qc(self, qc_stat, fq_type='se', about_qc='before', group=None):
        """

        """
        self.add_samples_info(qc_stat, fq_type=fq_type, about_qc=about_qc, group=group)
        quality_stat_before = qc_stat + "/qualityStat"
        self.add_gragh_info(quality_stat_before, "before")


    def add_after_qc(self, qc_stat, qc_dir=None,  fq_type='se', about_qc='after', group=None):
        self.add_samples_info(qc_stat, fq_type=fq_type, qc_dir=qc_dir, about_qc=about_qc, group=group)
        quality_stat_after = qc_stat + "/qualityStat"
        self.add_gragh_info(quality_stat_after, "after")
        qclen_id = self.add_qclen()
        self.add_qclen_detail(qclen_id, qc_stat, qc_dir)
        self.update_db_record('sg_qclen', qclen_id, status="end", main_id=qclen_id)

    @report_check
    def add_samples_info(self, qc_stat, qc_adapt=None, qc_dir=None, fq_type='se', about_qc='before', group=None):
        """
        :param qc_stat: 统计结果文件夹，即module.output_dir
        :param qc_adapt:去接头率文件，由于需求变动，可不传
        :param fq_type:测序类型
        :return:
        """
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        stat_file = qc_stat + "/fastq_stat.xls"
        '''
        dup_file = qc_stat + "/dup.xls"
        dup = ""
        dup_rate = {}
        '''
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
        '''
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
        '''
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
                if data["about_qc"] == "after":
                    with open(qc_dir + "/" + line[0] + "_qc_stat.xls", 'r') as f:
                        f.readline()
                        values = f.readline().split("\t")
                        data.update({
                            "useful_reads": values[6]
                        })

                num = num + 1
                if adapter:
                    if line[0] in adapt_rate:
                        data["adapt_rate"] = adapt_rate[line[0]]
                        data["about_qc"] = "after"
                data_list.append(data)

            if group:
                data_list.sort(key=lambda x: sample_list.index(x['old_name']))
            try:
                collection = self.db["sg_specimen"]
                result = collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入样品信息数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入样品信息数据成功")
                self.sample_ids = result.inserted_ids

        for i in result.inserted_ids:
            self.update_db_record('sg_specimen', i, status="end", main_id=i)

        sample_ids = [str(i) for i in result.inserted_ids]
        return sorted(sample_ids)

    @report_check
    def add_qclen(self, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'QcLengthStat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': "",
            # 'result_dir': result_dir,
            'status': 'end',
            'desc': '质控长度统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_qclen']
        qclen_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_qclen")
        return qclen_id

    @report_check
    def add_qclen_detail(self, qclen_id, qc_stat, qc_dir):
        if not isinstance(qclen_id, ObjectId):
            if isinstance(qclen_id, types.StringTypes):
                qclen_id = ObjectId(qclen_id)
            else:
                raise Exception('qclen_id必须为ObjectId对象或其对应的字符串！')

        stat_file = qc_stat + "/fastq_stat.xls"

        adapter = False
        if not os.path.exists(stat_file):
            raise Exception("%s文件不存在" % stat_file)

        with open(stat_file, "r") as f:
            data_list = []
            first_line = f.readline()
            num = 1
            if not re.match(r"#Sample_ID", first_line):
                raise Exception("%s文件类型不正确" % stat_file)
            for line in f:
                line = line.split()
                with open(qc_dir + "/" + line[0] + "_clean.length.txt", 'r') as f1:
                    #values = f.readline().split("\t")
                    data = {
                        "qclen_id": qclen_id,
                        "name":line[0],
                    }
                    len_stat = [{line.strip().split()[0]: int(line.strip().split()[1])} for line in f1.readlines()]
                    data.update({
                        "len_data": len_stat,
                    })
                data_list.append(data)

            try:
                collection = self.db["sg_qclen_detail"]
                result = collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入样品信息数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入样品信息数据成功")


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
            spname_spid = self.get_spname_spid()
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
                    total_base = int(line[12]) + int(line[13]) + int(line[14]) + int(line[15]) + int(line[16])
                    data = {
                        "project_sn": self.bind_object.sheet.project_sn,
                        "task_id": self.bind_object.sheet.id,
                        "specimen_name": sample_name,
                        "specimen_id": spname_spid[sample_name],
                        "type": site_type,
                        "about_qc": about_qc,
                        "column": int(line[0]),
                        "min": int(line[10]),
                        "max": int(line[11]),
                        "q1": int(line[6]),
                        "q3": int(line[8]),
                        "median": int(line[7]),
                        "error": 10 ** (float(line[5])/(-10)) * 100,
                        "A": int(line[12])/total_base * 100,
                        "T": int(line[15])/total_base * 100,
                        "C": int(line[13])/total_base * 100,
                        "G": int(line[14])/total_base * 100,
                        "N": int(line[16])/total_base * 100
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
    def get_spname_spid(self):
        if not self.sample_ids:
            raise Exception("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
        collection = self.db["sg_specimen"]
        spname_spid = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_spid[results['new_name']] = id_
        return spname_spid

    @report_check
    def add_specimen_group(self, file, productive_table=None):
        #spname_alias = self.get_spname_alias()
        category_names = list()
        specimen_names = list()
        group_dict = OrderedDict()
        if productive_table:
            productive_names = dict()
            mj_names = dict()
            with open(productive_table, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    items = line.strip().split("\t")
                    if len(items) >= 2:
                        productive_names[items[0]] = items[1]
                    if len(items) >= 3:
                        mj_names[items[0]] = items[2]
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
                if productive_table and sample in productive_names:
                    if productive_table and sample in mj_names:
                        col.update({"task_id": self.bind_object.sheet.id, "old_name": sample, "about_qc": "before"},
                                   {"$set": {"group": key, 'productive_name': productive_names[sample], 'mj_name': mj_names[sample]}})
                        col.update({"task_id": self.bind_object.sheet.id, "old_name": sample, "about_qc": "after"},
                                   {"$set": {"group": key, 'productive_name': productive_names[sample], 'mj_name': mj_names[sample]}})
                    else:
                        col.update({"task_id": self.bind_object.sheet.id, "old_name": sample, "about_qc": "before"},
                                   {"$set": {"group": key, 'productive_name': productive_names[sample]}})
                        col.update({"task_id": self.bind_object.sheet.id, "old_name": sample, "about_qc": "after"},
                                   {"$set": {"group": key, 'productive_name': productive_names[sample]}})
                else:
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
    def test_mongo(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "small_rna",
            "project_sn": "small_rna",
            "type": "workflow",
            "name": "small_rna.small_rna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.api_qc = wf.api.api("small_rna.small_rna_qc")

        test_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data_hsa/"

        after_qc_stat = test_dir + "before_qc"
        before_qc_stat = test_dir + "before_qc"
        quality_stat_after = after_qc_stat + "/qualityStat"
        quality_stat_before = before_qc_stat + "/qualityStat"  # 将qc前导表加于该处
        qc_dir = test_dir + "qc"
        fq_type = "SE"
        wf.api_qc.add_before_qc(before_qc_stat, fq_type=fq_type, about_qc="before")
        wf.api_qc.add_after_qc(after_qc_stat, qc_dir=qc_dir, fq_type=fq_type, about_qc="after")

        '''
        wf.group_id, wf.group_detail, wf.group_category = wf.api_qc.add_specimen_group(wf.option("group_table").prop["path"])
        wf.logger.info("group_detail为：" + str(wf.group_detail))
        wf.control_id, compare_detail = wf.api_qc.add_control_group(wf.option("control_file").prop["path"], wf.group_id)
        wf.compare_detail = compare_detail
        wf.api_qc.add_bam_path(wf.workflow_output)
        '''


if __name__ == '__main__':
    unittest.main()
