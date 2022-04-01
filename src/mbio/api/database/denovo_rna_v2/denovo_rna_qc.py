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
from mbio.api.database.denovo_rna_v2.api_base import ApiBase
import json
import glob
import re
import os
import unittest


class DenovoRnaQc(ApiBase):
    def __init__(self, bind_object):
        super(DenovoRnaQc, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'
        #self._db_name = Config().MONGODB + '_ref_rna'

    @report_check
    def add_samples_info(self, qc_stat, qc_adapt=None, fq_type='se', about_qc='before', group=None):
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
        dup_file = qc_stat + "/dup.xls"
        rfam_file=qc_stat+"/stat_results"
        rfam_rata={}
        if about_qc=="after":
            with open(rfam_file,"r") as rf:
                first_info_line=rf.readline()
                for info_line in rf.readlines():
                    info_line=info_line.split()
                    rfam_rata[info_line[0]]=float(info_line[-1])
        dup = ""
        dup_rate = {}
        adapter = False
        if not os.path.exists(stat_file):
            self.bind_object.set_error("%s文件不存在" , variables=( stat_file), code="52000701")
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
            if not re.match(r"#Sample_ID", first_line):
                self.bind_object.set_error("%s文件类型不正确" , variables=( stat_file), code="52000702")
            for line in f:
                line = line.split()
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "old_name": line[0],
                    "new_name": line[0],
                    "total_reads": int(line[1]),
                    "total_bases": int(line[2]),
                    "reads_with_ns": int(line[3]),
                    "n_reads_rate": float(line[4]),
                    "a_rate": float(line[5]),
                    "t_rate": float(line[6]),
                    "c_rate": float(line[7]),
                    "g_rate": float(line[8]),
                    "n_rate": float(line[9]),
                    "version":"v2",
                    "error_rate": float(line[10]),
                    "q20_rate": float(line[11]),
                    "q30_rate": float(line[12]),
                    "gc_rate": float(line[13]),
                    "about_qc": about_qc,
                    "desc": "",
                    "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    "type": fq_type, # 怎么得知待定
                    }
                if line[0] in rfam_rata:
                    data["rrna_ratio"] = rfam_rata[line[0]]
                if line[0] in dup_rate:
                    if dup == "se":
                        data["reads1_dup_rate"] = dup_rate[line[0]]
                    if dup == "pe":
                        data["reads1_dup_rate"] = dup_rate[line[0]][0]
                        data["reads2_dup_rate"] = dup_rate[line[0]][1]
                        data["paired_dup_rate"] = dup_rate[line[0]][2]
                if adapter:
                    if line[0] in adapt_rate:
                        data["adapt_rate"] = adapt_rate[line[0]]
                        data["about_qc"] = "after"
                data_list.append(data)
            if group:
                data_list.sort(key=lambda x: sample_list.index(x['old_name']))
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data_list)
            self.bind_object.logger.info("导入样品信息数据成功")
            self.sample_ids = result.inserted_ids
            for i in result.inserted_ids:
                self.update_db_record('sg_specimen', i, status="end", main_id=i)
        sample_ids = [str(i) for i in result.inserted_ids]
        return sorted(sample_ids)

    '''
    @report_check
    def add_samples_alias(self, alias_table, about_qc="before"):
        """
        :param alias_table: 别名表
        :param about_qc: 质控后的或是质控前的统计，质控前统计为before，指控后传after
        :return:
        """
        spname_id = self.get_spname_id()
        with open (alias_table,"r") as f:
            f.readline()
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    tmp = line.strip().split("\t")
                    sample_name = tmp[0]
                    alias_name = tmp[1]
                    sp_id = spname_id[sample_name]
                    col = self.db["sg_specimen"]
                    col.update({"task_id" : self.bind_object.sheet.id, "_id":sp_id, "about_qc":about_qc}, {"$set": {"alias": alias_name}})
        self.bind_object.logger.info("导入样本别名信息成功")
    '''

    @report_check
    def add_productive_name(self, samples=None, productive_table=None):
        collection = self.db["sg_specimen"]
        with open(productive_table, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                items = line.strip().split("\t")
                if len(items) >= 3 and items[0] in samples:
                    collection.update({"task_id": self.bind_object.sheet.id, "old_name": items[0]},
                                      {"$set": {"productive_name": items[1], "mj_name": items[2]}}, upsert=False,
                                      multi=True)
                elif len(items) >= 2 and items[0] in samples:
                    collection.update({"task_id": self.bind_object.sheet.id, "old_name": items[0]},
                                      {"$set": {"productive_name": items[1]}}, upsert=False, multi=True)

    @report_check
    def get_spname_spid(self):
        if not self.sample_ids:
            self.bind_object.set_error("样本id列表为空，请先调用add_samples_info产生sg_speciem的id", code="52000703")
        collection = self.db["sg_specimen"]
        spname_spid = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_spid[results['new_name']] = id_
        return spname_spid

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
            spname_spid = self.get_spname_spid()
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
        #spcecimen_names = list()
        #group_sorted = list()
        #for key in group_list.keys():
            #lst = group_list[key]
            #tmp_dict = dict()
            #for sample in lst:
                #sp_alias = spname_alias[sample]
                #tmp_dict[sp_alias] = sp_alias
                #tmp_dict[sample] = sample
                ##更新样品名的分组信息
                #col = self.db["sg_specimen"]
                #col.update({"task_id" : self.bind_object.sheet.id, "alias":sp_alias, "about_qc":"before"}, {"$set": {"group": key}})
                #col.update({"task_id" : self.bind_object.sheet.id, "alias":sp_alias, "about_qc":"after"}, {"$set": {"group": key}})
            #spcecimen_names.append(tmp_dict)
            #group_sorted.append(key)
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
            self.bind_object.set_error("样本id列表为空，请先调用add_samples_info产生sg_speciem的id", code="52000704")
        collection = self.db["sg_specimen"]
        spname_id = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_id[results['new_name']] = id_
        return spname_id

    def add_bam_path(self,dir_path):
        spname_spid = self.get_spname_spid()
        col = self.db["sg_specimen"]
        for spname in spname_spid:
            sp_id = spname_spid[spname]
            bam_path = dir_path + "/Bam/" + spname + ".bam"
            insert_data = {"bam_path": bam_path}
            col.update({"_id": sp_id}, {"$set": insert_data})

        # col = self.db["sg_specimen"]
        # results=col.find({"task_id":task_id,"about_qc":"after"})
        # for i in results:
        #     spname=i["old_name"]
        #     sp_id=i["_id"]
        #     bam_path = dir_path +"/"+ spname + ".bam"
        #     insert_data = {"bam_path": bam_path}
        #     col.update({"_id": sp_id}, {"$set": insert_data})

    def get_spname_alias(self):##获取name-alias键值对
        if not self.sample_ids:
            self.bind_object.set_error("样本id列表为空，请先调用add_samples_info产生sg_speciem的id", code="52000705")
        collection = self.db["sg_specimen"]
        spname_alias = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_alias[results['new_name']] = results['alias']
        return spname_alias


    def run1(self):
        self.bind_object.logger.info("开始质控导表")
        self.add_samples_info("/mnt/ilustre/users/sanger-dev/workspace/20190731/Single_hiseq_reads_stat_3763_6524/HiseqReadsStat/output", fq_type="PE",about_qc="before")
        self.add_gragh_info("/mnt/ilustre/users/sanger-dev/workspace/20190731/Single_hiseq_reads_stat_3763_6524/HiseqReadsStat/output/qualityStat", "before")
        self.add_samples_info("/mnt/ilustre/users/sanger-dev/workspace/20190731/Single_hiseq_reads_stat_8237_4023/HiseqReadsStat/output",fq_type="PE", about_qc="after")
        self.add_gragh_info("/mnt/ilustre/users/sanger-dev/workspace/20190731/Single_hiseq_reads_stat_8237_4023/HiseqReadsStat/output/qualityStat", "after")
        self.group_id, self.group_detail, self.group_category = self.add_specimen_group("/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/remote_input/group/group.txt")
        self.bind_object.logger.info(self.group_id)
        self.add_control_group("/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/remote_input/control/control.txt", self.group_id)

    def run3(self):
        self.add_bam_path("/mnt/ilustre/users/sanger-dev/workspace/20190814/Single_denovo_snp8413yyyy/Snp/bwa",task_id="denovo_rna_v2_upgrade")

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    def test(self):
        from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "denovo_rna_v2_upgrade",
            "project_sn": "denovo_rna_v2_upgrade",
            "type": "workflow",
            "name": "denovo_rna_v2.denovo_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = DenovoTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("denovo_rna_v2.denovo_rna_qc")
        wf.test_api.run3()

if __name__ == '__main__':
    unittest.main()

