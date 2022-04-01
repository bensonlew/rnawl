# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# date: 20191009

from __future__ import division

import datetime
import glob
import json
import os
import types
import unittest
from collections import OrderedDict

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId
from bson.son import SON

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Qc(ApiBase):
    def __init__(self, bind_object):
        super(Qc, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    @report_check
    def add_sample_info(self, sample_list, library, group_file=None, productive_table=None):
        """
        对样本信息进行导表
        """
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
        def sort_data_list(data_list, sample_order):
            data_dict = {data['old_name']: data for data in data_list}
            return [data_dict[sample] for sample in sample_order]

        if group_file:
            sample_order = list()
            for line in open(group_file):
                if line.strip() and line[0] != '#':
                    sample_order.append(line.strip().split('\t')[0])
        if not os.path.exists(sample_list):
            raise Exception("%s文件不存在" % sample_list)
        with open(sample_list, "r") as f:
            data_list = []
            samples = []
            for line in f:
                items = line.strip().split("\t")
                if len(items) >= 2 and items[1] != "":
                    sample_name = items[1]
                    if sample_name not in samples:
                        data = {
                            "project_sn": self.bind_object.sheet.project_sn,
                            "task_id": self.bind_object.sheet.id,
                            "library": library,
                            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            "old_name": sample_name,
                            "new_name": sample_name,
                            "group": "",
                            "params": "none",
                            "status": "end",
                            "desc": "样本信息表",
                            'version': 'v1'
                        }
                        if productive_table and sample_name in productive_names:
                            data.update({'productive_name': productive_names[sample_name]})
                        if productive_table and sample_name in mj_names:
                            data.update({'mj_name': mj_names[sample_name]})
                        samples.append(sample_name)
                        data_list.append(data)
        collection = self.db["specimen"]
        if group_file:
            self.bind_object.logger.debug('before sort data_list -> {}'.format(data_list))
            self.bind_object.logger.debug('sample_order -> {}'.format(sample_order))
            data_list = sort_data_list(data_list, sample_order)
            self.bind_object.logger.debug('after sort data_list -> {}'.format(data_list))
        result = collection.insert_many(data_list)
        self.bind_object.logger.info("导入%s文库样品信息数据成功" % library)
        self.sample_ids = result.inserted_ids
        for id in self.sample_ids:
            collection.update({"_id": id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": id}}, upsert=True)
        return samples

    @report_check
    def add_sample_group(self, group_file, library):
        category_names = list()
        specimen_names = list()
        group_dict = OrderedDict()
        with open(group_file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                group_dict.setdefault(tmp[1], list())
                if tmp[0] not in group_dict[tmp[1]]:
                    group_dict[tmp[1]].append(tmp[0])
        specimen_col = self.db["specimen"]
        for key in group_dict:
            category_names.append(key)
            specimen_names.append(group_dict[key])
            for sample in group_dict[key]:
                specimen_col.update({"task_id": self.bind_object.sheet.id, "old_name": sample, "library": library},
                                    {"$set": {"group": key}})

        data = {
            "task_id": self.bind_object.sheet.id,
            "library": library,
            "category_names": category_names,
            "specimen_names": specimen_names,
            "group_name": os.path.basename(group_file),
            "project_sn": self.bind_object.sheet.project_sn,
            "is_use": 1,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": "none",
            "status": "end",
            "desc": "样本分组信息表",
            'version': 'v1'
        }
        group_col = self.db["specimen_group"]
        try:
            group_id = group_col.insert_one(data).inserted_id
            group_col.update({"_id": group_id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": group_id}},
                             upsert=True)
        except Exception, e:
            self.bind_object.set_error("导入%s文库样本分组信息出错:%s" % (library, e))
        else:
            self.bind_object.logger.info("导入%s文库样本分组信息成功" % library)
            return group_id, specimen_names, category_names

    @report_check
    def add_group_compare(self, compare_file, library, group_id):
        con_list = list()
        with open(compare_file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                string = tmp[0] + "|" + tmp[1]
                con_list.append(string)
        group_col = self.db["specimen_group"]
        result = group_col.find_one({"_id": group_id})
        group_name = result["group_name"]
        data = {
            "task_id": self.bind_object.sheet.id,
            "project_sn": self.bind_object.sheet.project_sn,
            "library": library,
            "compare_group_name": group_name,
            "compare_names": json.dumps(con_list),
            "compare_category_name": "all",
            "specimen_group_id": str(group_id),
            "is_use": 1,
            "params": "none",
            "status": "end",
            "desc": "样本对照组信息表",
            'version': 'v1'
        }
        compare_col = self.db["specimen_group_compare"]
        try:
            com_id = compare_col.insert_one(SON(data)).inserted_id
            compare_col.update({"_id": com_id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": com_id}},
                               upsert=True)
        except Exception, e:
            self.bind_object.set_error("导入%s文库样本对照组信息出错:%s" % (library, e))
        else:
            self.bind_object.logger.info("导入%s文库样本对照组信息成功" % library)
            return com_id, con_list

    @report_check
    def add_qc(self, fq_type, library):
        insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            "library": library,
            'params': "none",
            'status': 'end',
            'desc': '质控统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'type': fq_type,
            'version': 'v1'
        }
        collection = self.db['qc']
        try:
            qc_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            self.bind_object.set_error("导入%s文库质控主表出错:%s" % (library, e))
        else:
            self.bind_object.logger.info("导入%s文库质控主表成功" % (library))
            collection.update({"_id": qc_id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": qc_id}},
                              upsert=True)
            return qc_id

    @report_check
    def add_qc_detail(self, qc_main_id, qc_stat_before, qc_stat_after, library, qc_result_dir=None, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        if not isinstance(qc_main_id, ObjectId):
            if isinstance(qc_main_id, types.StringTypes):
                qc_main_id = ObjectId(qc_main_id)
            else:
                raise Exception('qc_main_id必须为ObjectId对象或其对应的字符串！')
        before_stat_file = qc_stat_before + "/fastq_stat.xls"
        after_stat_file = qc_stat_after + "/fastq_stat.xls"
        if not os.path.exists(before_stat_file):
            raise Exception("%s文件不存在" % before_stat_file)
        if not os.path.exists(after_stat_file):
            raise Exception("%s文件不存在" % after_stat_file)
        data_list = []
        data = {}
        before_stat_pd = pd.read_table(before_stat_file, header=0, index_col=0)
        after_stat_pd = pd.read_table(after_stat_file, header=0, index_col=0)
        before_stat_pd_info = before_stat_pd[["Total_Reads", "Total_Bases"]]
        after_stat_pd_info = after_stat_pd[["Total_Reads", "Total_Bases", "Error%", "Q20%", "Q30%", "GC%"]]
        if library == "small":
            if qc_result_dir is None:
                raise Exception('smallrna文库必须提供质控结果目录！')
            small_qc_stat_files = glob.glob(qc_result_dir + "/*_qc_stat.xls")
            small_qc_stat_file = qc_result_dir + "/qc_stat.xls"
            with open(small_qc_stat_file, "w") as w:
                num = 0
                for file in small_qc_stat_files:
                    num += 1
                    with open(file, "r") as f:
                        if num > 1:
                            f.readline()
                        for line in f:
                            w.write(line)
            small_qc_pd = pd.read_table(small_qc_stat_file, header=0, index_col=0)
            small_qc_info = small_qc_pd[["Clean_reads"]]
            stat_final_pd = pd.concat([before_stat_pd_info, after_stat_pd_info, small_qc_info], axis=1)
        else:
            dup_file = qc_stat_after + "/dup.xls"
            rfam_file = qc_stat_after + "/stat_results"
            if not os.path.exists(dup_file):
                raise Exception("%s文件不存在" % dup_file)
            if not os.path.exists(rfam_file):
                raise Exception("%s文件不存在" % rfam_file)
            dup_pd = pd.read_table(dup_file, header=0, index_col=0)
            rfam_pd = pd.read_table(rfam_file, header=0, index_col=0)
            if 'PairedDup' in dup_pd.columns:
                fastq_type = "PE"
                dup_pd_info = dup_pd[["read1Dup", "read2Dup", "PairedDup"]]
            else:
                fastq_type = "SE"
                dup_pd_info = dup_pd[["readDup"]]
            rfam_pd_info = rfam_pd[["rRNA Ratio(%)"]]
            stat_final_pd = pd.concat([before_stat_pd_info, after_stat_pd_info, dup_pd_info, rfam_pd_info], axis=1)
        header = ['Sample', 'Raw Reads', 'Raw Bases', 'Clean reads', 'Clean bases', 'Error rate(%)', 'Q20(%)', 'Q30(%)',
                  'GC content(%)']
        if library == "small":
            header.extend(['Useful reads(18-32nt)'])
        else:
            if fastq_type == "SE":
                header.extend(['readDup', 'rRNA Ratio(%)'])
            else:
                header.extend(['Dup_R1', 'Dup_R2', 'Dup_Pair', 'rRNA Ratio(%)'])
        stat_final = qc_stat_after + "/fastq_stat_final.xls"
        with open(stat_final, "w") as w:
            w.write("\t".join(header) + "\n")
        stat_final_pd.to_csv(stat_final, header=False, index=True, sep='\t', mode='a')
        with open(stat_final, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "sample": line[0],
                    "raw_reads": int(line[1]),
                    "raw_bases": int(line[2]),
                    "clean_reads": int(line[3]),
                    "clean_bases": int(line[4]),
                    "error_rate": round(float(line[5]), 4),
                    "q20_rate": round(float(line[6]), 2),
                    "q30_rate": round(float(line[7]), 2),
                    "gc_rate": round(float(line[8]), 2),
                    "qc_id": qc_main_id,
                }
                if library == "small":
                    data.update({
                        "useful_reads": int(line[9]),
                    })
                    if qc_result_dir is None:
                        raise Exception('smallrna文库必须提供质控结果目录！')
                    with open(qc_result_dir + "/" + line[0] + "_clean.length.txt", 'r') as f1:
                        len_stat = [{line.strip().split()[0]: int(line.strip().split()[1])} for line in f1.readlines()]
                        data.update({
                            "len_data": len_stat,
                        })
                else:
                    if fastq_type == "SE":
                        data.update({
                            "dup": round(float(line[9]), 4) * 100,
                            "rrna_ratio": round(float(line[10]), 4)
                        })
                    else:
                        data.update({
                            "dup_r1": round(float(line[9]), 4) * 100,
                            "dup_r2": round(float(line[10]), 4) * 100,
                            "dup_pair": round(float(line[11]), 4) * 100,
                            "rrna_ratio": round(float(line[12]), 4)
                        })
                data_list.append(data)
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['sample']))
        try:
            collection = self.db["qc_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s文库样品质控信息数据出错:%s" % (library, e))
        else:
            self.bind_object.logger.info("导入%s文库样品质控信息数据成功" % library)

    @report_check
    def add_qc_graph(self, qc_main_id, qc_stat_dir, library, about_qc, qc_result_dir=None):
        if not isinstance(qc_main_id, ObjectId):
            if isinstance(qc_main_id, types.StringTypes):
                qc_main_id = ObjectId(qc_main_id)
            else:
                raise Exception('qc_main_id必须为ObjectId对象或其对应的字符串！')
        stat_files = sorted(glob.glob("{}/qualityStat/*".format(qc_stat_dir)))
        data_list = []
        for sf in stat_files:
            sample_name = os.path.basename(sf).split(".")[0]
            self.bind_object.logger.info('%s,%s' % (sf, sample_name))
            spname_spid = self.get_spname_spid()
            site = os.path.basename(sf).split(".")[1]
            if site == "l":
                site_type = "left"
            elif site == "r":
                site_type = "right"
            else:
                site_type = "single"
            with open(sf, "r") as f:
                f.readline()
                for line in f:
                    line = line.strip().split()
                    total = int(line[12]) + int(line[13]) + int(line[14]) + int(line[15]) + int(line[16])
                    data = {
                        # "project_sn": self.bind_object.sheet.project_sn,
                        # "task_id": self.bind_object.sheet.id,
                        "specimen_name": sample_name,
                        "specimen_id": spname_spid[sample_name],
                        "qc_id": qc_main_id,
                        "type": site_type,
                        "about_qc": about_qc,
                        "column": int(line[0]),
                        "min": int(line[10]),
                        "max": int(line[11]),
                        "q1": int(line[6]),
                        "q3": int(line[8]),
                        "median": int(line[7]),
                        "error": 10 ** (float(line[5]) / (-10)) * 100,
                        "A": int(line[12]) / total * 100,
                        "T": int(line[15]) / total * 100,
                        "C": int(line[13]) / total * 100,
                        "G": int(line[14]) / total * 100,
                        "N": int(line[16]) / total * 100,
                        # 'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    }
                    # if library == "small" and about_qc == "after":
                    #     if qc_result_dir is None:
                    #         raise Exception('smallrna文库必须提供质控结果目录！')
                    #     with open(qc_result_dir + "/" + sample_name + "_clean.length.txt", 'r') as f1:
                    #         len_stat = [{line.strip().split()[0]: int(line.strip().split()[1])} for line in f1.readlines()]
                    #         data.update({
                    #             "len_data": len_stat,
                    #         })
                    data_list.append(data)
        try:
            collection = self.db["qc_graph"]
            result = collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入{}文库样品{}质控画图数据信息出错:{}".format(library, about_qc, e))
        else:
            self.bind_object.logger.info("导入{}文库样品{}质控画图数据信息成功".format(library, about_qc))

    @report_check
    def get_spname_spid(self):
        if not self.sample_ids:
            raise Exception("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
        collection = self.db["specimen"]
        spname_spid = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_spid[results['new_name']] = id_
        return spname_spid

    @report_check
    def add_qc_detail_after(self, qc_id, stat_output_dir, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        before_stat_file = None
        data = list()
        document = dict()
        fastq_stat_df = pd.read_table(os.path.join(stat_output_dir, 'fastq_stat.xls'), header=0, index_col=0)
        fastq_stat_info_df = fastq_stat_df[['Total_Reads', 'Total_Bases', 'Error%', 'Q20%', 'Q30%', 'GC%']]
        dup_df = pd.read_table(os.path.join(stat_output_dir, 'dup.xls'), header=0, index_col=0)
        dup_info_df = dup_df[['read1Dup', 'read2Dup', 'PairedDup']]
        rfam_df = pd.read_table(os.path.join(stat_output_dir, 'stat_results'), header=0, index_col=0)
        rfam_info_df = rfam_df[['rRNA Ratio(%)']]
        fastq_stat_final_table = os.path.join(stat_output_dir, 'fastq_stat_final.xls')
        stat_df = pd.concat([fastq_stat_info_df, dup_info_df, rfam_info_df], axis=1).reset_index()
        stat_df.columns = ['Sample', 'Raw Reads', 'Raw Bases', 'Clean reads', 'Clean bases', 'Error rate(%)', 'Q20(%)',
                           'Q30(%)', 'GC content(%)', 'Dup_R1', 'Dup_R2', 'Dup_Pair', 'rRNA Ratio(%)']
        stat_df.to_csv(fastq_stat_final_table, sep='\t', index=False)
        with open(fastq_stat_final_table) as handle:
            handle.readline()
            for line in handle:
                eles = line.strip().split('\t')
                document = {
                    'qc_id': qc_id,
                    'sample': eles[0],
                    'raw_reads': int(eles[1]),
                    'raw_bases': int(eles[2]),
                    'clean_reads': int(eles[3]),
                    'clean_bases': int(eles[4]),
                    'error_rate': round(float(eles[5]), 4),
                    'q20_rate': round(float(eles[6]), 2),
                    'q30_rate': round(float(eles[7]), 2),
                    'gc_rate': round(float(eles[8]), 2),
                    'dup_r1': round(float(eles[9]), 4),
                    'dup_r2': round(float(eles[10]), 4),
                    'dup_pair': round(float(eles[11]), 4),
                    'rrna_ratio': round(float(eles[12]), 4)
                }
                data.append(document)
        if group:
            data.sort(key=lambda x: sample_list.index(x['sample']))
        self.create_db_table('qc_detail', data)

    def export_small_qc(self):
        ## smallRNA文库导表
        small_sample_list = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/remote_input/fastq_dir/rawdata/list.txt"  # 样本信息列表
        small_group_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/remote_input/group_table/group.txt"  # 分组方案表
        small_control_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/remote_input/control_file/control.txt"  # 对照组方案表
        self.add_sample_info(sample_list=small_sample_list, library="small")
        small_group_id, small_group_detail, small_group_category = self.add_sample_group(group_file=small_group_table,
                                                                                         library="small")
        small_control_id, small_control_detail = self.add_group_compare(compare_file=small_control_table,
                                                                        library="small", group_id=small_group_id)
        small_before_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/HiseqReadsStat/output"  # raw数据统计结果
        small_after_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/HiseqReadsStat1/output"  # clean数据统计结果
        qc_result = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/MirnaQc/output/clean_data"
        small_qc_id = self.add_qc(fq_type="SE", library="small")
        self.add_qc_detail(small_qc_id, small_before_qc_stat, small_after_qc_stat, "small", qc_result_dir=qc_result)
        self.add_qc_graph(small_qc_id, small_before_qc_stat, "small", "before")
        self.add_qc_graph(small_qc_id, small_after_qc_stat, "small", "after", qc_result_dir=qc_result)

    def export_long_qc(self):
        ## longRNA文库导表
        long_sample_list = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/fastq_dir/rawdata/list.txt"  # 样本信息列表
        long_group_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/group_table/group.txt"  # 分组方案表
        long_control_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/control_file/control.txt"  # 对照组方案表
        self.add_sample_info(sample_list=long_sample_list, library="long")
        long_group_id, long_group_detail, long_group_category = self.add_sample_group(group_file=long_group_table,
                                                                                      library="long")
        long_control_id, long_control_detail = self.add_group_compare(compare_file=long_control_table, library="long",
                                                                      group_id=long_group_id)
        long_before_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/HiseqReadsStat/output"  # raw数据统计结果
        long_after_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/HiseqReadsStat1/output"  # clean数据统计结果
        long_qc_id = self.add_qc(fq_type="PE", library="long")
        self.add_qc_detail(long_qc_id, long_before_qc_stat, long_after_qc_stat, "long")
        self.add_qc_graph(long_qc_id, long_before_qc_stat, "long", "before")
        self.add_qc_graph(long_qc_id, long_after_qc_stat, "long", "after")

    def export_circle_qc(self):
        ## circRNA文库导表
        circ_sample_list = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/fastq_dir/rawdata/list.txt"  # 样本信息列表
        circ_group_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/group_table/group.txt"  # 分组方案表
        circ_control_table = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/remote_input/control_file/control.txt"  # 对照组方案表
        self.add_sample_info(sample_list=circ_sample_list, library="circle")
        circle_group_id, circle_group_detail, circle_group_category = self.add_sample_group(group_file=circ_group_table,
                                                                                            library="circle")
        circle_control_id, circle_control_detail = self.add_group_compare(compare_file=circ_control_table,
                                                                          library="circle", group_id=circle_group_id)
        circle_before_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/HiseqReadsStat/output"  # raw数据统计结果
        circle_after_qc_stat = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/HiseqReadsStat1/output"  # clean数据统计结果
        circle_qc_id = self.add_qc(fq_type="PE", library="circle")
        self.add_qc_detail(circle_qc_id, circle_before_qc_stat, circle_after_qc_stat, "circle")
        self.add_qc_graph(circle_qc_id, circle_before_qc_stat, "circle", "before")
        self.add_qc_graph(circle_qc_id, circle_after_qc_stat, "circle", "after")

    def run1(self):
        group_table = "/mnt/ilustre/users/sanger-dev/workspace/20190717/Refrna_tsg_34883/remote_input/group_table/default_group.txt"
        control_table = "/mnt/ilustre/users/sanger-dev/workspace/20190717/Refrna_tsg_34883/remote_input/control_file/control_file.txt"
        group_id, group_detail, circle_group_category = self.add_sample_group(group_file=group_table, library="circle")
        control_id, control_detail = self.add_group_compare(compare_file=control_table, library="circle",
                                                            group_id=group_id)


################################################
class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test(test):
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.api_qc = wf.api.api("whole_transcriptome.qc")

        wf.api_qc.export_small_qc()
        # wf.api_qc.export_long_qc()
        # wf.api_qc.run1()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
