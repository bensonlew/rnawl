# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,shicaiping'

from __future__ import division

import datetime
import glob
import json
import os
import re
import unittest

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.ref_rna_v2.api_base import ApiBase


class Mapping(ApiBase):
    def __init__(self, bind_object):
        super(Mapping, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    @report_check
    def add_bam_path(self, dir_path):
        """
        将bam文件的路径插入sg_specimen表中，供可变剪切使用
        :param dir_path:传入的rnaseq_mapping的output_dir
        :return:
        """
        spname_spid = self.get_spname_spid()
        col = self.db["specimen"]
        for spname in spname_spid:
            sp_id = spname_spid[spname]
            bam_path = dir_path + "/Align/AlignBam/" + spname + ".bam"
            insert_data = {"bam_path": bam_path}
            col.update({"_id": sp_id}, {"$set": insert_data})

    @report_check
    def get_spname_spid(self):
        if not self.sample_ids:
            self.bind_object.set_error("样本id列表为空，请先调用add_samples_info产生specimen的id", code="53702006")
        collection = self.db["specimen"]
        spname_spid = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_spid[results['new_name']] = id_
        return spname_spid

    @report_check
    def add_mapping_stat(self, stat_file, library, method, sample_list_file=None, group=None):
        name = library + "rna_" + method + "_mapping"
        params = json.dumps(
            {'task_id': self.bind_object.sheet.id, 'submit_location': 'mapping', 'method': method, 'library': library},
            sort_keys=True)
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "params": params,
            "library": library,
            "method": method,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "start",
            "desc": '比对分析结果主表',
            "name": name,
            'version': 'v1'
        }
        try:
            collection = self.db['mapping']
            mapping_id = collection.insert_one(data).inserted_id
        except Exception, e:
            self.bind_object.set_error("导入%s文库比对信息主表出错:%s" % (library, e))
        else:
            self.bind_object.logger.info("导入%s文库比对信息主表成功" % library)
            collection.update({"_id": mapping_id, "task_id": self.bind_object.sheet.id},
                              {"$set": {"main_id": mapping_id, "status": "end"}}, upsert=True)

            if method.lower() == "hisat":
                self.add_hisat_mapping_stat(stat_file, mapping_id, group)
            elif method.lower() == "tophat":
                self.add_tophat_mapping_stat(stat_file, mapping_id, group)
            elif method.lower() == "star":
                self.add_star_mapping_stat(stat_file, mapping_id, group)
            else:
                if library == "small":

                    if group:
                        sample_list = list()
                        with open(group, 'r') as g:
                            for line in g.readlines():
                                if line.startswith('#'):
                                    continue
                                sample_list.append(line.strip().split('\t')[0])
                    else:
                        if os.path.exists(sample_list_file):
                            sample_list = list()
                            with open(sample_list_file, 'r') as f:
                                for line in f.readlines():
                                    sample = line.strip().split("\t")[1]
                                    if sample in sample_list:
                                        pass
                                    else:
                                        sample_list.append(sample)
                    collection.update({"_id": mapping_id, "task_id": self.bind_object.sheet.id},
                          {"$set": {"sample_list": sample_list}}, upsert=True)
                    self.add_small_mapping_stat(mapping_id=mapping_id, map_dir=stat_file, sample_list=sample_list, group=group)


    @report_check
    def add_rpkm_table(self, file_path, name=None, params='RSeQC-2.6.3', detail=True, library=None):
        """
        :param file_path: 文件夹，即~/MapAssessment/output/satur
        :param name: 主表的名称，可不传
        :param params: 参数，可不传
        :param detail: 如果需要主表跟detail表一起导入传true,否则false
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "saturation_origin",
            "library": library,
            "status": "end",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "curve_category": ['5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75',
                               '80', '85', '90', '95', '100'],
            "curve_specimen": {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)",
                               "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'version': 'v1'
        }
        collection = self.db["assessment_saturation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_rpkm_curve(file_path, inserted_id)
            self.update_db_record('assessment_saturation', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_rpkm_curve(self, rpkm_file, rpkm_id=None):
        """
        :param rpkm_file:文件夹，即~/MapAssessment/output/satur
        :param rpkm_id:主表id，必传
        :return:
        """
        rpkm_id = ObjectId(rpkm_id)
        curve_files = sorted(glob.glob("{}/*cluster_percent.xls".format(rpkm_file)))
        R_files = sorted(glob.glob("{}/*eRPKM.xls.saturation.R".format(rpkm_file)))
        sample_categaries = {}
        for rf in R_files:
            sample_name = os.path.basename(rf).split(".")[0][6:]
            categaries = self.add_satur_count(rf)
            sample_categaries[sample_name] = categaries
        curve_data = []
        for cf in curve_files:
            sample_name = os.path.basename(cf).split(".")[0][6:]
            with open(cf, "r") as f:
                line_list = []
                for line in f:
                    line = line.strip().split()
                    line.pop(0)
                    new_line = [float(x) for x in line]
                    line_list.append(new_line)
                data = {
                    "saturation_id": rpkm_id,
                    "specimen_name": sample_name,
                    "categories": sample_categaries[sample_name],
                    "column1": line_list[0],
                    "column2": line_list[1],
                    "column3": line_list[2],
                    "column4": line_list[3],
                    "column5": line_list[4],
                    "column6": line_list[5]
                }
                # self.bind_object.set_error data
                curve_data.append(data)
        try:
            collection = self.db["assessment_saturation_detail"]
            collection.insert_many(curve_data)
        except Exception, e:
            self.bind_object.set_error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_satur_count(self, count_r_file):
        categaries = {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)",
                      "column5": "[15-60)", "column6": ">=60"}
        with open(count_r_file, "r") as f:
            for line in f:
                if re.match(r"legend", line):
                    all_num = re.findall("num=[\d]*", line)
                    # print all_num
                    categaries = {
                        "column5": "[15-60)=" + all_num[4][4:],
                        "column4": "[3.5-15)=" + all_num[3][4:],
                        "column6": ">=60=" + all_num[5][4:],
                        "column1": "[0-0.3)=" + all_num[0][4:],
                        "column3": "[0.6-3.5)=" + all_num[2][4:],
                        "column2": "[0.3-0.6)=" + all_num[1][4:]
                    }
        return categaries

    @report_check
    def add_coverage_table(self, coverage, name=None, params='RSeQC-2.6.3', detail=True, library=None):
        """
        :param coverage: 文件夹，即~/MapAssessment/output/coverage
        :param name: 主表名称。可不传
        :param params:参数，可不传
        :param detail:如果需要主表跟detail表一起导入传true,否则false
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "coverage_origin",
            "status": "end",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "library": library,
            'version': 'v1'
        }
        collection = self.db["assessment_coverage"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_coverage_detail(coverage, inserted_id)
            self.update_db_record('assessment_coverage', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_coverage_detail(self, coverage, coverage_id=None):
        """
        :param coverage: 文件夹，即~/MapAssessment/output/coverage
        :param coverage_id: 主表ID，必传
        :return:
        """
        coverage_files = sorted(glob.glob("{}/*geneBodyCoverage.txt".format(coverage)))
        data_list = []
        for cf in coverage_files:
            sample_name = os.path.basename(cf).split(".")[0]
            if sample_name.startswith("coverage_"):
                sample_name = sample_name.split("coverage_")[1]
            # self.bind_object.set_error sample_name
            with open(cf, "r") as f:
                percent = f.readline().strip().split()
                value = f.next().strip().split()
                plot_value = {}
                for i in range(100):
                    plot_value[percent[i + 1]] = value[i + 1]
                data = {
                    "coverage_id": coverage_id,
                    "specimen_name": sample_name,
                    "plot_value": plot_value
                }
                data_list.append(data)
        try:
            collection = self.db["assessment_coverage_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_distribution_table(self, distribution=None, name=None, params='RSeQC-2.6.3', library=None, group=None):
        """
        :param distribution: 文件夹，即~/MapAssessment/output/distribution,不传的时候只导主表
        :param name: 主表名称。可不传
        :param params:参数，可不传
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "distribution_origin",
            "status": "end",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "library": library,
            'version': 'v1'
        }
        collection = self.db["assessment_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_distribution_detail(distribution, inserted_id, group)
            self.update_db_record('assessment_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_distribution_detail(self, distribution, distribution_id, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        data_list = []
        stat_files = sorted(glob.glob("{}/*reads_distribution.txt".format(distribution)))
        distributions = ["cds", "utr5", "utr3", "introns", "intergenic"]
        for fls in stat_files:
            sample_name = os.path.basename(fls).split(".")[0]
            values = []
            with open(fls, "r") as f:
                f.readline()
                f.next()
                f.next()
                f.next()
                f.next()
                data = {
                    "distribution_id": distribution_id,
                    "specimen_name": sample_name
                }
                for line in f:
                    if re.match(r"==", line):
                        continue
                    else:
                        line = line.strip().split()
                        values.append(float(line[2]))
                # print values
                # print sum(values[4:])
                values_new = values[:4]
                values_new.append(sum([values[6], values[9]]))
                # print values_new
                total = sum(values_new)
                for n, dis in enumerate(distributions):
                    # data[dis] = values_new[n]
                    data[dis] = str(values_new[n]) + "(" + str(float("%0.4f" % (values_new[n] / total)) * 100) + "%)"
                print data
                data_list.append(data)
            # print data_list
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["assessment_distribution_detail"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入reads distribution详情表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入reads distribution详情表成功")

    @report_check
    def add_chrom_distribution_table(self, distribution=None, name=None, params='RSeQC-2.6.3', library=None,
                                     sample_list_file=None, group=None):
        """
        :param distribution: 文件夹，即~/MapAssessment/output/chr_stat,不传的时候只导主表
        :param name: 主表名称。可不传
        :param params:参数，可不传
        :return:
        """
        # 改变染色体分布导表方式，原因在于需要对染色体丰度进行排序
        if library == "small":
            pass
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "distribution_origin",
            "status": "end",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "library": library,
            'version': 'v1'
        }
        collection = self.db["assessment_chrom_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_chrom_distribution_detail(distribution, inserted_id, library, sample_list_file, group)
            self.update_db_record('assessment_chrom_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_chrom_distribution_detail(self, chorm_stat, chorm_stat_id, library, sample_list_file=None, group=None):
        data_list = []
        distribution = set()
        specimen_names = []
        if group:
            sample_list_ = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list_.append(line.strip().split('\t')[0])
        if library == "small":
            if sample_list_file is None or not os.path.exists(sample_list_file):
                raise Exception("导入samll文库染色体reads分布详情表时必须提供样本列表")
            sample_list = list()
            with open(sample_list_file, 'r') as f:
                for line in f.readlines():
                    sample = line.strip().split("\t")[1]
                    if sample in sample_list:
                        pass
                    else:
                        sample_list.append(sample)
            for sample in sample_list:
                stat_file = glob.glob("{}/{}_map_stat.xls".format(chorm_stat, sample))[0]
                chrs = []
                specimen_names.append(sample)
                with open(stat_file, "r") as f:
                    f.readline()
                    f.readline()
                    f.readline()
                    data = {
                        "chrom_distribution_id": chorm_stat_id,
                        "specimen_name": sample,
                        "chr_values": []
                    }
                    for line in f:
                        values = {}
                        line = line.strip().split()
                        chrs.append(line[0])
                        distribution.add(line[0])
                        values["chr_name"] = line[0]
                        values["value"] = line[-4]
                        values["forword"] = line[-2]
                        values["reverse"] = line[-1]
                        data["chr_values"].append(values)
                    data_list.append(data)
        else:
            stat_files = sorted(glob.glob("{}/*.bam_chr_stat.xls".format(chorm_stat)))
            for fls in stat_files:
                chrs = []
                # chr_values = []
                sample_name = os.path.basename(fls).split(".")[0]
                specimen_names.append(sample_name)
                with open(fls, "r") as f:
                    data = {
                        "chrom_distribution_id": chorm_stat_id,
                        "specimen_name": sample_name,
                        "chr_values": []
                    }
                    for line in f:
                        values = {}
                        if re.match(r"#", line):
                            continue
                        else:
                            line = line.strip().split()
                            chrs.append(line[0])
                            distribution.add(line[0])
                            values["chr_name"] = line[0]
                            values["value"] = line[1]
                            data["chr_values"].append(values)
                    # print data
                    data_list.append(data)
                    # print chrs
        try:
            collection = self.db["assessment_chrom_distribution_detail"]
            result = collection.insert_many(data_list)
            main_collection = self.db["assessment_chrom_distribution"]
            if group:
                main_collection.update({"_id": ObjectId(chorm_stat_id)},
                                       {"$set": {"distribution": list(distribution), "specimen_name": sample_list_}})
            else:
                main_collection.update({"_id": ObjectId(chorm_stat_id)},
                                       {"$set": {"distribution": list(distribution), "specimen_name": specimen_names}})
        except Exception, e:
            self.bind_object.logger.info("导入assessment_chrom_distribution_detail出错:%s" % e)
        else:
            self.bind_object.logger.info("导入assessment_chrom_distribution_detail成功")

    @report_check
    def add_tophat_mapping_stat(self, stat_file, main_mapping_id, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        files = sorted(glob.glob("{}/*".format(stat_file)))
        data_list = []
        for fs in files:
            specimen_name = os.path.basename(fs).split(".")[0]
            print specimen_name
            f = open(fs, "r")
            data = {
                "specimen_name": specimen_name,
                "mapping_id": main_mapping_id
            }
            total_reads = 0
            map_reads = 0
            multiple = 0
            for line in f:
                if re.match(r"Left", line):
                    total_reads += int(f.next().split()[-1])
                    map_reads += int(f.next().split()[2])
                    multiple += int(f.next().split()[2])
                if re.match(r"Right", line):
                    total_reads += int(f.next().split()[-1])
                    map_reads += int(f.next().split()[2])
                    multiple += int(f.next().split()[2])
            print total_reads, map_reads, multiple, total_reads - multiple
            data["total_reads"] = total_reads
            # print "(" + str(float("%0.4f" % ( map_reads/total_reads)) * 100) + "%" + ")"
            data["mapping_reads"] = str(map_reads) + "(" + str(
                float("%0.4f" % (map_reads / total_reads)) * 100) + "%" + ")"
            data["multiple_mapped"] = str(multiple) + "(" + str(
                float("%0.4f" % (multiple / total_reads)) * 100) + "%" + ")"
            data["uniq_mapped"] = str(total_reads - multiple) + "(" + str(
                float("%0.4f" % ((map_reads - multiple) / total_reads)) * 100) + "%" + ")"
            # print data
            data_list.append(data)
            f.close()
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    @report_check
    def add_hisat_mapping_stat(self, stat_file, main_mapping_id, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        files = sorted(glob.glob("{}/*".format(stat_file)))
        data_list = []
        for fs in files:
            specimen_name = os.path.basename(fs).split(".")[0]
            values = []
            f = open(fs, "r")
            for line in f:
                if re.match(r" ", line):
                    line = line.split()
                    values.append(line[0])
            if len(values) == 13:
                total_reads = int(values[0]) * 2
                unmap_reads = int(values[-3])
                map_reads = int(total_reads) - int(unmap_reads)
                map_reads = str(map_reads) + "(" + str(
                    float("%0.4f" % (map_reads / int(total_reads))) * 100) + "%" + ")"
                uniq_mapped = int(values[2]) * 2 + int(values[6]) * 2 + int(values[-2])
                uniq_mapped = str(uniq_mapped) + "(" + str(
                    float("%0.4f" % (uniq_mapped / int(total_reads))) * 100) + "%" + ")"
                multi_mapped = int(values[3]) * 2 + int(values[-1])
                multi_mapped = str(multi_mapped) + "(" + str(
                    float("%0.4f" % (multi_mapped / int(total_reads))) * 100) + "%" + ")"
                print specimen_name
                print total_reads, map_reads, uniq_mapped, multi_mapped
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads,
                    "multiple_mapped": multi_mapped,
                    "uniq_mapped": uniq_mapped,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            elif len(values) == 4:
                total_reads = int(values[0])
                unmap_reads = int(values[1])
                map_reads = int(total_reads) - int(unmap_reads)
                uniq_mapped = int(values[2])
                multi_mapped = int(values[3])
                map_reads = str(map_reads) + "(" + str(
                    float("%0.4f" % (map_reads / int(total_reads))) * 100) + "%" + ")"
                uniq_mapped = str(uniq_mapped) + "(" + str(
                    float("%0.4f" % (uniq_mapped / int(total_reads))) * 100) + "%" + ")"
                multi_mapped = str(multi_mapped) + "(" + str(
                    float("%0.4f" % (multi_mapped / int(total_reads))) * 100) + "%" + ")"
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads,
                    "multiple_mapped": multi_mapped,
                    "uniq_mapped": uniq_mapped,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            f.close()
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    @report_check
    def add_star_mapping_stat(self, stat_file, main_mapping_id, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        files = sorted(glob.glob("{}/*".format(stat_file)))
        data_list = []
        for fs in files:
            specimen_name = os.path.basename(fs).split(".")[0]
            values = []
            f = open(fs, "r")
            for line in f:
                line = line.strip()
                if re.search("|", line):
                    line = line.split("|")
                    if len(line) > 1:
                        line_info = line[1].split()
                        line_info = line_info[0]
                        values.append(line_info)
            if values[5] > 150:
                total_reads = int(values[4]) * 2
                if not values[-7].endswith("%"):
                    multi_mapped = int(values[-7]) * 2 + int(values[-9]) * 2
                else:
                    multi_mapped = int(values[-10]) * 2 + int(values[-12]) * 2
                uniq_mapped = int(values[6]) * 2
                map_reads = multi_mapped + uniq_mapped
                map_reads_rate = float("%0.4f" % (map_reads / total_reads)) * 100
                uniq_mapped_rate = float("%0.4f" % (uniq_mapped / total_reads)) * 100
                multi_mapped_rate = float("%0.4f" % (multi_mapped / total_reads)) * 100
                multi_mapped_f = str(multi_mapped) + "(" + str(multi_mapped_rate) + "%" + ")"
                map_reads_f = str(map_reads) + "\t" + "(" + str(map_reads_rate) + "%" + ")"
                uniq_mapped_f = str(uniq_mapped) + "\t" + "(" + str(uniq_mapped_rate) + "%" + ")"
                print specimen_name
                print total_reads, map_reads, uniq_mapped, multi_mapped
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads_f,
                    "multiple_mapped": multi_mapped_f,
                    "uniq_mapped": uniq_mapped_f,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            else:
                total_reads = int(values[4])
                if not values[-7].endswith("%"):
                    multi_mapped = int(values[-7]) + int(values[-9])
                else:
                    multi_mapped = int(values[-10]) + int(values[-12])
                uniq_mapped = int(values[6])
                map_reads = multi_mapped + uniq_mapped
                map_reads_rate = float("%0.4f" % (map_reads / total_reads)) * 100
                uniq_mapped_rate = float("%0.4f" % (uniq_mapped / total_reads)) * 100
                multi_mapped_rate = float("%0.4f" % (multi_mapped / total_reads)) * 100
                multi_mapped_f = str(multi_mapped) + "(" + str(multi_mapped_rate) + "%" + ")"
                map_reads_f = str(map_reads) + "\t" + "(" + str(map_reads_rate) + "%" + ")"
                uniq_mapped_f = str(uniq_mapped) + "\t" + "(" + str(uniq_mapped_rate) + "%" + ")"
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads_f,
                    "multiple_mapped": multi_mapped_f,
                    "uniq_mapped": uniq_mapped_f,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            f.close()
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    @report_check
    def add_small_mapping_stat(self, mapping_id=None, map_dir=None, sample_list=None, group=None):
        if group:
            sample_list_ = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list_.append(line.strip().split('\t')[0])
        if not isinstance(mapping_id, ObjectId):
            if isinstance(mapping_id, types.StringTypes):
                mapping_id = ObjectId(mapping_id)
            else:
                raise Exception('mapping_id必须为ObjectId对象或其对应的字符串！')
        data_stat_list = []
        all_stat = pd.read_table(map_dir + "/All_map_stat.xls", header=0, index_col=0,
                                 usecols=['Chromosome', 'Total_num'])
        with open(map_dir + "/Genome_map_stat.xls", 'w') as f:
            f.write("Sample\tTotal reads\tMapped reads\tMapped reads(+)\tMapped reads(-)\n")

            all_stat = all_stat[2:]
            all_stat['ref'] = all_stat.index
            all_stat.rename(str.lower, axis='columns', inplace=True)
            for sample in sample_list:
                sample_stat = pd.read_table(map_dir + "/" + sample + "_map_stat.xls", header=0, index_col=0,
                                            usecols=['Chromosome', 'Total_num', 'Forword.1', 'Reverse.1'])[2:]
                all_stat[sample + '_total'] = sample_stat['Total_num'].map(int)
                all_stat[sample + '_for'] = sample_stat['Forword.1'].map(int)
                all_stat[sample + '_rev'] = sample_stat['Reverse.1'].map(int)

                stat2 = pd.read_table(map_dir + "/" + sample + "_map_stat.xls", header=0, index_col=0)[:2]
                ratio = float(stat2.iloc[1]['Total_num']) / float(stat2.iloc[0]['Total_num'])
                data = dict({
                    "specimen_name": sample,
                    "total_reads": int(stat2.iloc[0]['Total_num']),
                    "mapped_reads": int(stat2.iloc[1]['Total_num']),
                    "mapped_reads_for": int(stat2.iloc[1]['Forword.1']),
                    "mapped_reads_rec": int(stat2.iloc[1]['Reverse.1']),
                    "mapping_id": mapping_id,
                    "ratio": ratio
                })
                data_stat_list.append(data)

                f.write("{}\t{}\t{}\t{}\t{}\n".format(sample,
                                                      int(stat2.iloc[0]['Total_num']),
                                                      int(stat2.iloc[1]['Total_num']),
                                                      stat2.iloc[1]['Forword.1'],
                                                      stat2.iloc[1]['Reverse.1']))
        all_stat['mapping_id'] = mapping_id
        all_stat['type'] = "graph"
        all_stat = all_stat.fillna(0)

        def float2int(x):
            if type(x) == "float":
                return int(x)
            else:
                return x

        all_stat = all_stat.applymap(float2int)
        columns = []
        columns.extend([x + "_total" for x in sample_list])
        columns_rename = []
        columns_rename.extend(sample_list)
        all_stat_choose = all_stat.loc[:, columns]
        rename_dict = dict(zip(columns, columns_rename))
        all_stat_choose = all_stat_choose.rename(columns=rename_dict)
        all_stat_choose.to_csv(map_dir + "/Chro_map_stat.xls", sep="\t", index=True)
        row_dict_list = all_stat.to_dict('records')
        if group:
            data_stat_list.sort(key=lambda x: sample_list_.index(x['specimen_name']))
        try:
            self.create_db_table('mapping_detail', row_dict_list)
            self.create_db_table('mapping_stat', data_stat_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    def export_mapping(self):
        ## 导smallrna文库比对结果
        stat_file = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/MapperAndStat/output"
        sample_list_file = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Smallrna_tsg_35702/MirnaQc/output/clean_data/list.txt"
        library = "small"
        method = "bowtie"
        self.add_mapping_stat(stat_file=stat_file, library=library, method=method, sample_list_file=sample_list_file)
        # 改变染色体分布导表方式，原因在于需要对染色体丰度进行排序
        # self.add_chrom_distribution_table(distribution=stat_file, params='none', library=library,
        #                                  sample_list_file=sample_list_file)

        ## 导longrna文库比对结果
        stat_file = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/RnaseqMapping/output/stat"
        satur_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/saturation"
        coverage_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/coverage"
        chrstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/chr_stat"
        bamstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/distribution"
        library = "long"
        method = "hisat"
        # self.add_mapping_stat(stat_file=stat_file, library=library, method=method)
        # self.add_rpkm_table(satur_dir, params='RSeQC-2.6.3', detail=True, library=library)
        self.add_coverage_table(coverage_dir, params='RSeQC-2.6.3', detail=True, library=library)
        # self.add_distribution_table(distribution=bamstat_dir, params='RSeQC-2.6.3', library=library)
        # self.add_chrom_distribution_table(distribution=chrstat_dir, params='RSeQC-2.6.3', library=library)

        ## 导circle文库比对结果
        stat_file = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/RnaseqMapping/output/stat"
        satur_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/saturation"
        coverage_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/coverage"
        chrstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/chr_stat"
        bamstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191008/LncRna_tsg_35703/MapAssessment/output/distribution"
        library = "circle"
        method = "hisat"
        # self.add_mapping_stat(stat_file=stat_file, library=library, method=method)
        # self.add_rpkm_table(satur_dir, params='RSeQC-2.6.3', detail=True, library=library)
        self.add_coverage_table(coverage_dir, params='RSeQC-2.6.3', detail=True, library=library)
        # self.add_distribution_table(distribution=bamstat_dir, params='RSeQC-2.6.3', library=library)
        # self.add_chrom_distribution_table(distribution=chrstat_dir, params='RSeQC-2.6.3', library=library)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("whole_transcriptome.mapping")
        wf.test_api.export_mapping()


if __name__ == '__main__':
    unittest.main()
