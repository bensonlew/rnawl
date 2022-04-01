# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# last modified by zhangyitong at 20200903

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

from mbio.api.database.medical_transcriptome.api_base import ApiBase


class Mapping(ApiBase):
    def __init__(self, bind_object):
        super(Mapping, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'

    @report_check
    def add_bam_path(self, dir_path):
        """
        将bam文件的路径插入sg_specimen表中，供可变剪切使用
        :param dir_path:传入的rnaseq_mapping的output_dir
        :return:
        """
        spname_spid = self.get_spname_spid()
        col = self.db["sg_specimen"]
        for spname in spname_spid:
            sp_id = spname_spid[spname]
            bam_path = dir_path + "/Align/AlignBam/" + spname + ".bam"
            insert_data = {"bam_path": bam_path}
            col.update({"_id": sp_id}, {"$set": insert_data})

    @report_check
    def get_spname_spid(self):
        if not self.sample_ids:
            self.bind_object.set_error("样本id列表为空，请先调用add_samples_info产生specimen的id", code="53702006")
        collection = self.db["sg_specimen"]
        spname_spid = {}
        for id_ in self.sample_ids:
            results = collection.find_one({"_id": id_})
            spname_spid[results['new_name']] = id_
        return spname_spid

    @report_check
    def add_mapping_stat(self, stat_file, method, group, sample_list_file=None):
        name = "rna_" + method + "_mapping"
        bam_path = os.path.dirname(stat_file + "/bam")
        params = json.dumps(
            {'task_id': self.bind_object.sheet.id, 'submit_location': 'mapping', 'method': method},
            sort_keys=True)
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "params": params,
            "method": method,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "start",
            "desc": '比对分析结果主表',
            "name": name,
            'version': 'v1',
            'bam_path': bam_path
        }
        try:
            collection = self.db['sg_mapping']
            mapping_id = collection.insert_one(data).inserted_id
        except Exception, e:
            self.bind_object.set_error("导入比对信息主表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入比对信息主表成功")
            collection.update({"_id": mapping_id, "task_id": self.bind_object.sheet.id},
                              {"$set": {"main_id": mapping_id, "status": "end"}}, upsert=True)

            if method.lower() == "hisat":
                self.add_hisat_mapping_stat(stat_file, mapping_id, group)
            elif method.lower() == "tophat":
                self.add_tophat_mapping_stat(stat_file, mapping_id, group)
            elif method.lower() == "star":
                self.add_star_mapping_stat(stat_file, mapping_id, group)

    @report_check
    def add_rpkm_table(self, file_path, group, name=None, params='RSeQC-2.6.3', detail=True):
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
            "status": "end",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "curve_category": ['5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75',
                               '80', '85', '90', '95', '100'],
            "curve_specimen": {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)",
                               "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'version': 'v1'
        }
        collection = self.db["sg_assessment_saturation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_rpkm_curve(file_path, group, inserted_id)
            self.update_db_record('sg_assessment_saturation', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_rpkm_curve(self, rpkm_file, group, rpkm_id=None):
        """
        :param rpkm_file:文件夹，即~/MapAssessment/output/satur
        :param rpkm_id:主表id，必传
        :return:
        """
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
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
        curve_data.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_assessment_saturation_detail"]
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
    def add_coverage_table(self, coverage, group, name=None, params='RSeQC-2.6.3', detail=True):
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
            'version': 'v1'
        }
        collection = self.db["sg_assessment_coverage"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_coverage_detail(coverage, group, inserted_id)
            self.update_db_record('sg_assessment_coverage', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_coverage_detail(self, coverage, group, coverage_id=None):
        """
        :param coverage: 文件夹，即~/MapAssessment/output/coverage
        :param coverage_id: 主表ID，必传
        :return:
        """
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
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
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_assessment_coverage_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_distribution_table(self, distribution=None, group=None, name=None, params='RSeQC-2.6.3'):
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
            'version': 'v1'
        }
        collection = self.db["sg_assessment_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_distribution_detail(distribution, group, inserted_id)
            self.update_db_record('sg_assessment_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_distribution_detail(self, distribution, group, distribution_id):
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
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_assessment_distribution_detail"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入reads distribution详情表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入reads distribution详情表成功")

    @report_check
    def add_chrom_distribution_table(self, distribution=None, group=None, name=None, params='RSeQC-2.6.3'):
        """
        :param distribution: 文件夹，即~/MapAssessment/output/chr_stat,不传的时候只导主表
        :param name: 主表名称。可不传
        :param params:参数，可不传
        :return:
        """
        # 改变染色体分布导表方式，原因在于需要对染色体丰度进行排序
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "distribution_origin",
            "status": "end",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'version': 'v1'
        }
        collection = self.db["sg_assessment_chrom_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_chrom_distribution_detail(distribution, group, inserted_id)
            self.update_db_record('sg_assessment_chrom_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_chrom_distribution_detail(self, chorm_stat, group, chorm_stat_id):
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
        data_list = []
        distribution = set()
        specimen_names = []
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
                data_list.append(data)
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_assessment_chrom_distribution_detail"]
            result = collection.insert_many(data_list)
            main_collection = self.db["sg_assessment_chrom_distribution"]
            main_collection.update({"_id": ObjectId(chorm_stat_id)},
                                   {"$set": {"distribution": list(distribution), "specimen_names": sample_list}})
        except Exception, e:
            self.bind_object.logger.info("导入assessment_chrom_distribution_detail出错:%s" % e)
        else:
            self.bind_object.logger.info("导入assessment_chrom_distribution_detail成功")

    @report_check
    def add_tophat_mapping_stat(self, stat_file, main_mapping_id, group):
        files = sorted(glob.glob("{}/*".format(stat_file)))
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
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
            data["mapping_reads"] = str(map_reads)
            data["mapping_reads_percent"] = str(float("%0.4f" % (map_reads / total_reads)) * 100) + "%"
            data["multiple_mapped"] = str(multiple)
            data["multiple_mapped_percent"] = str(float("%0.4f" % (multiple / total_reads)) * 100) + "%"
            data["uniq_mapped"] = str(total_reads - multiple)
            data["uniq_mapped_percent"] = str(float("%0.4f" % ((map_reads - multiple) / total_reads)) * 100) + "%"
            # print data
            data_list.append(data)
            f.close()
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    @report_check
    def add_hisat_mapping_stat(self, stat_file, main_mapping_id, group):
        files = sorted(glob.glob("{}/*".format(stat_file)))
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
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
                map_reads = str(map_reads)
                map_reads_percent = str(float("%0.4f" % (int(map_reads) / int(total_reads))) * 100) + "%"
                uniq_mapped = int(values[2]) * 2 + int(values[6]) * 2 + int(values[-2])
                uniq_mapped = str(uniq_mapped)
                uniq_mapped_percent = str(float("%0.4f" % (int(uniq_mapped) / int(total_reads))) * 100) + "%"
                multi_mapped = int(values[3]) * 2 + int(values[-1])
                multi_mapped = str(multi_mapped)
                multi_mapped_percent = str(float("%0.4f" % (int(multi_mapped) / int(total_reads))) * 100) + "%"
                print specimen_name
                print total_reads, map_reads, uniq_mapped, multi_mapped
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads,
                    "mapping_reads_percent": map_reads_percent,
                    "multiple_mapped": multi_mapped,
                    "multiple_mapped_percent": multi_mapped_percent,
                    "uniq_mapped": uniq_mapped,
                    "uniq_mapped_percent": uniq_mapped_percent,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            elif len(values) == 4:
                total_reads = int(values[0])
                unmap_reads = int(values[1])
                map_reads = int(total_reads) - int(unmap_reads)
                uniq_mapped = int(values[2])
                multi_mapped = int(values[3])
                map_reads = str(map_reads)
                map_reads_percent = str(float("%0.4f" % (int(map_reads) / int(total_reads))) * 100) + "%"
                uniq_mapped = str(uniq_mapped)
                uniq_mapped_percent = str(float("%0.4f" % (int(uniq_mapped) / int(total_reads))) * 100) + "%"
                multi_mapped = str(multi_mapped)
                multi_mapped_percent = str(float("%0.4f" % (int(multi_mapped) / int(total_reads))) * 100) + "%"
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads,
                    "mapping_reads_percent": map_reads_percent,
                    "multiple_mapped": multi_mapped,
                    "multiple_mapped_percent": multi_mapped_percent,
                    "uniq_mapped": uniq_mapped,
                    "uniq_mapped_percent": uniq_mapped_percent,
                    "mapping_id": main_mapping_id,
                }
                data_list.append(data)
            f.close()
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    @report_check
    def add_star_mapping_stat(self, stat_file, main_mapping_id, group):
        files = sorted(glob.glob("{}/*".format(stat_file)))
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
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
                map_reads_rate = float("%0.4f" % (int(map_reads) / total_reads)) * 100
                uniq_mapped_rate = float("%0.4f" % (int(uniq_mapped) / total_reads)) * 100
                multi_mapped_rate = float("%0.4f" % (int(multi_mapped) / total_reads)) * 100
                multi_mapped_f = str(multi_mapped)
                map_reads_f = str(map_reads)
                uniq_mapped_f = str(uniq_mapped)
                print specimen_name
                print total_reads, map_reads, uniq_mapped, multi_mapped
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads_f,
                    "mapping_reads_percent": map_reads_rate,
                    "multiple_mapped": multi_mapped_f,
                    "multiple_mapped_percent": multi_mapped_rate,
                    "uniq_mapped": uniq_mapped_f,
                    "uniq_mapped_percent": uniq_mapped_rate,
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
                map_reads_rate = float("%0.4f" % (int(map_reads) / total_reads)) * 100
                uniq_mapped_rate = float("%0.4f" % (int(uniq_mapped) / total_reads)) * 100
                multi_mapped_rate = float("%0.4f" % (int(multi_mapped) / total_reads)) * 100
                multi_mapped_f = str(multi_mapped)
                map_reads_f = str(map_reads)
                uniq_mapped_f = str(uniq_mapped)
                data = {
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads_f,
                    "mapping_reads_percent": map_reads_rate,
                    "multiple_mapped": multi_mapped_f,
                    "multiple_mapped_percent": multi_mapped_rate,
                    "uniq_mapped": uniq_mapped_f,
                    "uniq_mapped_percent": uniq_mapped_rate,
                    "mapping_id": main_mapping_id
                }
                data_list.append(data)
            f.close()
        data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_mapping_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入比对结果详情表出错:%s" % e)
        else:
            print("导入比对结果详情表成功")

    def export_mapping(self):
        stat_file = "/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/RnaseqMapping/output/stat"
        satur_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/MapAssessment/output/saturation"
        coverage_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/MapAssessment/output/coverage"
        chrstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/MapAssessment/output/chr_stat"
        bamstat_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/MapAssessment/output/distribution"
        method = "hisat"
        self.add_mapping_stat(stat_file=stat_file, method=method)
        self.add_rpkm_table(satur_dir, params='RSeQC-2.6.3', detail=True)
        self.add_coverage_table(coverage_dir, params='RSeQC-2.6.3', detail=True)
        self.add_distribution_table(distribution=bamstat_dir, params='RSeQC-2.6.3')
        self.add_chrom_distribution_table(distribution=chrstat_dir, params='RSeQC-2.6.3')


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "medical_transcriptome",
            "project_sn": "medical_transcriptome",
            "type": "workflow",
            "name": "medical_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.mapping")
        wf.test_api.export_mapping()


if __name__ == '__main__':
    unittest.main()
