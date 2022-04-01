# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last modified by shicaiping at 20180509
from __future__ import division
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
from bson.son import SON
import json
import glob
import re
import os
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from collections import OrderedDict, defaultdict
import unittest

class ProkRnaQc(ApiBase):
    def __init__(self, bind_object):
        super(ProkRnaQc, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    @report_check
    def add_samples_info(self, qc_stat_raw, qc_stat_clean, list_txt, qc_adapt=None, fq_type='pe', group=None):
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
        stat_file_raw = qc_stat_raw + "/fastq_stat.xls"
        stat_file_clean = qc_stat_clean + "/fastq_stat.xls"
        dup_file_raw = qc_stat_raw + "/dup.xls"
        dup_file_clean = qc_stat_clean + "/dup.xls"
        dup = ""
        dup_rate = {}
        adapter = False
        if not os.path.exists(stat_file_raw):
            raise Exception("%s文件不存在" % stat_file_raw)
        if not os.path.exists(stat_file_clean):
            raise Exception("%s文件不存在" % stat_file_clean)
        if not os.path.exists(list_txt):
            raise Exception("%s文件不存在" % list_txt)
        if qc_adapt is not None:
            adapt_rate = {}
            adapter = True
            with open(qc_adapt, "r") as a:
                for line in a:
                    line = line.split()
                    adapt_rate[line[0]] = line[1]
        if os.path.exists(dup_file_raw):
            with open(dup_file_raw, "r") as d:
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
            with open(dup_file_clean, "r") as d:
                col_num = len(d.readline().split())
                if col_num == 4:
                    for line in d:
                        line = line.split()
                        dup_rate[line[0]] += [float(line[1]), float(line[2]), float(line[3])]
                if col_num == 2:
                    for line in d:
                        line = line.split()
                        dup_rate[line[0]] = list(dup_rate[line[0]]).append(line[1])
        stat_info = OrderedDict()
        with open(stat_file_raw, "r") as f:
            data_list = []
            first_line = f.readline()
            if not re.match(r"#Sample_ID", first_line):
                raise Exception("%s文件类型不正确" % stat_file_raw)
            for line in f:
                line = line.strip().split()
                stat_info[line[0]] = line

        with open(stat_file_clean, "r") as f:
            data_list = []
            first_line = f.readline()
            if not re.match(r"#Sample_ID", first_line):
                raise Exception("%s文件类型不正确" % stat_file_raw)
            for line in f:
                line = line.strip().split()
                stat_info[line[0]] += line[1:]
        with open(list_txt, 'r') as list_r:
            sample2file = defaultdict(list)
            for line in list_r.readlines():
                line = line.strip().split('\t')
                sample2file[line[1]].append(line[0])
        for line in stat_info:
            read_length = int(int(stat_info[line][2])/int(stat_info[line][1]))
            if read_length == 151:
                read_length = 150
            read_length = str(read_length) + ' x ' + '2' if fq_type.lower() == 'pe' else '1'
            data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "old_name": stat_info[line][0],
                    "new_name": stat_info[line][0],
                    "raw_reads": int(stat_info[line][1]),
                    "raw_bases": int(stat_info[line][2]),
                    "raw_ns": int(stat_info[line][3]),
                    "raw_n": float(stat_info[line][4]),
                    "raw_a": float(stat_info[line][5]),
                    "raw_t": float(stat_info[line][6]),
                    "raw_c": float(stat_info[line][7]),
                    "raw_g": float(stat_info[line][8]),
                    "raw_n_rate": float(stat_info[line][9]),
                    "raw_error": float(stat_info[line][10]),
                    "raw_q20": float(stat_info[line][11]),
                    "raw_q30": float(stat_info[line][12]),
                    "raw_gc": float(stat_info[line][13]),
                    "clean_reads": int(stat_info[line][14]),
                    "clean_bases": int(stat_info[line][15]),
                    "clean_ns": int(stat_info[line][16]),
                    "clean_n": float(stat_info[line][17]),
                    "clean_a": float(stat_info[line][18]),
                    "clean_t": float(stat_info[line][19]),
                    "clean_c": float(stat_info[line][20]),
                    "clean_g": float(stat_info[line][21]),
                    "clean_n_rate": float(stat_info[line][22]),
                    "clean_error": float(stat_info[line][23]),
                    "clean_q20": float(stat_info[line][24]),
                    "clean_q30": float(stat_info[line][25]),
                    "gc_rate": float(stat_info[line][2]),
                    # "about_qc": about_qc,
                    "desc": "",
                    "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    "type": fq_type,   # 怎么得知待定
                    "raw_file": ','.join(sample2file[stat_info[line][0]]),
                    "read_length": read_length
                    }
            if line in dup_rate:
                if dup == "se":
                    data["reads1_dup_rate"] = dup_rate[line]
                if dup == "pe":
                    data["reads1_dup_rate"] = list(dup_rate[line][0]) + list(dup_rate[line[0]][3])
                    data["reads2_dup_rate"] = list(dup_rate[line][1]) + list(dup_rate[line[0]][4])
                    data["paired_dup_rate"] = list(dup_rate[line][1]) + list(dup_rate[line[0]][4])
            if adapter:
                if line in adapt_rate:
                    data["adapt_rate"] = adapt_rate[line]
                    data["about_qc"] = "after"
            data_list.append(data)
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['old_name']))
        try:
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入样品信息数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品信息数据成功")
            self.sample_ids = result.inserted_ids
        for i in result.inserted_ids:
            self.update_db_record('sg_specimen', i, status="end", main_id=i, fq_type=fq_type)
        sample_ids = [str(i) for i in result.inserted_ids]
        return sorted(sample_ids)

    @report_check
    def add_productive_name(self, samples=None, productive_table=None):
        collection = self.db["sg_specimen"]
        with open(productive_table, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                items = line.strip().split("\t")
                if len(items) >= 2 and items[0] in samples:
                    collection.update(
                        {"task_id": self.bind_object.sheet.id, "old_name": items[0]},
                        {"$set": {"productive_name": items[1]}}, upsert=False, multi=True)

    @report_check
    def add_gragh_info(self, quality_stat, about_qc="before"):
        """
        :param quality_stat: ouput_dir里一个叫qualityStat的文件夹，即~/output_dir/qualityStat
        :param about_qc:指控后的或是质控前的统计，质控前统计为before，指控后传after
        :return:
        """
        stat_files = sorted(glob.glob("{}/*".format(quality_stat)))
        data_list = []
        for sf in stat_files:
            sample_name = os.path.basename(sf).split(".")[0]
            self.bind_object.logger.info('%s,%s' % (sf, sample_name))
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
                        #"project_sn": self.bind_object.sheet.project_sn,
                        #"task_id": self.bind_object.sheet.id,
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
                        "N": int(line[16])/total_base * 100,
                        #'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    }
                    data_list.append(data)
        try:
            collection = self.db["sg_specimen_graphic"]
            result = collection.insert_many(data_list)
            # for i in result.inserted_ids:
            #     self.update_db_record('sg_specimen_graphic', i, status="end", main_id=i)
        except Exception as e:
            self.bind_object.set_error("导入样品画图数据信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品画图数据信息成功")

    @report_check
    def add_specimen_group(self, file):
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
                col.update({"task_id" : self.bind_object.sheet.id, "old_name": sample}, {"$set": {"group": key}}, upsert=True)
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
                tmp = line.strip().split()
                if len(tmp) > 1:
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
        except Exception as e:
            self.bind_object.set_error("导入样本对照组信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样本对照组信息成功")
            return com_id, con_list

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
            bam_path = os.path.join(dir_path, "Align/AlignBam", spname + ".bam")
            insert_data = {"bam_path": bam_path}
            col.update({"_id": sp_id}, {"$set": insert_data})

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
    def add_rfam_stat(self, stat_file, stat_file2=None, group=None):
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
            f = open(fs, "r")
            f.readline()
            items = f.readline().strip().split("\t")
            specimen_name = items[0]
            r_rna_ratio = items[1]
            data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": self.bind_object.sheet.id,
                "params": "rfam",
                "specimen_name": specimen_name,
                "status": "end",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "r_rna_ratio": r_rna_ratio
                }
            print stat_file2
            if stat_file2:
                with open(stat_file2 + '/{}.stat'.format(specimen_name), 'r') as f:
                    lines = f.readlines()
                    if 'bam_sort_core' not in lines[-2]:
                        data.update({
                            "r_rna_ratio2": lines[-2].split(" ")[0].rstrip("%")
                        })
                    else:
                        data.update({
                            "r_rna_ratio2": lines[-3].split(" ")[0].rstrip("%")
                        })
            data_list.append(data)
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_specimen_rfam"]
            result = collection.insert_many(data_list)
            for i in result.inserted_ids:
                self.update_db_record('sg_specimen_rfam', i, status="end", main_id=i)

        except Exception as e:
            print("导入rfam注释信息出错:%s" % e)
        else:
            print("导入rfam注释信息成功")

    @report_check
    def add_mapping_stat(self, stat_file, type="genome"):
        """
        :param stat_file: 比对结果统计文件，在mapassessment module的ouput_dir中，即~/MapAssessment/output/bam_stat.xls
        :param type: 比对到基因组的传genome，比对到基因的传gene
        :return:
        """
        with open(stat_file, "r") as f:
            data_list = []
            f.readline()
            for line in f:
                line = line.strip().split()
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "type": type,
                    "specimen_name": line[0],
                    "total_reads": line[1],
                    "mapping_reads": line[2] + "(" + str(float("%0.4f" % (float(line[2])/float(line[1])))*100) + "%" + ")",
                    "multiple_mapped": line[3] + "(" + str(float("%0.4f" % (float(line[3])/float(line[1])))*100) + "%" + ")",
                    "uniq_mapped": line[4] + "(" + str(float("%0.4f" % (float(line[4])/float(line[1])))*100) + "%" + ")",
                    "map_to_up": line[5] + "(" + str(float("%0.4f" % (float(line[5])/float(line[1])))*100) + "%" + ")",
                    "map_to_down": line[6] + "(" + str(float("%0.4f" % (float(line[6])/float(line[1])))*100) + "%" + ")",
                    "non_splice_reads": line[7] + "(" + str(float("%0.4f" % (float(line[7])/float(line[1])))*100) + "%" + ")",
                    "splice_reads": line[8] + "(" + str(float("%0.4f" % (float(line[8])/float(line[1])))*100) + "%" + ")",
                    "mapping_rate": str(float("%0.4f" % (float(line[2])/float(line[1])))*100) + "%",
                    "multiple_rate": str(float("%0.4f" % (float(line[3])/float(line[1])))*100) + "%",
                    "uniq_rate": str(float("%0.4f" % (float(line[4])/float(line[1])))*100) + "%",
                    'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                }
                data_list.append(data)
        try:
            collection = self.db["sg_specimen_mapping"]
            result = collection.insert_many(data_list)
            for i in result.inserted_ids:
                self.update_db_record('sg_specimen_mapping', i, status="end", main_id=i)
        except Exception as e:
            self.bind_object.set_error("导入比对结果统计信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入比对结果统计信息成功")

    @report_check
    def add_rpkm_table(self, file_path, name=None, params='RSeQC-2.6.3', detail=True):
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
            "curve_category": ['5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85', '90', '95', '100'],
            "curve_specimen": {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_assessment_saturation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_rpkm_curve(file_path, inserted_id)
            self.update_db_record('sg_assessment_saturation', inserted_id, status="end", main_id=inserted_id)
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
                    new_line = [ float(x) for x in line ]
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
            collection = self.db["sg_assessment_saturation_curve"]
            collection.insert_many(curve_data)
        except Exception as e:
            self.bind_object.set_error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_satur_count(self, count_r_file):
        categaries = {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"}
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
    def add_coverage_table(self, coverage, name=None, params='RSeQC-2.6.3', detail=True):
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
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_assessment_coverage"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_coverage_detail(coverage, inserted_id)
            self.update_db_record('sg_assessment_coverage', inserted_id, status="end", main_id=inserted_id)
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
            sample_name = os.path.basename(cf).split(".")[0][9:]
            # self.bind_object.set_error sample_name

            with open(cf, "r") as f:
                percent = f.readline().strip().split()
                value = f.next().strip().split()
                plot_value = {}
                for i in range(100):
                    plot_value[percent[i+1]] = value[i+1]
                data = {
                    "coverage_id": coverage_id,
                    "specimen_name": sample_name,
                    "plot_value": plot_value
                }
                data_list.append(data)
        try:
            collection = self.db["sg_assessment_coverage_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_distribution_table(self, distribution=None, name=None, params='RSeQC-2.6.3'):
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
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_assessment_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_distribution_detail(distribution, inserted_id)
            self.update_db_record('sg_assessment_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_distribution_detail(self, distribution, distribution_id):
        data_list = []
        stat_files = sorted(glob.glob("{}/*reads_distribution.txt".format(distribution)))
        for fls in stat_files:
            sample_name = os.path.basename(fls).split(".")[0]
            with open(fls, "r") as f:
                intergenic = 0
                total = 0
                cds = 0
                utr3 = 0
                utr5 = 0
                introns = 0
                for line in f:
                    if line.startswith("Total Reads"):
                        pass
                    elif line.startswith("Total Tags"):
                        total = float(line.strip().split()[2])
                    elif line.startswith("Total Assigned Tags"):
                        pass
                    elif line.startswith("====="):
                        pass
                    elif line.startswith("Group"):
                        pass
                    elif line.startswith("CDS_Exons"):
                        cds = float(line.strip().split()[2])
                    elif line.startswith("5'UTR_Exons"):
                        utr5 = float(line.strip().split()[2])
                    elif line.startswith("3'UTR_Exons"):
                        utr3 = float(line.strip().split()[2])
                    elif line.startswith("introns"):
                        introns = float(line.strip().split()[2])
                    else:
                        intergenic += float(line.strip().split()[2])
                if total == 0:
                    total = cds + utr5 + utr3 + introns + intergenic

                data = {
                    "distribution_id": distribution_id,
                    "specimen_name": sample_name,
                    "cds": str(cds) + "(" + str(float("%0.4f" % (cds / total)) * 100) + ")",
                    "utr3": str(utr3) + "(" + str(float("%0.4f" % (utr3 / total)) * 100) + ")",
                    "utr5": str(utr3) + "(" + str(float("%0.4f" % (utr5 / total)) * 100) + ")",
                    "intergenic": str(intergenic) + "(" + str(float("%0.4f" % (intergenic / total)) * 100) + ")",
                    "introns": str(introns) + "(" + str(float("%0.4f" % (introns / total)) * 100) + ")",
                }
                print data
                data_list.append(data)
                self.db["sg_specimen_mapping"].update({"task_id": self.bind_object.sheet.id, "specimen_name": sample_name}, {"$set": {"cds": str(int(float(data["cds"].split('(')[0]))),"cds_ratio": data["cds"].split('(')[1].strip(')')}})
        try:
            collection = self.db["sg_assessment_distribution_detail"]
            result = collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入reads distribution详情表出错:%s" % e)
        else:
            self.bind_object.logger.info("导入reads distribution详情表成功")

    @report_check
    def add_chrom_distribution_table(self, distribution=None, name=None, params='RSeQC-2.6.3'):
        """
        :param distribution: 文件夹，即~/MapAssessment/output/chr_stat,不传的时候只导主表
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
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_assessment_chrom_distribution"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if distribution:
            self.add_chorm_distribution_detail(distribution, inserted_id)
            self.update_db_record('sg_assessment_chrom_distribution', inserted_id, status="end", main_id=inserted_id)
        return inserted_id

    @report_check
    def add_chorm_distribution_detail(self, chorm_stat, chorm_stat_id):
        stat_files = sorted(glob.glob("{}/*.bam_chr_stat.xls".format(chorm_stat)))
        data_list = []
        distribution = set()
        specimen_names = []
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
            collection = self.db["sg_assessment_chrom_distribution_detail"]
            result = collection.insert_many(data_list)
            main_collection = self.db["sg_assessment_chrom_distribution"]
            main_collection.update({"_id": ObjectId(chorm_stat_id)}, {"$set": {"distribution": list(distribution), "specimen_name": specimen_names}})
        except Exception as e:
            self.bind_object.logger.info("导入sg_assessment_chrom_distribution_detail出错:%s" % e)
        else:
            self.bind_object.logger.info("导入sg_assessment_chrom_distribution_detail成功")

    @report_check
    def add_tophat_mapping_stat(self, stat_file):
        files = sorted(glob.glob("{}/*".format(stat_file)))
        data_list = []
        for fs in files:
            specimen_name = os.path.basename(fs).split(".")[0]
            print(specimen_name)
            f = open(fs, "r")
            data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_name": specimen_name,
                    "params": "tophat",
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            total_reads = 0
            map_reads = 0
            multiple = 0
            for line in f:
                # print map_reads
                if re.match(r"Left", line):
                    total_reads += int(f.next().split()[-1])
                    map_reads += int(f.next().split()[2])
                    multiple += int(f.next().split()[2])
                if re.match(r"Right", line):
                    total_reads += int(f.next().split()[-1])
                    map_reads += int(f.next().split()[2])
                    multiple += int(f.next().split()[2])
            print(total_reads, map_reads, multiple, total_reads - multiple)
            data["total_reads"] = total_reads
            # print "(" + str(float("%0.4f" % ( map_reads/total_reads)) * 100) + "%" + ")"
            data["mapping_reads"] = str(map_reads) + "(" + str(float("%0.4f" % ( map_reads/total_reads)) * 100) + "%" + ")"
            data["multiple_mapped"] = str(multiple) + "(" + str(float("%0.4f" % ( multiple/total_reads)) * 100) + "%" + ")"
            data["uniq_mapped"] = str(total_reads - multiple) + "(" + str(float("%0.4f" % ((map_reads - multiple)/total_reads)) * 100) + "%" + ")"
            # print data
            data_list.append(data)
            f.close()
        try:
            collection = self.db["sg_specimen_mapping"]
            result = collection.insert_many(data_list)
            for i in result.inserted_ids:
                self.update_db_record('sg_specimen_mapping', i, status="end", main_id=i)
        except Exception as e:
            print("导入比对结果统计信息出错:%s" % e)
        else:
            print("导入比对结果统计信息成功")

    @report_check
    def add_bowtie2_mapping_stat(self, stat_file, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        files = sorted(glob.glob("{}/*".format(stat_file)))
        # print files
        data_list = []
        for fs in files:
            specimen_name = os.path.basename(fs).split(".")[0]
            # print specimen_name
            values = []
            f = open(fs, "r")
            for line in f:
                # print line
                if re.match(r" ", line):
                    line = line.split()
                    values.append(line[0])
            if len(values) == 13:
                total_reads = int(values[0])*2
                unmap_reads = int(values[-3])
                map_reads = int(total_reads) - int(unmap_reads)
                map_reads = str(map_reads) + "(" + str(float("%0.4f" % (map_reads/int(total_reads))) * 100) + ")"
                uniq_mapped = int(values[2])*2 + int(values[6])*2 + int(values[-2])
                uniq_mapped = str(uniq_mapped) + "(" + str(float("%0.4f" % (uniq_mapped/int(total_reads))) * 100) + ")"
                multi_mapped = int(values[3])*2 + int(values[-1])
                multi_mapped = str(multi_mapped) + "(" + str(float("%0.4f" % (multi_mapped/int(total_reads))) * 100) + ")"
                print(specimen_name)
                print(total_reads, map_reads, uniq_mapped, multi_mapped)
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "params": "hisat",
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads.split('(')[0],
                    "mapping_reads_ratio": map_reads.split('(')[1].strip(')'),
                    "unmapping_reads": unmap_reads,
                    "unmapping_reads_ratio": str(float("%0.4f" % (unmap_reads/int(total_reads))) * 100),
                    "multiple_mapped": multi_mapped.split('(')[0],
                    "multiple_mapped_ratio": multi_mapped.split('(')[1].strip(')'),
                    "uniq_mapped": uniq_mapped.split('(')[0],
                    "uniq_mapped_ratio": uniq_mapped.split('(')[1].strip(')'),
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                data_list.append(data)
            elif len(values) == 4:
                total_reads = int(values[0])
                unmap_reads = int(values[1])
                map_reads = int(total_reads) - int(unmap_reads)
                uniq_mapped = int(values[2])
                multi_mapped = int(values[3])
                map_reads = str(map_reads) + "(" + str(float("%0.4f" % (map_reads/int(total_reads))) * 100) + ")"
                uniq_mapped = str(uniq_mapped) + "(" + str(float("%0.4f" % (uniq_mapped/int(total_reads))) * 100) + ")"
                multi_mapped = str(multi_mapped) + "(" + str(float("%0.4f" % (multi_mapped/int(total_reads))) * 100) + ")"
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "params": "hisat",
                    "specimen_name": specimen_name,
                    "total_reads": total_reads,
                    "mapping_reads": map_reads.split('(')[0],
                    "mapping_reads_ratio": map_reads.split('(')[1].strip(')'),
                    "multiple_mapped": multi_mapped.split('(')[0],
                    "multiple_mapped_ratio": multi_mapped.split('(')[1].strip(')'),
                    "uniq_mapped": uniq_mapped.split('(')[0],
                    "uniq_mapped_ratio": uniq_mapped.split('(')[1].strip(')'),
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                data_list.append(data)
            f.close()
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['specimen_name']))
        try:
            collection = self.db["sg_specimen_mapping"]
            result = collection.insert_many(data_list)
            for i in result.inserted_ids:
                self.update_db_record('sg_specimen_mapping', i, status="end", main_id=i)
        except Exception as e:
            print("导入比对结果统计信息出错:%s" % e)
        else:
            print("导入比对结果统计信息成功")

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "prok_rna_srna__",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "prok_rna_srna",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "prok_rna.prokrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ProkrnaTestApiWorkflow(wsheet)

        stat_file_raw = '/mnt/ilustre/users/sanger-dev/workspace/20180807/Single_HiseqReadsStat_4159/HiseqReadsStat/output'
        stat_file_clean = '/mnt/ilustre/users/sanger-dev/workspace/20180806/Single_HiseqReadsStat_1333/HiseqReadsStat/output'
        bowtie2_stat = '/mnt/ilustre/users/sanger-dev/workspace/20180809/Single_RnaseqMapping_3447/RnaseqMapping/output/stat'
        file_path = '/mnt/ilustre/users/sanger-dev/workspace/20180809/Single_MapAssessment_9826/MapAssessment/output/saturation'
        coverage = '/mnt/ilustre/users/sanger-dev/workspace/20180809/Single_MapAssessment_9826/MapAssessment/output/coverage'
        distribution = '/mnt/ilustre/users/sanger-dev/workspace/20180809/Single_MapAssessment_9826/MapAssessment/output/distribution'
        chrom_distribution = '/mnt/ilustre/users/sanger-dev/workspace/20180809/Single_MapAssessment_9826/MapAssessment/output/chr_stat'
        list_txt = '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/data/list.txt'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("prok_rna.prok_rna_qc")
        params = {
            "name": 'try_all',
            "database": 'all',
            "software": 'RNAplex,RNAhybrid',
        }

        wf.test_api.add_samples_info(stat_file_raw, stat_file_clean, list_txt)
        wf.test_api.add_gragh_info(stat_file_raw + '/qualityStat', about_qc="before")
        wf.test_api.add_gragh_info(stat_file_clean + '/qualityStat', about_qc="after")
        wf.test_api.add_bowtie2_mapping_stat(bowtie2_stat)
        wf.test_api.add_rpkm_table(file_path)
        wf.test_api.add_coverage_table(coverage)
        wf.test_api.add_distribution_table(distribution)
        wf.test_api.add_chrom_distribution_table(chrom_distribution)

if __name__ == '__main__':
    unittest.main()
