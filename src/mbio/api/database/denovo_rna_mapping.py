# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from bson.objectid import ObjectId
from cStringIO import StringIO
import bson.binary
import datetime
import pandas
import numpy
import json
import glob
import re
import os


class DenovoRnaMapping(Base):
    def __init__(self, bind_object):
        super(DenovoRnaMapping, self).__init__(bind_object)
        self._project_type = 'ref_rna'

    @report_check
    def add_mapping_stat(self, stat_file):
        with open(stat_file, "r") as f:
            data_list = []
            f.readline()
            for line in f:
                line = line.strip().split()
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_name": line[0],
                    "mapping_reads": line[1],
                    "mapping_rate": str(float("%0.4f" % (int(line[2])/int(line[1])))*100) + "%",
                    "multiple_mapped": line[3],
                    "multiple_rate": str(float("%0.4f" % (int(line[3])/int(line[1])))*100) + "%",
                    "uniq_mapped": line[4],
                    "uniq_rate": str(float("%0.4f" % (int(line[4])/int(line[1])))*100) + "%"
                }
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_specimen_mapping"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入比对结果统计信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入比对结果统计信息成功")

    @report_check
    def add_rpkm_table(self, file_path, name=None, params=None, detail=True):
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "saturation_origin",
            "status": "start",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "curve_specimen": {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_rpkm"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_rpkm_detail(file_path, inserted_id)
            self.add_rpkm_box(file_path, inserted_id)
            self.add_rpkm_curve(file_path, inserted_id)
        return inserted_id

    @report_check
    def add_rpkm_detail(self, rpkm_file, rpkm_id=None):
        rpkm_id = ObjectId(rpkm_id)
        rpkm_tables = glob.glob("{}/*eRPKM.xls".format(rpkm_file))
        rpkm_detail = []
        for rt in rpkm_tables:
            sample_name = os.path.basename(rt).split(".")[0][6:]
            with open(rt, "r") as f:
                column_name = f.readline().strip().split()[6:]
                for line in f:
                    line = line.strip().split()
                    col_num = 6
                    data = {
                        "rpkm_id": rpkm_id,
                        "specimen_name": sample_name,
                        # "specimen_id" : get_sample_id(sample_name, "after"),
                        "transcript_id": line[0]
                    }
                    for col in column_name:
                        data["{}".format(col)] = line[col_num]
                        col_num += 1
                    rpkm_detail.append(data)
        try:
            collection = self.db["sg_denovo_rpkm_detail"]
            collection.insert_many(rpkm_detail)
        except Exception, e:
            self.bind_object.logger.error("导入rpkm detail出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm detail成功")

    def add_rpkm_box(self, rpkm_file, rpkm_id=None):
        rpkm_plot = glob.glob("{}/*saturation.r".format(rpkm_file))
        rpkm_box = []
        for rp in rpkm_plot:
            sample_name = os.path.basename(rp).split(".")[0][6:]
            sam_dict = {}
            with open(rp, "r") as f:
                for line in f:
                    if re.match(r"S", line):
                        split_line = line.strip().split("=c(")
                        sampling_percent = split_line[0][1:]
                        box_data = split_line[1][:-1].split(",")
                        box_data = map(float, box_data)
                        data = pandas.DataFrame({"name": box_data})
                        result = data.boxplot(return_type='dict')
                        qualities = {
                            "q1": result['whiskers'][0].get_data()[1][0],
                            "mean": numpy.mean(box_data, axis=0),
                            "q3": result['whiskers'][1].get_data()[1][0],
                            "max": max(box_data),
                            "min": min(box_data)
                            }
                        if sampling_percent in sam_dict:
                            sam_dict[sampling_percent].append(qualities)
                        else:
                            sam_dict[sampling_percent] = [qualities]
            for sam in sam_dict:
                data = {
                    "rpkm_id": rpkm_id,
                    "specimen_name": sample_name,
                    "sampling": sam,
                    "Q1": sam_dict[sam][0],
                    "Q2": sam_dict[sam][1],
                    "Q3": sam_dict[sam][2],
                    "Q4": sam_dict[sam][3]
                }
                rpkm_box.append(data)
        try:
            collection = self.db["sg_denovo_rpkm_box"]
            collection.insert_many(rpkm_box)
        except Exception, e:
            self.bind_object.logger.error("导入rpkm箱线图数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm箱线图数据")

    @report_check
    def add_rpkm_curve(self, rpkm_file, rpkm_id=None):
        rpkm_id = ObjectId(rpkm_id)
        curve_files = glob.glob("{}/*cluster_percent.xls".format(rpkm_file))
        rpkm_pdf = glob.glob("{}/*.pdf".format(rpkm_file))
        erpkm = glob.glob("{}/*.eRPKM.xls".format(rpkm_file))
        # curve_category = []
        with open(erpkm[0], "r") as f:
            category_line = f.readline().strip().split("\t")[6:]
            curve_category = [i[:-1] for i in category_line]
            # print(category_line)
        curve_data = []
        for cf in curve_files:
            sample_name = os.path.basename(cf).split(".")[0][6:]
            with open(cf, "r") as f:
                line_list = []
                for line in f:
                    line = line.strip().split()
                    line.pop(0)
                    line_list.append(line)
                data = {
                    "rpkm_id": rpkm_id,
                    "specimen_name": sample_name,
                    "column1": line_list[0],
                    "column2": line_list[1],
                    "column3": line_list[2],
                    "column4": line_list[3],
                    "column5": line_list[4],
                    "column6": line_list[5]
                }
                # self.bind_object.logger.error data
                curve_data.append(data)

        for rp in rpkm_pdf:
            sample_name = os.path.basename(rp).split(".")[0][6:]
            # print sample_name
            with open(rp, 'rb') as s:
                box_id = StringIO(s.read())
                box_id = bson.binary.Binary(box_id.getvalue())
                for cd in curve_data:
                    if cd["specimen_name"] == sample_name:
                        cd["box_id"] = box_id
        # print curve_data
        try:
            collection = self.db["sg_denovo_rpkm_curve"]
            main_collection = self.db["sg_denovo_rpkm"]
            collection.insert_many(curve_data)
            main_collection.update({"_id": ObjectId(rpkm_id)}, {"$set": {"curve_category": curve_category, "status": "end"}})
        except Exception, e:
            self.bind_object.logger.error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_coverage_table(self, coverage, name=None, params=None, detail=True):
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "coverage_origin",
            "status": "start",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_coverage"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_coverage_detail(coverage, inserted_id)
        return inserted_id

    @report_check
    def add_coverage_detail(self, coverage, coverage_id=None):
        coverage_files = glob.glob("{}/*geneBodyCoverage.txt".format(coverage))
        data_list = []
        for cf in coverage_files:
            sample_name = os.path.basename(cf).split(".")[0][9:]
            # self.bind_object.logger.error sample_name
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
            collection = self.db["sg_denovo_coverage_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入rpkm曲线数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入rpkm曲线数据成功")

    @report_check
    def add_duplication_table(self, dup, name=None, params=None, detail=True):
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": name if name else "duplication_origin",
            "status": "start",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_duplicate"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            self.add_duplication_detail(dup, inserted_id)
        return inserted_id

    @report_check
    def add_duplication_detail(self, dup, dup_id=None):
        dup_seqs = glob.glob("{}/*seq.DupRate.xls".format(dup))
        dup_poss = glob.glob("{}/*pos.DupRate.xls".format(dup))
        data_list = []
        for dp in dup_seqs:
            sample_name = os.path.basename(dp).split(".")[0][4:]
            with open(dp, "r") as f:
                f.readline()
                map_point = {}
                for line in f:
                    line = line.strip().split()
                    map_point[line[0]] = line[1]
                data = {
                    "dup_id": dup_id,
                    "specimen_name": sample_name,
                    "map_point": map_point
                }
                data_list.append(data)

        for ds in dup_poss:
            sample_name = os.path.basename(ds).split(".")[0][4:]
            with open(ds, "r") as f:
                f.readline()
                seq_point = {}
                for line in f:
                    line = line.strip().split()
                    seq_point[line[0]] = line[1]
            for dl in data_list:
                if sample_name == dl["specimen_name"]:
                    dl["seq_point"] = seq_point
        try:
            collection = self.db["sg_denovo_duplicate_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入冗余分析数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入冗余分析数据成功")

    @report_check
    def add_correlation_table(self, correlation, name=None, params=None, express_id=None, detail=True, seq_type=None):
        correlation_tree = glob.glob("{}/*.tre".format(correlation))
        with open(correlation_tree[0], "r") as t:
            correlation_tree = t.readline().strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', correlation_tree)
            tree_list = [i[1] for i in raw_samp]
        # if not params:
        #     params = dict()
        #     params['express_id'] = express_id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "type": seq_type,
            "name": name if name else "correlation_origin_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "status": "start",
            "desc": "",
            "correlation_tree": correlation_tree,
            "tree_list": tree_list,
            # "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_correlation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if detail:
            # self.add_correlation_detail(correlation, inserted_id)
            pca_file = os.path.join(correlation, 'pca_importance.xls')
            pca_rotation = os.path.join(correlation, 'pca_rotation.xls')
            site_file = os.path.join(correlation, 'pca_sites.xls')
            self.add_correlation_detail(collection=correlation, correlation_id=inserted_id, updata_tree=True)
            if os.path.exists(pca_file):
                self.add_pca(pca_file=pca_file, correlation_id=inserted_id)
                self.add_pca_rotation(input_file=pca_rotation, db_name='sg_denovo_correlation_pca_rotation', correlation_id=inserted_id)
                self.add_pca_rotation(input_file=site_file, db_name='sg_denovo_correlation_pca_sites', correlation_id=inserted_id)
        return inserted_id

    @report_check
    def add_correlation_detail(self, collection, correlation_id=None, updata_tree=False):
        correlation_id = ObjectId(correlation_id)
        correlation_matrix = collection + "/correlation_matrix.xls"
        data_list = []
        with open(correlation_matrix, "r") as m:
            samples = m.readline().strip().split()
            for line in m:
                data = {
                    "correlation_id": correlation_id
                }
                line = line.strip().split()
                for i, s in enumerate(samples):
                    data["specimen_name"] = line[0]
                    data[s] = line[i+1]
                # self.bind_object.logger.error data
                data_list.append(data)
        if updata_tree:
            col_tree = collection + "/corr_col.tre"
            row_tree = collection + "/corr_row.tre"
            f = open(col_tree, "r")
            r = open(row_tree, "r")
            tree_col = f.readline().strip()
            tree_row = r.readline().strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', tree_col)
            tree_list = [i[1] for i in raw_samp]
            f.close()
            r.close()
            collection_first = self.db['sg_denovo_correlation']
            collection_first.update({"_id": ObjectId(correlation_id)}, {"$set": {"correlation_tree": tree_row, "row_tree": tree_row, "col_tree": tree_col, "tree_list": tree_list}})
        try:
            collection = self.db["sg_denovo_correlation_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入相关系数分析数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入相关系数分析数据成功")

    @report_check
    def add_pca(self, pca_file, correlation_id=None):
        data_list = []
        correlation_id = ObjectId(correlation_id)
        with open(pca_file, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "correlation_id": correlation_id,
                    line[0]: line[1]
                }
                data_list.append(data)
        try:
            collection = self.db["sg_denovo_correlation_pca"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入sg_denovo_correlation_pca数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入sg_denovo_correlation_pca数据成功")

    @report_check
    def add_pca_rotation(self, input_file, db_name, correlation_id=None):
        data_list = []
        correlation_id = ObjectId(correlation_id)
        if db_name == "sg_denovo_correlation_pca_rotation":
            col_name = "gene_id"
        else:
            col_name = "species_name"
        with open(input_file, "r") as f:
            pcas = f.readline().strip().split("\t")[1:]
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "correlation_id": correlation_id,
                    col_name: line[0]
                }
                for n, p in enumerate(pcas):
                    data[p] = line[n + 1]
                data_list.append(data)
        try:
            collection = self.db[db_name]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s数据出错:%s" % db_name, e)
        else:
            self.bind_object.logger.info("导入%s数据成功" % db_name)
