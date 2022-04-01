# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.25

from biocluster.api.database.base import Base, report_check
from api_base import ApiBase
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import os
import re
import json
import math


class Qtl(ApiBase):
    def __init__(self, bind_object):
        """
        QTL定位分析导表
        """
        super(Qtl, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def add_sg_qtl(self, project_sn, task_id, params=None, name=None):
        """
        sg_qtl
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": "origin_qtl_analysis",
            "param": params if params else "null",
            "desc": "qtl定位主表"
        }
        main_id = self.db["sg_qtl"].insert_one(insert_data).inserted_id
        self.db["sg_qtl"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def update_sg_qtl(self, qtl_id, parent_source):
        """
        更新sg_qtl的oarent_source
        """
        qtl_id = self.check_objectid(qtl_id)
        self.db["sg_qtl"].update_one({"main_id": qtl_id}, {"$set": {"parent_source": parent_source}})
        self.bind_object.logger.info("更新sg_qtl的parent_source成功")

    def add_sg_qtl_result(self, qtl_id, qtl_path, data_type):
        """
        sg_qtl_result
        qtl_path: qtl_result.xls
        """
        qtl_id = self.check_objectid(qtl_id)
        self.check_exists(qtl_path)
        data_list = []
        marker_ids = {}
        trait_list = []
        with open(qtl_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "qtl_id": qtl_id,
                    "data_type": data_type,
                    "trait": item[0],
                    "chr": item[1],
                    "positon": round(float(item[2]), 4),
                    "lod": round(float(item[3]), 4),
                    "r2": round(float(item[4]), 4),
                    "start": round(float(item[5]), 4),
                    "end": round(float(item[6]), 4),
                    "marker_id": item[7],
                    "marker_start": item[8],
                    "marker_end": item[9]
                }
                data_list.append(insert_data)
                trait_list.append(item[0])
                if item[7] not in marker_ids.keys():
                    marker_ids[item[7]] = []
                marker_ids[item[7]].append(item[0])
        if data_list:
            self.db["sg_qtl_result"].insert_many(data_list)
        data_list1 = []
        for marker_id in marker_ids.keys():
            insert_data = {
                "qtl_id": qtl_id,
                "data_type": data_type,
                "marker_id": marker_id,
                "trait_number": len(marker_ids[marker_id]),
                "trait": ",".join(marker_ids[marker_id])
            }
            data_list1.append(insert_data)
        if data_list1:
            self.db["sg_qtl_pleiotropy"].insert_many(data_list1)
        trait_list = list(set(trait_list))
        self.db["sg_qtl"].update_one({"main_id": qtl_id}, {"$set": {"trait_list": trait_list}})
        return trait_list

    def add_sg_qtl_region(self, qtl_id, qtl_dir, data_type):
        """
        sg_qtl_region
        """
        qtl_id = self.check_objectid(qtl_id)
        self.check_exists(qtl_dir)
        for f in os.listdir(qtl_dir):
            if f.endswith(".assocation.xls"):
                trit_asso = os.path.join(qtl_dir, f)
                data_list = []
                with open(trit_asso, "r") as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        insert_data = {
                            "qtl_id": qtl_id,
                            "data_type": data_type,
                            "trait": item[0],
                            "chr": item[1],
                            "pos1": item[2],
                            "ref": item[3],
                            "alt": item[4],
                            "type": item[5],
                            "f1_genotype": item[6],
                            "f1_depth": item[7],
                            "m1_genotype": item[8],
                            "m1_depth": item[9]
                        }
                        data_list.append(insert_data)
                if data_list:
                    self.col_insert_data("sg_qtl_region", data_list)
                else:
                    self.bind_object.logger.info("定位到的区域在pop.final.vcf里没有结果,请检查")

    def add_sg_box_region(self, qtl_id, qtl_dir, data_type):
        """
        关联区域详情箱线图
        """
        qtl_id = self.check_objectid(qtl_id)
        self.check_exists(qtl_dir)
        result = self.db["sg_qtl"].find_one({"_id": qtl_id})
        if not result:
            raise Exception("没找到sg_qtl记录，请检查")
        task_id = result["task_id"]
        mark_list, mark_id = [], ""
        for f in os.listdir(qtl_dir):
            if f.endswith(".box.xls"):
                data_list = []
                box_path = os.path.join(qtl_dir, f)
                trit = f.split(".box.xls")[0]
                with open(box_path, "r") as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        if item[0] not in mark_list:
                            mark_list.append(item[0])
                            try:
                                result = self.db["sg_qtl_region"].find_one({"qtl_id": qtl_id, "trait": trit, "chr": item[0].split("_")[0], "pos1": item[0].split("_")[1]})
                                origin_id = result["_id"]
                            except:
                                continue
                            if mark_id != "":
                                self.db["sg_boxplot"].update_one({"_id": boxplot_id}, {"$set": {"categories": categories}})
                            categories = []
                            mark_id = item[0]
                            boxplot_id = self.add_sg_box(origin_id, task_id, mark_id, categories, data_type)
                        box_data = []
                        # box_data.append(float(item[6]) if item[6] == "nan" else None)
                        # box_data.append(float(item[5]) if item[5] == "nan" else None)
                        # box_data.append(float(item[4]) if item[4] == "nan" else None)
                        # box_data.append(float(item[3]) if item[3] == "nan" else None)
                        # box_data.append(float(item[2]) if item[2] == "nan" else None)
                        box_data.append(None if math.isnan(float(item[6])) else float(item[6]))
                        box_data.append(None if math.isnan(float(item[5])) else float(item[5]))
                        box_data.append(None if math.isnan(float(item[4])) else float(item[4]))
                        box_data.append(None if math.isnan(float(item[3])) else float(item[3]))
                        box_data.append(None if math.isnan(float(item[2])) else float(item[2]))
                        categories.append(item[1])
                        try:
                            scatter_data_ = json.loads(item[7])
                            scatter_data = {}
                            for key in scatter_data_:
                                if math.isnan(scatter_data_[key]):
                                    scatter_data[key] = None
                                else:
                                    scatter_data[key] = scatter_data_[key]
                        except:
                            scatter_data = {}
                        insert_data = {
                            "boxplot_id": boxplot_id,
                            "name": item[1],
                            "value": box_data,
                            "scatter_data": scatter_data
                        }
                        data_list.append(insert_data)
                if not data_list:
                    self.bind_object.logger.info("箱线图数据为空")
                else:
                    self.col_insert_data("sg_boxplot_detail", data_list)

    def add_sg_box(self, origin_id, task_id, mark_id, categories, data_type):
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "origin_id": origin_id,
            "task_id": task_id,
            "data_type": data_type,
            "name": mark_id,
            "categories": categories,
            "type": 2,
            "location": "qtl_box",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        boxplot_id = self.db["sg_boxplot"].insert_one(insert_data).inserted_id
        return boxplot_id

    def add_sg_manhattan(self, qtl_id, trait_list, qtl_dir, data_type):
        """
        sg_manhattan, sg_manhattan_detail
        """
        qtl_id = self.check_objectid(qtl_id)
        location = "qtl_result"
        chr_data, trit_lod = {}, {}
        data_list = []
        for trait in trait_list:
            qtl_csv = os.path.join(qtl_dir, trait + ".qtl.result")
            scan_csv = os.path.join(qtl_dir, trait + ".scan.csv")
            self.check_exists(qtl_csv)
            self.check_exists(scan_csv)
            with open(qtl_csv, "r") as f:
                lines = f.readlines()
                item = lines[1].strip().split("\t")
                pm1 = float(item[5])
                trit_lod[trait] = float(item[5])
            chr_list, loca_list = [], []
            trit_value = []
            chr_data[trait] = {}
            trait_data = {}
            chr_pos = {}
            with open(scan_csv, "r") as f:
                lines = f.readlines()
                for i in range(1, len(lines)):
                    item = lines[i].strip().split("\t")
                    chr = item[0].split("_")[0].split('"')[1]
                    if chr not in chr_list:
                        if i != 1:
                            item1 = lines[i-1].strip().split("\t")
                            start_pos = float(item1[2])
                            chr_ = item1[0].split("_")[0].split('"')[1]
                            chr_pos[chr_] = start_pos
                        chr_list.append(chr)
                        chr_data[trait][chr] = []
                        trait_data[chr] = []
                    chr_pos[chr] = float(item[2])
                    chr_data[trait][chr].append([float(item[2]), float(item[3])])
                    trait_data[chr].append([float(item[2]), float(item[3])])
            chr_list = self.sort_chr(chr_list)
            insert_data = {
                "origin_id": qtl_id,
                "data_type": data_type,
                "name": trait,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "location": "qtl_result",
                "type": 2,
                "chr_list": chr_list,
                "threshold_name": "LOD",
                "threshold_value": pm1
            }
            main_id = self.db["sg_manhattan"].insert_one(insert_data).inserted_id
            for chr in chr_list:
                insert_data = {
                    "manhattan_id":  main_id,
                    "pos": chr_pos[chr],
                    "name": chr,
                    "value": trait_data[chr]
                }
                self.db["sg_manhattan_detail"].insert_one(insert_data)

    def add_sg_distribution(self, qtl_id, total_map, qtl_dir, type=None, data_type="total"):
        """
        一因多效图
        """
        qtl_id = self.check_objectid(qtl_id)
        self.check_exists(total_map)
        self.check_exists(qtl_dir)
        result = self.db["sg_qtl"].find_one({"_id": qtl_id})
        if not result:
            raise Exception("没找到sg_qtl记录，请检查")
        task_id = result["task_id"]
        # chr_list = []
        # data_dict, snp_dict = {}, {}
        # if type == "map":
        #     with open(total_map, "r") as f:
        #         for line in f:
        #             if line.startswith("group"):
        #                 continue
        #             item = line.strip().split("\t")
        #             mark_id = item[0]
        #             chr = item[0].split("_")[0]
        #             if chr not in chr_list:
        #                 chr_list.append(chr)
        #                 data_dict[chr] = []
        #                 snp_dict[chr] = []
        #             data_dict[chr].append(round(float(item[2]), 4))
        # else:
        #     with open(total_map, "r") as f:
        #         head = f.readline()
        #         for line in f:
        #             item = line.strip().split(",")
        #             mark_id = item[0]
        #             chr = item[0].split("_")[0]
        #             if chr not in chr_list:
        #                 chr_list.append(chr)
        #                 data_dict[chr] = []
        #                 snp_dict[chr] = []
        #             data_dict[chr].append(round(float(item[2]), 4))
        # for f in os.listdir(qtl_dir):
        #     if f.endswith(".qtl.result"):
        #         qtl_result = os.path.join(qtl_dir, f)
        #         with open(qtl_result, "r") as f:
        #             trit = os.path.basename(qtl_result).split(".qtl.result")[0]
        #             lines = f.readlines()
        #             for line in lines[1:]:
        #                 item = line.strip().split("\t")
        #                 chr = item[0].split("_")[0]
        #                 pos = round(float(item[2]), 4)
        #                 if chr not in snp_dict.keys():
        #                     continue
        #                 snp_dict[chr].append([pos, trit])
        # start, end = 0, 0
        # for chr in data_dict:
        #     if data_dict[chr][-1] > end:
        #         end = data_dict[chr][-1]
        # end = math.ceil(end)
        # name = ""
        # location = "qtl_distribution"
        # chr_list = self.sort_chr(chr_list)
        # main_id = self.sg_distribution(task_id, qtl_id, name, "2", start, end, location, chr_list, data_type)
        # for chr in chr_list:
        #     snp_value = snp_dict[chr]
        #     self.sg_distribution_detail(main_id, chr, data_dict[chr], snp_value)
        # self.db["sg_qtl"].update_one({"_id": qtl_id}, {"$set": {"chr_list": chr_list}})
        lds, chrs_, scas_, all_value = [], [], [], []
        dis_min, dis_max = 10000000000000000, 0
        with open(total_map, 'r') as r:
            for line in r:
                temp = line.strip().split('\t')
                if temp[0] == 'group':
                    lds.append(int(temp[1]))
                    all_value.append({temp[0]: int(temp[1])})
                else:
                    temp_ = temp[0].split("_")[0]
                    if re.match(r"^chr.*", temp_):
                        if temp_ not in chrs_:
                            chrs_.append(int(temp_[3:]))
                    elif re.match(r'^sca.*', temp_):
                        if temp_ not in scas_:
                            scas_.append(int(temp_[3:]))
                    if float("%.4f" % float(temp[1])) > dis_max:
                        dis_max = float("%.4f" % float(temp[1]))
                    if float("%.4f" % float(temp[1])) < dis_min:
                        dis_min = float("%.4f" % float(temp[1]))
                    all_value.append({temp[0]: temp[1]})
        cat_lds = list(set(lds))
        cat_lds.sort()
        chrs_.sort()
        scas_.sort()
        chrs = list(set(chrs_))
        chr_sca = self.set_chr_sca(chrs, [])
        name = ""
        location = "qtl_distribution"
        main_id = self.sg_distribution(task_id, qtl_id, name, "2", dis_min, dis_max, location, ["LG"+str(i) for i in cat_lds], data_type)
        value_dict, marker_dict, snp_dict = {}, {}, {}
        for m in range(0, len(lds)):
            value = []
            tem_value = []
            marker_value = []
            if m == len(lds) - 1:
                if {"group": lds[m]} in all_value:
                    tem_value = all_value[all_value.index({"group": lds[m]}) + 1:]
            else:
                if {"group": lds[m]} in all_value and {"group": lds[m + 1]} in all_value:
                    tem_value = all_value[all_value.index({"group": lds[m]}) + 1: all_value.index(
                        {"group": lds[m + 1]})]
            if tem_value:
                for n in tem_value:
                    value.append([n.keys()[0].split("_")[0], float("%.4f" % float(n.values()[0]))])
                    marker_value.append(n.keys()[0])
            value_dict[lds[m]] = value
            marker_dict[lds[m]] = marker_value
            snp_dict[lds[m]] = []
            m += 1
        for f in os.listdir(qtl_dir):
            if f.endswith(".qtl.result"):
                qtl_result = os.path.join(qtl_dir, f)
                with open(qtl_result, "r") as f:
                    trit = os.path.basename(qtl_result).split(".qtl.result")[0]
                    lines = f.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        pos = round(float(item[2]), 4)
                        for lg in cat_lds:
                            if item[0] in marker_dict[lg]:
                                snp_dict[lg].append([pos, trit])
        for n in cat_lds:
            if n in value_dict.keys():
                self.sg_distribution_detail(main_id, 'LG{}'.format(n), value_dict[n], snp_dict[n])
        self.db["sg_qtl"].update_one({"_id": qtl_id}, {"$set": {"chr_list": chr_sca, "lgs": ["LG"+str(i) for i in cat_lds]}})

    def sort_chr(self, chr_lists):
        """
        对chr/sca进行排序
        """
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:   # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        return final_list

    def set_chr_sca(self, chrs_, scas_):
        """
        合并chr与sca，完成排序
        add by hongdong@20180705
        :return:
        """
        chr_sca = []
        for m in chrs_:
            chr_sca.append("chr{}".format(m))
        if scas_:
            for n in scas_:
                chr_sca.append("sca{}".format(n))
        return chr_sca


if __name__ == "__main__":
    a = Qtl(None)
    project_sn = "gmap_test"
    task_id = "tsanger_30729"
    # qtl_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/qtl/3-15.xls"
    # qtl_id = a.add_sg_qtl(project_sn, task_id)
    # a.add_sg_qtl_result(qtl_id, qtl_path)
    # qtl_id = "5b30ad7ea4e1af06ffe7b1da"
    # qtl_id = "5b30ad7ea4e1af06ffe7b1d1"
    # a.add_sg_qtl_region(qtl_id, qtl_path)
    # total_map = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/qtl/total.map"
    # qtl_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180709/Single_qtl_analysis6/QtlAnalysis/output"
    # qtl_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/qtl"
    # a.add_sg_distribution(qtl_id, task_id, total_map, qtl_dir)
    # trait_list = ["GL", "GW"]
    # qtl_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180713/Single_qtl_analysis1/QtlAnalysis/output"
    # qtl_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180718/Single_qtl_analysis2/QtlAnalysis/output"
    # a.add_sg_manhattan(qtl_id, trait_list, qtl_dir)
    # qtl_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180713/Single_qtl_analysis1/QtlAnalysis/output/test"
    # a.add_sg_qtl_region(qtl_id, qtl_dir)
    # a.add_sg_box_region(qtl_id, task_id, qtl_dir)
    # qtl_id = "5b5fed97a4e1af643ac083c6"
    # total_map = "/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/LinkageGrouping/LinkageGrouping_by_hand_20180730_180925508156/evalutaion/total.csv"
    qtl_dir = "/mnt/lustre/users/sanger/workspace/20181227/Single_i-sanger_145158_QtlAnalysis_1227165611637150/QtlAnalysis/output"
    type = "csv"
    data_type = "total"
    # a.add_sg_distribution(qtl_id=qtl_id, total_map=total_map, qtl_dir=qtl_dir, type=type, data_type=data_type)
    qtl_id = "5c2493abec02cc67917a59a7"
    a.add_sg_box_region(qtl_id, qtl_dir, data_type)
