# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.14

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import os
import json
import re


class BinnerCalculate(Base):
    def __init__(self, bind_object):
        """
        Bin Marker分析导表
        """
        super(BinnerCalculate, self).__init__(bind_object)
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

    def add_sg_binmarker(self, project_sn, task_id, params=None, name=None):
        """
        sg_binmarker
        """
        if params:
            new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            new_params = "null"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_bin_marker",
            "params": new_params,
            "desc": "BinMarker分析主表"
        }
        main_id = self.db["sg_binmarker"].insert_one(insert_data).inserted_id
        self.db["sg_binmarker"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def add_sg_binmarker_bin(self, binmarker_id, bin_stat):
        """
        sg_binmarker_bin
        bin_stat: bin_stat.xls
        """
        binmarker_id = self.check_objectid(binmarker_id)
        self.check_exists(bin_stat)
        data_list = []
        chr_list = []
        sca_list = []
        with open(bin_stat, "r") as f:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            if len(head) == 4:
                type = ["snp", "indel"]
            elif head[2].startswith("SNP"):
                type = ["snp"]
            else:
                type = ["indel"]
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "binmarker_id": binmarker_id,
                    "chr": item[0],
                    "bin_num": int(item[1])
                }
                if len(type) == 2:
                    insert_data["snp_num"] = float(item[2])
                    insert_data["indel_num"] = float(item[3])
                elif type[0] == "snp":
                    insert_data["snp_num"] = float(item[2])
                else:
                    insert_data["indel_num"] = float(item[2])
                if item[0].lower().startswith("chr"):
                    chr_list.append(insert_data)
                elif item[0].lower().startswith("sca"):
                    sca_list.append(insert_data)
            if len(chr_list) > 0:
                try:
<<<<<<< HEAD
                    chr_list.sort(key=lambda i: int(re.findall("\d+", i["chr"])[0]))
                except:
                    chr_list.sort()
                data_list.extend(chr_list)
            if len(sca_list) > 0:
                try:
                    sca_list.sort(key=lambda i: int(re.findall("\d+", i["chr"])[0]))
                except:
                    sca_list.sort()
=======
                    chr_list.sort(key=lambda i: int(re.findall("\d+", i)[0]))
                except:
                    pass
                data_list.extend(chr_list)
            if len(sca_list) > 0:
                try:
                    sca_list.sort(key=lambda i: int(re.findall("\d+", i)[0]))
                except:
                    pass
>>>>>>> hd_s3_1024
                data_list.extend(sca_list)
        self.db["sg_binmarker_bin"].insert_many(data_list)
        self.db["sg_binmarker"].update({"main_id": binmarker_id}, {"$set": {"type": type}})

    def add_sg_binmarker_var(self, binmarker_id, bin_info):
        """
        sg_binmarker_var
        bin_info: bin_info.xls
        """
        binmarker_id = self.check_objectid(binmarker_id)
        self.check_exists(bin_info)
        data_list = []
        with open(bin_info, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "binmarker_id": binmarker_id,
                    "marker_id": item[0],
                    "chr": item[1],
                    "pos": item[2],
                    "variant_num": int(item[3]),
                    "type": item[4],
                    "nind": int(item[5]),
                    "nmiss": int(item[6]),
                    "geno": item[7],
                    "sgeno": item[8],
                    "signif": self.change_type(item[9])
                }
                try:
                    insert_data["segretion"] = item[10]
                except:
                    insert_data["segretion"] = "--"
                data_list.append(insert_data)
        self.db["sg_binmarker_var"].insert_many(data_list)

    def change_type(self, data):
        try:
            temp = float(data)
        except:
            temp = data
        return temp

    def add_sg_binmarker_var_detail(self, binmarker_id, bin_pos):
        """
        sg_binmarker_var_detail
        """
        binmarker_id = self.check_objectid(binmarker_id)
        self.check_exists(bin_pos)
        data_list = []
        marker_id = ""
        with open(bin_pos, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[0] != marker_id:
                    marker_id = item[0]
                    result = self.db["sg_binmarker_var"].find_one({"binmarker_id": binmarker_id, "marker_id": item[0]})
                    if not result:
                        raise Exception("没找到marker_id为{}对应的标记信息".format(item[0]))
                var_id = result["_id"]
                insert_data = {
                    "var_id": var_id,
                    "marker_id": item[0],
                    "chr": item[1],
                    "pos": item[2],
                    "ref": item[3],
                    "alt": item[4]
                }
                data_list.append(insert_data)
        if data_list:
            self.db["sg_binmarker_var_detail"].insert_many(data_list)
        else:
            self.bind_object.logger.info("bin_pos.xls结果为空")

    def update_sg_binmarker_path(self, binmarker_id, filter_marker_path):
        """
        更新sg_binmarker表里的filter_marker_path、detail_info_path路径
        filter_marker_path: Total.bin.marker
        filter_marker_path: pop.filtered.detail.info
        """
        binmarker_id = self.check_objectid(binmarker_id)
        result = self.db["sg_binmarker"].find_one({"main_id": binmarker_id})
        if not result:
            raise Exception("没有在表sg_binmarker中找到main_id：{}的结果，请检查".format(binmarker_id))
        params = json.loads(result["params"])
        matrix_id = params["matrix_id"]
        matrix_id = self.check_objectid(matrix_id)
        marker_result = self.db["sg_marker"].find_one({"main_id": matrix_id})
        if not marker_result:
            raise Exception("没有在表sg_marker中找到main_id：{}的结果，请检查".format(matrix_id))
        detail_info_path = marker_result["detail_info_path"]
        self.db["sg_binmarker"].update_one({"main_id": binmarker_id}, {"$set": {"detail_info_path":\
                                            detail_info_path, "filtered_marker_path": filter_marker_path}})

if __name__ == "__main__":
    a = BinnerCalculate(None)
    project_sn = "gmap_test"
    task_id = "tsanger_30729"
    # bin_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180612/Single_binner_calculate2/BinnerCalculate/output/bin_stat.xls"
    # bin_info = "/mnt/ilustre/users/sanger-dev/workspace/20180612/Single_binner_calculate2/BinnerCalculate/output/bin_info.xls"
    # bin_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180621/Single_bin_marker2/BinnerCalculate/output/bin_stat.xls"
    # bin_info = "/mnt/ilustre/users/sanger-dev/workspace/20180621/Single_bin_marker2/BinnerCalculate/output/bin_info.xls"
    # bin_pos = "/mnt/ilustre/users/sanger-dev/workspace/20180621/Single_bin_marker2/BinnerCalculate/output/bin_pos.xls"
    # binmarker_id = a.add_sg_binmarker(project_sn, task_id)
    # a.add_sg_binmarker_bin(binmarker_id, bin_stat)
    # a.add_sg_binmarker_var(binmarker_id, bin_info)
    # a.add_sg_binmarker_var_detail(binmarker_id, bin_pos)
    binmarker_id = "5b4ed523a4e1af55d4b8d7aa"
    filter_marker_path = "rerewrweset/files/gmap_test/gmap_test/gmap_test_1/hongdong/0716_bin/03.binner/Total.bin.marker"
    a.update_sg_binmarker_path(binmarker_id, filter_marker_path)
