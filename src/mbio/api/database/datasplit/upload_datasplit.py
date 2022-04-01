# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20190318

import os
import json
import types
import datetime
from bson.objectid import ObjectId
from biocluster.config import Config
from biocluster.api.database.base import Base, report_check

class UploadDatasplit(Base):
    """
    数据上传导表、信息检查
    """
    def __init__(self, bind_object):
        super(UploadDatasplit, self).__init__(bind_object)
        self._project_type = "datasplit"

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, types.StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, path):
        """
        检查文件是否存在,返回正确的路径
        """
        if os.path.exists(path):
            return path
        # elif os.path.exists(path.replace("ilustre", "clustre")):
        #     return path.replace("ilustre", "clustre")
        path = self.new_split_path(path)
        if os.path.exists(path):
            return path
        else:
            raise Exception('{}所指定的路径不存在，请检查'.format(path))

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path

    def check_lib_info(self, seq_board, lib_info):
        """
        检查上传的文库信息是否和系统里的文库信息一致，若是不一致，则不能进行导入
        """
        lib_info = self.check_exists(lib_info)
        results = self.db["sg_board"].find({"board_number": seq_board})
        if results != 1:
            return "在sg_board表里找到板:%s的记录不为一条,请检查" % seq_board
        board_id = results[0]["_id"]
        with open(lib_info, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lane_name = item[2]
                lib_name = item[0]
                lane_result = self.db["sg_board_lane"].find_one({"board_id": board_id, "lane_name": lane_name})
                if not lane_result:
                    return "没有在sg_board_lane表里找到板:%s对应的lane:%s,请检查" % seq_board, lane_name
                lane_id = lane_result["_id"]
                lib_results = self.db["sg_board_library"].find({"board_id": board_id, "lane_id": lane_id, "library_number": lib_name})
                if lib_results.count() != 1:
                    return "在sg_board_library表里找到板:%s对应的文库记录不为一条,请检查" % seq_board, lib_name
                if item[7] != lib_results[0]["library_type"]:
                    return "文库:%s的文库类型和系统里的信息不匹配，请检查" % lib_name
                if item[12] != lib_results[0]["index"]:
                    return "文库:%s的index和系统里的信息不匹配，请检查" % lib_name
                if item[13] != lib_results[0]["index_seq"]:
                    return "文库:%s的index序列和系统里的信息不匹配，请检查" % lib_name
                for f in item[16].split(";"):
                    self.check_exists(f)

    def check_sample_info(self, seq_board, sample_info):
        """
        检查上传的样本信息是否和系统里的信息一致
        """
        sample_info = self.check_exists(sample_info)
        results = self.db["sg_board"].find({"board_number": seq_board})
        if results != 1:
            return "在sg_board表里找到板:%s的记录不为一条,请检查" % seq_board
        board_id = results[0]["_id"]
        with open(sample_info, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lane_name = item[2]
                lib_name = item[0]
                lane_result = self.db["sg_board_lane"].find_one({"board_id": board_id, "lane_name": lane_name})
                if not lane_result:
                    return "没有在sg_board_lane表里找到板:%s对应的lane:%s,请检查" % seq_board, lane_name
                lane_id = lane_result["_id"]
                lib_results = self.db["sg_board_library"].find({"board_id": board_id, "lane_id": lane_id, "library_number": lib_name})
                if lib_results.count() != 1:
                    return "在sg_board_library表里找到板:%s对应的文库记录不为一条,请检查" % seq_board, lib_name
                if item[7] != lib_results[0]["library_type"]:
                    return "文库:%s的文库类型和系统里的信息不匹配，请检查" % lib_name
                for f in item[16].split(";"):
                    self.check_exists(f)
