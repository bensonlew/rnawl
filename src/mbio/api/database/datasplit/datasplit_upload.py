# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20210621
"""
蛋白代谢、三代数据上传导表
"""

import os
import re
import json
import types
import datetime
from bson.objectid import ObjectId
from collections import defaultdict
from biocluster.config import Config
from mbio.files.medical.html import HtmlFile
from biocluster.api.database.base import Base, report_check

class DatasplitUpload(Base):
    def __init__(self, bind_object):
        super(DatasplitUpload, self).__init__(bind_object)
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
        检查文件是否存在
        """
        if not os.path.exists(path):
            raise Exception("{}所指定的路径不存在，请检查".format(path))
        else:
            return True

    def update_coll_status_pacbio(self, status, lane_library_ids, desc=""):
        """
        更新sg_specimen_s3_pacbio表的status
        """
        lane_library_ids = lane_library_ids.split(";")
        for lane_library_id in lane_library_ids:
            query_dict = {"lane_library_id": lane_library_id}
            update_dict = {
                "status": status,
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "desc": desc
            }
            self.db["sg_specimen_s3_pacbio"].update(query_dict, {"$set": update_dict})
        self.bind_object.logger.info("更新sg_specimen_s3_pacbio表的status为{}".format(status))

    def update_coll_path_pacbio(self, project_table,md5sum=None):
        """
        更新sg_specimen_s3_pacbio表的work_raw_path、work_clean_path、raw_path、clean_path
        """
        self.check_exists(project_table)
        with open(project_table, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                query_dict = {"lane_library_id": item[1]}
                raw_bytes, clean_bytes = "", ""
                raw_md5sum,clean_md5sum = "",""
                if os.path.exists(item[-4]):
                    raw_bytes = str(os.path.getsize(item[-4]))
                    raw_md5sum = md5sum[item[-2]][os.path.basename(item[-2])]
                if os.path.exists(item[-3]):
                    clean_bytes = str(os.path.getsize(item[-3]))
                    clean_md5sum = md5sum[item[-1]][os.path.basename(item[-1])]
                update_dict = {
                    "work_raw_path": item[-4],
                    "work_clean_path": item[-3],
                    "raw_path": item[-2],
                    "clean_path": item[-1],
                    "raw_bytes": raw_bytes,
                    "clean_bytes": clean_bytes,
                    "clean_md5sum":clean_md5sum,
                    "raw_md5sum":raw_md5sum,
                    "status": "end",
                    "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "desc": "样本上传成功"
                }
                self.db["sg_specimen_s3_pacbio"].update(query_dict, {"$set": update_dict})
        self.bind_object.logger.info("更新sg_specimen_s3_pacbio表的path和status")

    def update_coll_status_protein(self, status, project_sn, task_sn, samples, desc=""):
        """
        更新sg_specimen_s3_protein表的status
        """
        samples_list = samples.split(";")
        results = self.db["sg_specimen_s3_protein"].find({"contract_sn": project_sn, "task_sn": task_sn})
        for result in results:
            if result["sample_name"] not in samples_list:
                continue
            query_dict = {"_id": result["_id"]}
            update_dict = {
                "status": status,
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "desc": desc
            }
            self.db["sg_specimen_s3_protein"].update(query_dict, {"$set": update_dict})
        self.bind_object.logger.info("更新sg_specimen_s3_protein表的status为{}".format(status))

    def update_coll_path_protein(self, project_sn, task_sn, sample_path, work_dir, s3_dir):
        """
        更新sg_specimen_s3_protein表的sample_ts_path、sample_s3_path
        """
        self.sample_info = {}
        with open(sample_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                path = work_dir
                s3_path = s3_dir if s3_dir.endswith("/") else s3_dir + "/"
                search_lib_dir, s3_search_lib_dir = "", ""
                if len(item) > 1 and item[1] != "":
                    path = item[1]
                if len(item) > 2 and item[2] != "":
                    s3_path = item[2]
                if len(item) > 3 and item[3] != "":
                    search_lib_dir = item[3]
                if len(item) > 4 and item[4] != "":
                    s3_search_lib_dir = item[4] if item[4].endswith("/") else item[4] + "/"
                self.sample_info[item[0]] = {
                    "path": path,
                    "s3_path": s3_path,
                    "search_lib_dir": search_lib_dir,
                    "s3_search_lib_dir": s3_search_lib_dir
                }
        results = self.db["sg_specimen_s3_protein"].find({"contract_sn": project_sn, "task_sn": task_sn})
        for result in results:
            query_dict = {"_id": result["_id"]}
            if result["sample_name"] in self.sample_info.keys():
                update_dict = {
                    "sample_ts_path": self.sample_info[result["sample_name"]]["path"],
                    "sample_s3_path": self.sample_info[result["sample_name"]]["s3_path"],
                    "search_lib_path": self.sample_info[result["sample_name"]]["s3_search_lib_dir"],
                    "search_lib_ts_path": self.sample_info[result["sample_name"]]["search_lib_dir"],
                    "status": "end",
                    "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "desc": "样本上传成功"
                }
            elif result["status"] == "start":
                    update_dict = {
                        "status": "end",
                        "sample_ts_path": work_dir,
                        "sample_s3_path": s3_dir if s3_dir.endswith("/") else s3_dir + "/",
                        "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "desc": "样本"+result["sample_name"]+"没有找到对应的s3路径"
                    }
                    self.bind_object.logger.info("样本{}没有找到对应的s3路径".format(result["sample_name"]))
            else:
                # update_dict = {"status": "end"}
                continue
            self.db["sg_specimen_s3_protein"].update(query_dict, {"$set": update_dict})
        self.bind_object.logger.info("更新sg_specimen_s3_protein表的path成功")

    def update_coll_qc_path_protein(self, project_sn, task_sn, qc_table):
        """
        更新sg_specimen_s3_protein_qc表的qc_ts_path、qc_s3_path
        """
        data_list = []
        with open(qc_table, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if len(item) > 1:
                    insert_data = {
                        "contract_sn": project_sn,
                        "task_sn": task_sn,
                        "sample_name": item[2],
                        "qc_ts_path": item[3],
                        "qc_s3_path": item[4],
                        "status": "end",
                        "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "desc": "qc文件上传成功"
                    }
                    data_list.append(insert_data)
        if len(data_list) > 0:
            self.db["sg_specimen_s3_protein_qc"].insert_many(data_list)
            self.bind_object.logger.info("更新ssg_specimen_s3_protein_qc表的path成功")

if __name__ == "__main__":
    a = DatasplitUpload(None)
    a.update_coll_path_pacbio(project_table="/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210826/Single_upload_pacbio_302900_20210826_111522/UploadPacbio/output/project_table.path.xls")
