# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20180624

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

class DatasplitMerge(Base):
    def __init__(self, bind_object):
        super(DatasplitMerge, self).__init__(bind_object)
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

    def update_sg_split_specimen_merge(self, merge_id, output_dir, s3_upload_dir, operation_type):
        """
        更新合并样本的信息
        """
        self.check_exists(output_dir)
        merge_id = self.check_objectid(merge_id)
        md5sum_file = os.path.join(output_dir,"md5sum.txt")
        md5sum_dict = {}
        if os.path.exists(md5sum_file):
            with open(md5sum_file,"r") as mf:
                while 1:
                    line = mf.readline()
                    if not line:
                        break
                    fd = line.rstrip().split("  ")
                    md5sum_dict[fd[1]] = fd[0]
        raw_path, raw75_path, clean_path, raw_bytes, raw75_bytes, clean_bytes, raw_md5sum, clean_md5sum, raw75_md5sum = [], [], [], [], [], [], [], [], []
        for f in os.listdir(output_dir):
            path = os.path.join(output_dir, f)
            s3_path = os.path.join(s3_upload_dir, f)
            if f.endswith("R1.raw.fastq.gz"):
                raw_md5sum.insert(0,md5sum_dict[f])
                raw_path.insert(0, s3_path)
                raw_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("R2.raw.fastq.gz"):
                raw_md5sum.append(md5sum_dict[f])
                raw_path.append(s3_path)
                raw_bytes.append(str(os.path.getsize(path)))
            elif f.endswith("clean.1.fastq.gz"):
                clean_md5sum.insert(0, md5sum_dict[f])
                clean_path.insert(0, s3_path)
                clean_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("clean.2.fastq.gz"):
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith(".clean.fastq.gz"):  # 多样性
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith("clean.1.fastq.gz"):  # 多样性
                clean_md5sum.insert(0, md5sum_dict[f])
                clean_path.insert(0, s3_path)
                clean_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("clean.2.fastq.gz"):  # 多样性
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith(".fasta.gz"):
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith(".fq.gz"):
                raw75_md5sum.append(md5sum_dict[f])
                raw75_path.append(s3_path)
                raw75_bytes.append(str(os.path.getsize(path)))
        self.bind_object.logger.info(";".join(clean_path))
        if operation_type == "merge":
            update_dict = {
                "raw75_md5sum":";".join(raw75_md5sum),
                "raw_md5sum":";".join(raw_md5sum),
                "clean_md5sum":";".join(clean_md5sum),
                "raw_path": ";".join(raw_path),
                "raw75_path": ";".join(raw75_path),
                "clean_path": ";".join(clean_path),
                "raw_bytes": ";".join(raw_bytes),
                "raw75_bytes": ";".join(raw75_bytes),
                "clean_bytes": ";".join(clean_bytes),
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            self.db["sg_split_specimen_merge"].update({"_id": merge_id}, {"$set": update_dict})
        elif operation_type == "rename":
            result = self.db["sg_split_specimen_rename"].find_one({"_id": merge_id})
            update_dict = {
                "re_raw_path": result["raw_path"],
                "re_clean_path": ";".join(clean_path),
                "re_raw_bytes": result["raw_bytes"],
                "re_clean_bytes": ";".join(clean_bytes),
                "re_raw_md5sum":result["raw_md5sum"] if result.has_key("raw_md5sum") else "",
                "re_clean_md5sum":";".join(clean_md5sum),
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            self.db["sg_split_specimen_rename"].update({"_id": merge_id}, {"$set": update_dict})
        self.bind_object.logger.info("导表结束")

if __name__ == "__main__":
    a = DatasplitMerge(None)
    # merge_id =
    # a.update_sg_split_specimen_merge(merge_id, self.output_dir, s3_output_dir, self.option("operation_type"))
