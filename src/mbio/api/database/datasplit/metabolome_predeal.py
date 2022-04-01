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
from biocluster.api.database.base import Base, report_check

class MetabolomePredeal(Base):
    def __init__(self, bind_object):
        super(MetabolomePredeal, self).__init__(bind_object)
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

    def update_sg_predeal(self, predeal_id, output_dir, s3_output_dir, status, desc, log):
        """
        更新sg_predeal表
        """
        md5sum_file = os.path.join(output_dir,"md5sum.txt")
        md5sum_dict = {}
        with open(md5sum_file,"r") as m5:
            while 1:
                line = m5.readline()
                if not line:
                    break 
                fd = line.rstrip().split("  ")
                md5sum_dict[fd[1]] = fd[0]
        predeal_id = self.check_objectid(predeal_id)
        pos_path, neg_path, group_path, compare_path, deal_pos_path, deal_neg_path = [''] * 6
        pos_md5sum, neg_md5sum, group_md5sum, compare_md5sum, deal_pos_md5sum, deal_neg_md5sum = [''] * 6
        s3_pos_path, s3_neg_path, s3_group_path, s3_compare_path, s3_deal_pos_path, s3_deal_neg_path = [''] * 6
        for file in os.listdir(output_dir):
            if file in ['group.txt', 'group_auto.txt']:
                group_path = os.path.join(output_dir, file)
                s3_group_path = os.path.join(s3_output_dir, file)
                group_md5sum = md5sum_dict[file]
            if file in ['control.txt', 'control_auto.txt', 'compare.txt']:
                compare_path = os.path.join(output_dir, file)
                s3_compare_path = os.path.join(s3_output_dir, file)
                compare_md5sum = md5sum_dict[file]
            if "Dealed_pos_measurement.xls" in file:
                deal_pos_path = os.path.join(output_dir, file)
                s3_deal_pos_path = os.path.join(s3_output_dir, file)
                deal_pos_md5sum = md5sum_dict[file]
            elif "Dealed_neg_measurement.xls" in file:
                deal_neg_path = os.path.join(output_dir, file)
                s3_deal_neg_path = os.path.join(s3_output_dir, file)
                deal_neg_md5sum = md5sum_dict[file]
            elif "pos_measurement.xls" in file:
                pos_path = os.path.join(output_dir, file)
                s3_pos_path = os.path.join(s3_output_dir, file)
                pos_md5sum = md5sum_dict[file]
            elif "neg_measurement.xls" in file:
                neg_path = os.path.join(output_dir, file)
                s3_neg_path = os.path.join(s3_output_dir, file)
                neg_md5sum = md5sum_dict[file]
        query_dict = {"_id": predeal_id}
        update_dict = {
            "pos_path": s3_pos_path,
            "neg_path": s3_neg_path,
            "group_path": s3_group_path,
            "compare_path": s3_compare_path,
            "oreder_pos_path": s3_deal_pos_path,
            "oreder_neg_path": s3_deal_neg_path,
            "group_md5sum": group_md5sum,
            "compare_md5sum": compare_md5sum,
            "oreder_neg_md5sum": deal_neg_md5sum,
            "pos_md5sum": pos_md5sum,
            "neg_md5sum": neg_md5sum,
            "oreder_pos_md5sum":deal_pos_md5sum,
            "pos_size": os.path.getsize(pos_path) if os.path.exists(pos_path) else 0,
            "neg_size": os.path.getsize(neg_path) if os.path.exists(neg_path) else 0,
            "group_size": os.path.getsize(group_path) if os.path.exists(group_path) else 0,
            "compare_size": os.path.getsize(compare_path) if os.path.exists(compare_path) else 0,
            "oreder_pos_size": os.path.getsize(deal_pos_path) if os.path.exists(deal_pos_path) else 0,
            "oreder_neg_size": os.path.getsize(deal_neg_path) if os.path.exists(deal_neg_path) else 0,
            "status": status,
            "desc": desc,
            "log": log,
            "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        self.db["sg_predeal"].update(query_dict, {"$set": update_dict})
        self.bind_object.logger.info("更新sg_predeal表成功")


if __name__ == "__main__":
    a = MetabolomePredeal(None)
    a.update_sg_predeal(project_table="/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210826/Single_upload_pacbio_302900_20210826_111522/UploadPacbio/output/project_table.path.xls")
