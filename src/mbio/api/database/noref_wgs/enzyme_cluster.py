# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20181220

import os
import json
import types
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check


class EnzymeCluster(Base):
    """
    无参WGS导表：酶切片段聚类分析
    """
    def __init__(self, bind_object):
        super(EnzymeCluster, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"

    def check_exists(self, file):
        """
        检查file文件是否存在
        """
        if not os.path.exists(file):
            raise Exception("文件：%s不存在，请检查" % file)
        return True

    def check_objectid(self, id):
        """
        检查id是否是mongo的_id
        """
        if not isinstance(id, ObjectId):
            if isinstance(id, types.StringTypes):
                id = ObjectId(id)
            else:
                raise Exception("id:%s必须为ObjectID或其对应的字符串" % id)
        return id

    def add_sg_cluster_tag(self, task_id, params=None):
        """
        sg_cluster_tag
        """
        if params:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            params = {}
        insert_data = {
            "task_id": task_id,
            "name": "TagStat",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "desc": "tag信息统计主表"
        }
        main_id = self.db["sg_cluster_tag"].insert_one(insert_data).inserted_id
        self.db["sg_cluster_tag"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        print "sg_cluster_tag success!"
        return main_id

    def add_sg_cluster_tag_stat(self, cluster_tag_id, tag_path, type="stacks"):
        """
        sg_cluster_tag_stat
        """
        self.check_exists(tag_path)
        cluster_tag_id = self.check_objectid(cluster_tag_id)
        data_list = []
        with open(tag_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "cluster_tag_id": cluster_tag_id,
                    "specimen_id": item[0]
                }
                if type == 'stacks':
                    insert_data.update({
                        "tag_num": int(item[1]),
                        "average_length": float(item[4]),
                        "average_depth": float(item[3]),
                        "tag_cover_5": round(float(item[-2]) / float(item[1]), 4) * 100,
                        "tag_cover_10": round(float(item[-1]) / float(item[1]), 4) * 100
                    })
                else:
                    insert_data.update({
                        "tag_num": int(item[1]),
                        "average_length": float(item[2]),
                        "average_depth": float(item[3]),
                        "tag_cover_5": float(item[-2]),
                        "tag_cover_10": float(item[-1])
                    })
                data_list.append(insert_data)
        self.db["sg_cluster_tag_stat"].insert_many(data_list)
        print "sg_cluster_tag_stat success!"

    def add_sg_tag_curve(self, task_id, cluster_tag_id, specimen_id, depth_path):
        """
        tag深度分布导表
        """
        self.check_exists(depth_path)
        cluster_tag_id = self.check_objectid(cluster_tag_id)
        categories, value_list = [], []
        with open(depth_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[0] == specimen_id:
                    categories.append(item[1])
                    value_list.append(int(item[2]))
        curve_id = self.add_sg_curve(task_id, cluster_tag_id, specimen_id, categories, "tag_curve", specimen_id, types=1)
        self.add_sg_curve_detail(curve_id, specimen_id, value_list)

    def add_sg_cluster_consensus(self, task_id, params=None):
        """
        sg_cluster_consensus
        """
        insert_data = {
            "task_id": task_id,
            "name": "ClusterAnalysis",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "desc": "聚类分析主表"
        }
        main_id = self.db["sg_cluster_consensus"].insert_one(insert_data).inserted_id
        self.db["sg_cluster_consensus"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        print "sg_cluster_consensus success!"
        return main_id

    def add_sg_cluster_consensus_stat(self, cluster_consensus_id, consensus_path):
        """
        sg_cluster_consensus_stat
        """
        self.check_exists(consensus_path)
        cluster_consensus_id = self.check_objectid(cluster_consensus_id)
        data_list = []
        with open(consensus_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "cluster_consensus_id": cluster_consensus_id,
                    "name": item[0],
                    "value": item[1]
                }
                data_list.append(insert_data)
        self.db["sg_cluster_consensus_stat"].insert_many(data_list)
        print "sg_cluster_consensus_stat success!"

    def add_sg_cluster_bar(self, cluster_consensus_id, bar_path):
        """
        Consensus聚类覆盖度分布图
        """
        self.check_exists(bar_path)
        cluster_consensus_id = self.check_objectid(cluster_consensus_id)
        categories, value = [], []
        with open(bar_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                categories.append(item[0])
                value.append(float(item[1]))
        bar_id = self.sg_bar(task_id, cluster_consensus_id, "consensus_bar", categories, "1", "consensus_bar")
        self.sg_bar_detail(bar_id, "consensus_bar", value)

    def add_sg_curve(self, task_id, origin_id, name, categories, location, sample, types=1):
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "other_attr": "",
            "ext_name": {"type": "sample", "title": sample},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        coll = self.db["sg_curve"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_curve success!"
        return main_id

    def add_sg_curve_detail(self, curve_id, name, value):
        curve_id = self.check_objectid(curve_id)
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": value
        }
        coll = self.db["sg_curve_detail"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_curve_detail success!"

    def sg_bar(self, task_id, origin_id, name, categories, types, location=None, other_attr=None, ext_name=None):
        origin_id = self.check_objectid(origin_id)
        insert_data = {
            "task_id": task_id,
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location if location else "",
            "other_attr": other_attr if other_attr else "",
            "ext_name": ext_name if ext_name else "",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        main_id = self.db["sg_bar"].insert_one(insert_data).inserted_id
        return main_id

    def sg_bar_detail(self, bar_id, name, value):
        bar_id = self.check_objectid(bar_id)
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": value
        }
        self.db["sg_bar_detail"].insert_one(insert_data)

if __name__ == "__main__":
    a = EnzymeCluster(None)
    task_id = "noref_test1"
    cluster_tag_id = a.add_sg_cluster_tag(task_id)
    tag_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/tags_stat.1.txt"
    a.add_sg_cluster_tag_stat(cluster_tag_id, tag_path)
    depth_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/DEP.txt"
    specimen_id = "1080"
    a.add_sg_tag_curve(task_id, cluster_tag_id, specimen_id, depth_path)
    specimen_id = "1145"
    a.add_sg_tag_curve(task_id, cluster_tag_id, specimen_id, depth_path)
    specimen_id = "c29875"
    a.add_sg_tag_curve(task_id, cluster_tag_id, specimen_id, depth_path)
    specimen_id = "1828"
    a.add_sg_tag_curve(task_id, cluster_tag_id, specimen_id, depth_path)
    cluster_consensus_id = a.add_sg_cluster_consensus(task_id)
    consensus_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/consensus_stat.txt"
    a.add_sg_cluster_consensus_stat(cluster_consensus_id, consensus_path)
    bar_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/bar.txt"
    a.add_sg_cluster_bar(cluster_consensus_id, bar_path)
