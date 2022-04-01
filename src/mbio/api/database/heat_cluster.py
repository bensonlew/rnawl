# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# last modified by guhaidong 20171109

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import re
# from biocluster.config import Config


class HeatCluster(Base):
    def __init__(self, bind_object):
        super(HeatCluster, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    def create_newick_table(self, params, linkage, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51003501")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录", code="51003502")
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "物种聚类树"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "table_type": "otu",
            "tree_type": "otu_cluster",
            "hcluster_method": linkage,
            "status": "end",
            "table_id": from_otu_table,
            "desc": desc,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "cluster_newick",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db['sg_newick_tree']
        inserted_id = collection.insert_one(insert_data).inserted_id
        #collection.update({"_id":inserted_id},{"$set": {"main_id" :inserted_id}})
        return inserted_id

    @report_check
    def update_newick(self, path, newick_id):
        if not isinstance(newick_id, ObjectId):
            if isinstance(newick_id, StringTypes):
                newick_id = ObjectId(newick_id)
            else:
                self.bind_object.set_error("newick_id必须为ObjectId对象或其对应的字符串!", code="51003503")
        with open(path, 'rb') as r:
            line = r.readline().strip('\r\n')
            samples = self.bind_object.hcluster.option("newicktree").prop['sample']
        new_samples = list()
        for sp in samples:
            new_samples.append(re.split("\s*;\s*", sp).pop(-1))
        insert_data = {
            "value": line,
            "samples": new_samples
        }
        collection = self.db['sg_newick_tree']
        try:
            # collection.find_one_and_update({"_id": newick_id}, {"$set": insert_data})
            collection.update({"_id": newick_id}, {"$set": insert_data})
        except Exception as e:
            self.bind_object.logger.error("导入newick tree%s信息出错:%s" % (path, e))
        else:
            self.bind_object.logger.info("导入newick tree%s信息成功!" % path)
