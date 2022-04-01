# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20170926
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
import re
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id


class HclusterTree(Base):
    def __init__(self, bind_object):
        super(HclusterTree, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_hcluster_tree(self, file_path, main=False, tree_id=None, task_id=None, anno_type=None, params=None,
                          update_dist_id=None, name=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "Hcluster_Origin",
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'anno_type': anno_type
            }
            collection = self.db['hcluster_tree']
            # 将主表名称写在这里
            tree_id = collection.insert_one(insert_data).inserted_id
        else:
            if tree_id is None:
                self.bind_object.set_error("main为False时需提供tree_id!", code="52801201")
            if not isinstance(tree_id, ObjectId):
                tree_id = ObjectId(tree_id)
        if update_dist_id:
            self.update_dist(update_dist_id, tree_id)
        collection = self.db['hcluster_tree']
        tree_info= collection.find_one({'_id': tree_id})
        task_name2id = tree_info["task_id"]
        self.sample_2_id = name2id(task_name2id, type="task")
        self.bind_object.logger.info(self.sample_2_id)
        with open(file_path, 'r') as f:
            specimen_tree = f.readline()
            specimen = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', specimen_tree)
            for i in range(len(specimen)):
                sp = specimen[i][1]
                sp_id = self.sample_2_id[sp]
                specimen_tree = specimen_tree.replace("(" + sp + ":", "(" + str(sp_id) + ":")
                specimen_tree = specimen_tree.replace("," + sp + ":", "," + str(sp_id) + ":")
            try:
                collection = self.db["hcluster_tree"]
                collection.update_one({"_id": ObjectId(tree_id)}, {"$set": {"specimen_tree": specimen_tree}})
            except Exception as e:
                self.bind_object.logger.error("导入hcluster tree%s信息出错:%s" % (file_path, e))
                self.bind_object.set_error("导入hcluster tree信息出错", code="52801202")
            else:
                self.bind_object.logger.info("导入hcluster tree%s信息成功!" % file_path)
        return tree_id

    @report_check
    def update_dist(self, distance_id, tree_id):
        """
        从newick树更新距离矩阵结果的主表的newick_tree_id
        """
        self.db['specimen_distance'].update_one({'_id': ObjectId(distance_id)},
                                                {'$set': {'hcluster_tree_id': ObjectId(tree_id)}})
