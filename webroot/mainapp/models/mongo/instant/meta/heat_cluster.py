# -*- coding: utf-8 -*-
# __author__ = 'xuting'

from bson.objectid import ObjectId
import datetime
import re
from types import StringTypes
from mainapp.models.mongo.core.base import Base
from mbio.files.meta.beta_diversity.newick_tree import NewickTreeFile
# from biocluster.config import Config


class HeatClusterMongo(Base):
    def __init__(self, bind_object):
        super(HeatClusterMongo, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self._params = self.PackParams()

    def PackParams(self):
        data = self.bind_object.data
        params = dict()
        params['otu_id'] = str(data.otu_id)
        params['level_id'] = str(data.level_id)
        params['group_id'] = str(data.group_id)
        my_sp = re.split(',', data.specimen_ids)
        my_sp.sort()
        params['specimen_ids'] = ','.join(my_sp)
        params["submit_location"] = data.submit_location
        params["task_type"] = data.task_type
        params = self.SortDict(params)
        return params

    def CreateTreeTable(self, filePath, name=None):
        data = self.bind_object.data
        if data.otu_id != 0 and not isinstance(data.otu_id, ObjectId):
            if isinstance(data.otu_id, StringTypes):
                data.otu_id = ObjectId(data.otu_id)
            else:
                raise Exception("传入的otu_id必须为ObjectId对象或其对应的字符串")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": data.otu_id})
        if not result:
            raise Exception("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(data.otu_id)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        T = NewickTreeFile()
        T.set_path(filePath)
        T.get_info()
        samples = T.prop['sample']
        with open(filePath, 'rb') as r:
            line = r.readline().strip('\r\n')
        insertData = {
            "project_sn": project_sn,
            "task_id": task_id,
            "table_type": "dist",
            "tree_type": "cluster",
            "value": line,
            "samples": samples,
            "hcluster_method": data.linkage,
            "status": "end",
            "desc": "热图物种聚类树",
            "name": name if name else "cluster_newick",
            "params": self._params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db['sg_newick_tree']
        insertedId = collection.insert_one(insertData).inserted_id
        return insertedId
