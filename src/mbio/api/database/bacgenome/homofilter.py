# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last modify: 2019.04.20

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId
import types


class Homofilter(Base):
    def __init__(self, bind_object):
        super(Homofilter, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_homofilter(self, params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "核心和特有基因分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Homofilter_Origin"
        }
        collection = self.db["compare_homofilter"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_homofilter_detail(self, inserted_id, detail_file):
        inserted_id = self.check_id(inserted_id)
        data_list = []
        detail_table = pd.read_table(detail_file, sep='\t', header=0)
        if len(detail_table) < 1:
            return
        for i in range(len(detail_table)):
            data = {
                "homofilter_id": inserted_id,
                "cluster_id": detail_table["ClusterID"][i],
                "gene_id": detail_table["Gene ID"][i],
                "sample": str(detail_table["Sample Name"][i]),
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["compare_homofilter_detail"]
            collection.insert_many(data_list)
            main_collection = self.db['compare_homofilter']
            main_collection.update({'_id':inserted_id},{"$set":{"main_id": inserted_id}})
        except Exception, e:
            self.bind_object.set_error("导入%s结果表出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % detail_file)


    @report_check
    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！')
        return object_id



