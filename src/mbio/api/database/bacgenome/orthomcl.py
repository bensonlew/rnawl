# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last modify: 2019.04.15

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId
import types
import re


class Orthomcl(Base):
    def __init__(self, bind_object):
        super(Orthomcl, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_homology(self, params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "同源分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Orthomcl_Origin"
        }
        collection = self.db["compare_homology"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_homology_detail(self, inserted_id, samples, stat_file, up_path):
        inserted_id = self.check_id(inserted_id)
        data_list = []
        stat = pd.read_table(stat_file, sep='\t', header=0)
        samples = samples.split(",")
        self.update_stat_file(inserted_id, up_path)
        if len(stat) < 1:
            return
        for i in range(len(stat)):
            data = {
                "homology_id": inserted_id,
                "cluster_id": stat["#ClusterID"][i],
            }
            for j in samples:
                genes = stat[j][i]
                data[j] = genes
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["compare_homology_detail"]
            collection.insert_many(data_list)
            main_collection = self.db['compare_homology']
            main_collection.update({'_id':inserted_id},{"$set":{"main_id": inserted_id}})

        except Exception, e:
            self.bind_object.set_error("导入%s结果表出错:%s" % (stat_file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % stat_file)

    @report_check
    def add_homology_venn(self, inserted_id, venn_file):
        inserted_id = self.check_id(inserted_id)
        data_list = []
        detail = pd.read_table(venn_file, sep='\t', header=0)
        if len(detail) < 1:
            return
        for i in range(len(detail)):
            lable = re.sub('\s*only$','',detail["#Lables"][i])
            lable = '&'.join(sorted([ spl.strip() for spl in lable.split('&')]))
            data = {
                "homology_id": inserted_id,
                "lable": lable,
                "family_number":  detail["Family_number"][i],
                "gene_number": detail["Gene_number"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["compare_homology_venn"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s结果表出错:%s" % (venn_file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % venn_file)

    @report_check
    def update_stat_file(self, main_table_id, stat_file):
        self.db['compare_homology'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'stat_file': stat_file,'main_id':ObjectId(main_table_id)}})

    @report_check
    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！')
        return object_id



