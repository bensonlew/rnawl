# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last modify: 2019.04.15

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class Methylation(Base):
    def __init__(self, bind_object):
        super(Methylation, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_methy(self, params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "甲基化分析",
            "params":  "{soft:smrtanalysis}" if not params else json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Methylation_Origin"
        }
        collection = self.db["methylation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_methy_stat(self, inserted_id, sample, stat_file):
        data_list = []
        stat = pd.read_csv(stat_file)
        if len(stat) < 1:
            return
        for i in range(len(stat)):
            data = {
                "methy_id": inserted_id,
                "specimen_id": sample,
                "motif": stat["motifString"][i],
                "center": stat["centerPos"][i],
                "modif_type": stat["modificationType"][i],
                "fraction":  float(stat["fraction"][i]),
                "n_detected": stat["nDetected"][i],
                "n_genome": stat["nGenome"][i],
                "group_tag": stat["groupTag"][i],
                "partner_motif": stat["partnerMotifString"][i],
                "mean_score": stat["meanScore"][i],
                "mean_ipd_ratio": stat["meanIpdRatio"][i],
                "mean_coverage": stat["meanCoverage"][i],
                "objective_score": stat["objectiveScore"][i],

            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["methylation_stat"]
            collection.insert_many(data_list)
            main_col = self.db['methylation']
            main_col.update({"_id":inserted_id},{"$set":{"main_id": inserted_id}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (stat_file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % stat_file)

    @report_check
    def add_methy_detail(self, inserted_id,sample, detail_file):
        data_list = []
        detail = pd.read_table(detail_file, sep='\t', header=0)
        if len(detail) < 1:
            return
        for i in range(len(detail)):
            data = {
                "methy_id": inserted_id,
                "specimen_id": sample,
                "location": detail["Location"][i],
                "strand":  detail["Strand"][i],
                "start": detail["Start"][i],
                "end": detail["End"][i],
                "motif": detail["motifString"][i],
                "score": detail["Score"][i],
                "modif_type": detail["modificationType"][i],
                "coverage": detail["Coverage"][i],
                "ipd_ratio": detail["IPDRatio"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["methylation_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (detail_file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % detail_file)

