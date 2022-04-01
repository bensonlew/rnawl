# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class PredictRepeat(Base):
    def __init__(self, bind_object):
        super(PredictRepeat, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_predict_repeat(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Repeat预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "RepeatPredict_Origin",
            "genome_id": genome_id,
        }
        collection = self.db["predict_repeat"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_predict_repeat_detail(self, inserted_id, genome_id, predict_gff=None):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        num = len(gff)
        total_base = gff["Length"].sum()
        main_collection = self.db["predict_repeat"]
        task_id = main_collection.find_one({"_id": ObjectId(inserted_id)})['task_id']
        main_gene = self.db["predict_cds"]
        cds_id = main_gene.find_one({"task_id": task_id, "genome_id": genome_id})['_id']
        detail_gene = self.db["assembly_detail"]
        genome_total_base = detail_gene.find_one({"genome_id": genome_id})['scaf_base']
        # self.bind_object.logger.info(genome_total_base)
        in_genome = float(total_base) / float(genome_total_base)
        gff["gff_id"] = (gff["Attributes"].str.split(";", expand=True)[0]).str.split("=", expand=True)[1]
        gff["period_size"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split("=", expand=True)[1].astype(int)
        gff["copy_no"] = (gff["Attributes"].str.split(";", expand=True)[2]).str.split("=", expand=True)[1].astype(float)
        gff["percent_matches"] = (gff["Attributes"].str.split(";", expand=True)[3]).str.split("=", expand=True)[1].astype(int)
        gff["percent_indels"] = (gff["Attributes"].str.split(";", expand=True)[4]).str.split("=", expand=True)[1].astype(int)
        gff["consensus"] = (gff["Attributes"].str.split(";", expand=True)[5]).str.split("=", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "genome_id": genome_id,
                "location": gff["Sequence Name"][i],
                "repeat": gff["gff_id"][i],
                "strand": gff["Strand"][i],
                "start": gff["Start"][i],
                "end": gff["End"][i],
                "len": gff["Length"][i],
                "period_size": gff["period_size"][i],
                "copy_num": gff["copy_no"][i],
                "match": gff["percent_matches"][i],
                "indel": gff["percent_indels"][i],
                "consensus": gff["consensus"][i],
                "type": "Tandem Repeat",
            }
            data_son = SON(data)
            data_list.append(data_son)
        data = {
            "predict_id": ObjectId(inserted_id),
            "genome_id": genome_id,
            "num": num,
            "total_base": total_base,
            "in_genome": in_genome,
            "type": "Tandem Repeat",
        }
        try:
            collection = self.db["predict_repeat_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["predict_repeat_stat"]
            collection2.insert_one(data)
            main_collection = self.db["predict_repeat"]
            main_collection.update({'_id': ObjectId(inserted_id)},{'$set':{'main_id':ObjectId(inserted_id),
                                                                        }})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

