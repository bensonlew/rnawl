# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class RepeatPredict(Base):
    def __init__(self, bind_object):
        super(RepeatPredict, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_repeat_predict(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Repeat预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "RepeatPredict_Origin"
        }
        collection = self.db["repeat_predict"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_repeat_predict_detail(self, inserted_id, specimen_id, predict_gff):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        num = len(gff)
        total_base = gff["Length"].sum()
        main_collection = self.db["repeat_predict"]
        task_id = main_collection.find_one({"_id": ObjectId(inserted_id)})['task_id']
        main_gene = self.db["gene_predict"]

        #cds_id = main_gene.find_one({"task_id": task_id})['_id']
        pre_all = main_gene.find({"task_id": task_id})
        detail_gene = self.db["gene_predict_specimen"]
        genome_total_base = ''
        for g_p in pre_all:
            cds_id = g_p['_id']
            g_p_detail = detail_gene.find_one({"specimen_id": specimen_id, "predict_id": cds_id})
            if g_p_detail:
                genome_total_base = g_p_detail['total_base']
                break

        # self.bind_object.logger.info(genome_total_base)
        if genome_total_base not in ['','-']:
            in_genome = float(total_base) / float(genome_total_base)
        else:
            in_genome = '-'
        gff["gff_id"] = (gff["Attributes"].str.split(";", expand=True)[0]).str.split("=", expand=True)[1]
        gff["period_size"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split("=", expand=True)[1].astype(
            int)
        gff["copy_no"] = (gff["Attributes"].str.split(";", expand=True)[2]).str.split("=", expand=True)[1].astype(float)
        gff["percent_matches"] = (gff["Attributes"].str.split(";", expand=True)[3]).str.split("=", expand=True)[
            1].astype(int)
        gff["percent_indels"] = (gff["Attributes"].str.split(";", expand=True)[4]).str.split("=", expand=True)[
            1].astype(int)
        gff["consensus"] = (gff["Attributes"].str.split(";", expand=True)[5]).str.split("=", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "location": gff["Sequence Name"][i],
                "repeat": gff["gff_id"][i],
                "strand": gff["Strand"][i],
                "start": gff["Start"][i],
                "end": gff["End"][i],
                "len": gff["Length"][i],
                "period_size": gff["period_size"][i],
                "copy_no": gff["copy_no"][i],
                "percent_matches": gff["percent_matches"][i],
                "percent_indels": gff["percent_indels"][i],
                "consensus": gff["consensus"][i],
                "type": "TandemRepeat",
            }
            data_son = SON(data)
            data_list.append(data_son)
        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "num": num,
            "total_base": total_base,
            "in_genome": in_genome,
            "type": "TandemRepeat",
        }
        try:
            collection = self.db["repeat_predict_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["repeat_predict_specimen"]
            collection2.insert_one(data)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)
