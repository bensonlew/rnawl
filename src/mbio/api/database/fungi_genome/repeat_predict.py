# -*- coding: utf-8 -*-
# __author__ = 'zhujuan zouguanqing'
# last modify: 2018.06.12

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId
import re


class RepeatPredict(Base):
    def __init__(self, bind_object):
        super(RepeatPredict, self).__init__(bind_object)
        self._project_type = "fungigenome"

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
        return inserted_id

    @report_check
    def add_repeat_predict_detail(self, inserted_id, specimen_id, predict_gff, predict_tbl):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        #num = len(gff)
        #total_base = gff["Length"].sum()
        main_collection = self.db["repeat_predict"]
        task_id = main_collection.find_one({"_id": ObjectId(inserted_id)})['task_id']
        #main_gene = self.db["gene_predict"]
        #cds_id = main_gene.find_one({"task_id": task_id})['_id']
        detail_gene = self.db["gene_predict_specimen"]
        #genome_total_base = detail_gene.find_one({"specimen_id": specimen_id, "predict_id": ObjectId(cds_id)})[
        #    'total_base']
        # self.bind_object.logger.info(genome_total_base)
        #in_genome = float(total_base) / float(genome_total_base)

        gff["gff_id"] = (gff["Attributes"].str.split(";", expand=True)[0]).str.split("=", expand=True)[1]
        gff["class"] = (gff["Attributes"].str.split(";",expand=True)[2]).str.split("=", expand=True)[1]
        gff["match"] = ((gff["Attributes"].str.split(";",expand=True)[1]).str.split(" ",expand=True)[0]).str.split("=",expand=True)[1]
        gff["r_start"] = (gff["Attributes"].str.split(";",expand=True)[1]).str.split(" ",expand=True)[1]
        gff["r_end"] = (gff["Attributes"].str.split(";",expand=True)[1]).str.split(" ",expand=True)[2]
        #gff["period_size"] = (gff["Attributes"].str.split(";", expand=True)[1]).str.split("=", expand=True)[1].astype(
        #    int)
        #gff["copy_no"] = (gff["Attributes"].str.split(";", expand=True)[2]).str.split("=", expand=True)[1].astype(float)
        #gff["percent_matches"] = (gff["Attributes"].str.split(";", expand=True)[3]).str.split("=", expand=True)[
        #   1].astype(int)
        #gff["percent_indels"] = (gff["Attributes"].str.split(";", expand=True)[4]).str.split("=", expand=True)[
        #    1].astype(int)
        #gff["consensus"] = (gff["Attributes"].str.split(";", expand=True)[5]).str.split("=", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "location": gff["Sequence Name"][i],
                "repeat": gff["gff_id"][i],
                "strand": gff["Strand"][i],
                "start": gff["Begin"][i],
                "end": gff["End"][i],
                #"len": gff["Length"][i],
                "len": int(gff["End"][i])-int(gff["Begin"][i])+1,
                #"period_size": gff["period_size"][i],
                #"copy_no": gff["copy_no"][i],
                #"percent_matches": gff["percent_matches"][i],
                #"percent_indels": gff["percent_indels"][i],
                #"consensus": gff["consensus"][i],
                "type": gff["class"][i],
                "match": gff["match"][i],
                "r_start": int(gff["r_start"][i]),
                "r_end" : int(gff["r_end"][i]),
                "r_len" : int(gff["r_end"][i])-int(gff["r_start"][i])+1

            }
            data_son = SON(data)
            data_list.append(data_son)

        with open(predict_tbl) as fr:
            for line in fr:
                if line.find("bases masked:") != -1 :
                    bases_masked = re.findall("bases\smasked:\s+(\d+)\s*bp\D+([\d|\.]*)",line)
                elif line.find("SINEs:") != -1 :
                    sines = re.findall("SINEs:\D+(\d+)\D+",line)
                elif line.find("LINEs:")!= -1 :
                    lines = re.findall("LINEs:\D+(\d+)\D+",line)
                elif line.find("LTR elements:")!= -1 :
                    ltr = re.findall("LTR\selements:\D+(\d+)\D+",line)
                elif line.find("DNA elements:")!= -1 :
                    dna = re.findall("DNA\selements:\D+(\d+)\D+",line)
                elif line.find("Unclassified:")!= -1 :
                    unclass = re.findall("Unclassified:\D+(\d+)\D+",line)
                elif line.find("Small RNA") !=-1:
                    srna = re.findall("Small RNA:\D+(\d+)\D+",line)
                elif line.find("Satellites:")!= -1 :
                    satellite = re.findall("Satellites:\D+(\d+)\D+",line)
                elif line.find("Simple repeats:")!= -1 :
                    simple_rep = re.findall("Simple repeats:\D+(\d+)\D+",line)
                elif line.find("Low complexity:")!= -1 :
                    low_comp = re.findall("Low complexity:\D+(\d+)\D+",line)




        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "total_base": int(bases_masked[0][0]),
            "in_genome": float(bases_masked[0][1]),
            "type": "Interspersed",
            "sine": sines[0] ,
            "line": lines[0] ,
            "ltr": ltr[0] ,
            "dna": dna[0],
            "unclass": unclass[0],
            "srna" :srna[0],
            "satellite": satellite[0],
            "simple_rep":simple_rep[0],
            "low_comp":low_comp[0]

        }

        try:
            collection = self.db["repeat_predict_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["repeat_predict_specimen"]
            collection2.insert_one(data)
            main_collection = self.db["repeat_predict"]
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})

        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)
