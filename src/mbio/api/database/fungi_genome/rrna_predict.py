# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json, os
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId


class RrnaPredict(Base):
    def __init__(self, bind_object):
        super(RrnaPredict, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_rrna_predict(self, file_path1, file_path2, file_path3, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "rRNA预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "file_path": [file_path1, file_path2, file_path3],
            "name": "rRNAPredict_Origin"
        }
        collection = self.db["rrna_predict"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_rrna_predict_detail(self, inserted_id, specimen_id, predict_gff):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        rtype = gff["Attributes"].str.split(";", expand=True)[0]
        num = len(rtype)
        rrna = Counter(rtype)
        s18 = rrna["Name=18S_rRNA"]
        s28 = rrna["Name=28S_rRNA"]
        s5 = rrna["Name=5S_rRNA"]
        s5_8 = rrna["Name=5_8S_rRNA"]
        gff["type"] = rtype.str.split("=", expand=True)[1]

        gff["location"] = gff["Sequence id"].str.split("_", expand=True)[0]
        gff["gff_id"] = gff["Sequence id"].str.split("_", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            type_tmp = gff["type"][i]
            if type_tmp == "5_8S_rRNA":
                type_tmp = "5.8S_rRNA"
            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id": gff["Gene ID"][i],
                "location": gff["location"][i],
                "rrna": gff["gff_id"][i],
                "strand": gff["Strand"][i],
                "start": gff["Start"][i],
                "end": gff["End"][i],
                "type": type_tmp
            }
            data_son = SON(data)
            data_list.append(data_son)
        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "num": num,
            "18s": s18,
            "28s": s28,
            "5s": s5,
            "5_8s":s5_8
        }
        try:
            collection = self.db["rrna_predict_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["rrna_predict_specimen"]
            collection2.insert_one(data)

            main_collection = self.db["rrna_predict"]
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})

        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

    @report_check
    def add_rrna_predict_seq(self, inserted_id, specimen_id, fnn):
        if not os.path.exists(fnn):
            self.bind_object.logger.info("不存在该文件！")
            return
        data_list = []
        gene_list = []
        fnnseq = {}
        with open(fnn, "r") as f:
            lines = f.readlines()
            gene_id = ""
            for line in lines:
                line = line.strip()
                if line.startswith(">"):
                    gene_id = (line.split(" ")[0]).split(">")[1]
                    seq_info = line
                    gene_list.append(gene_id)
                    fnnseq[gene_id] = [seq_info]
                else:
                    seq_fnn = line
                    fnnseq[gene_id].append(seq_fnn)
        for gene in gene_list:
            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id": gene,
                "seq_info": fnnseq[gene][0],
                "seq_fnn": fnnseq[gene][1],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["rrna_predict_seq"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (fnn, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % fnn)
