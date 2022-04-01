# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181227


from biocluster.api.database.base import Base, report_check
import datetime
import json, os
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId


class PredictRrna(Base):
    def __init__(self, bind_object):
        super(PredictRrna, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_predict_rrna(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "rRNA预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "rRNAPredict_Origin",
            "genome_id": genome_id,
        }
        collection = self.db["predict_rrna"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_predict_rrna_detail(self, inserted_id, genome_id, predict_gff=None):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        rtype = gff["Attributes"].str.split(";", expand=True)[0]
        num = len(rtype)
        rrna = Counter(rtype)
        s16 = rrna["Name=16S_rRNA"]
        s23 = rrna["Name=23S_rRNA"]
        s5 = rrna["Name=5S_rRNA"]
        gff["type"] = rtype.str.split("=", expand=True)[1]
        gff["location"] = gff["Sequence id"].str.split("_", expand=True)[0]
        gff["gff_id"] = gff["Sequence id"].str.split("_", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "genome_id": genome_id,
                "gene_id": gff["Gene ID"][i],
                "location": gff["location"][i],
                "rrna": gff["gff_id"][i],
                "strand": gff["Strand"][i],
                "start": gff["Start"][i],
                "end": gff["End"][i],
                "type": gff["type"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        data = {
            "predict_id": ObjectId(inserted_id),
            "genome_id": genome_id,
            "num": num,
            "16s": s16,
            "23s": s23,
            "5s": s5
        }
        try:
            collection = self.db["predict_rrna_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["predict_rrna_stat"]
            collection2.insert_one(data)

        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

    @report_check
    def add_predict_rrna_seq(self, inserted_id, genome_id, fnn=None, predict_gff=None, rrna_path=None):
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
                line = line.strip('\r\n')
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
                "genome_id": genome_id,
                "gene_id": gene,
                "seq_info": fnnseq[gene][0],
                "seq_fnn": fnnseq[gene][1],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["predict_rrna_seq"]
            collection.insert_many(data_list)
            main_collection = self.db["predict_rrna"]
            if predict_gff != None:
                main_collection.update({'_id': ObjectId(inserted_id)},{'$set':{'main_id':ObjectId(inserted_id),
                                                                                "fnn_path": predict_gff}})
            else:
                main_collection.update({'_id': ObjectId(inserted_id)},{'$set':{'main_id':ObjectId(inserted_id),
                                                                               }})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (fnn, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % fnn)
