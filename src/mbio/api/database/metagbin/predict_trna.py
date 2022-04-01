# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json, re, os
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId


class PredictTrna(Base):
    def __init__(self, bind_object):
        super(PredictTrna, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_predict_trna(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "tRNA预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "file_path": "",
            "amino_list": "",
            "genome_id": genome_id,
            "name": name if name else "tRNAPredict_Origin",
        }
        collection = self.db["predict_trna"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_predict_trna_detail(self, inserted_id, genome_id, predict_gff=None):
        gff = pd.read_table(predict_gff, sep='\t', header=0)
        if len(gff) < 1:
            return
        ttype = gff["tRNA Type"]
        num = len(ttype)
        type_num = len(set(ttype))
        amino = sorted(dict(Counter(ttype)).items(),key=lambda items:items[0])
        gff["location"] = gff["Sequence id"].str.split("_", expand=True)[0]
        gff["gff_id"] = gff["Sequence id"].str.split("_", expand=True)[1]
        data_list = []
        for i in range(len(gff)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "genome_id": genome_id,
                "gene_id": gff["Gene ID"][i],
                "location": gff["location"][i],
                "trna": gff["gff_id"][i],
                "start": gff["Start"][i],
                "end": gff["End"][i],
                "type": gff["tRNA Type"][i],
                "anti": gff["Anti Codon"][i],
                "score": gff["Score"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        data = {
            "predict_id": ObjectId(inserted_id),
            "genome_id": genome_id,
            "num": num,
            "type_num": type_num,
            "amino": amino
        }
        try:
            collection = self.db["predict_trna_detail"]
            collection.insert_many(data_list)
            collection2 = self.db["predict_trna_stat"]
            collection2.insert_one(data)
            main_collection = self.db["predict_trna"]
            amino_list = sorted(list(set(
                list(set(ttype)) + main_collection.find_one({"_id": ObjectId(inserted_id)})['amino_list'].split(','))))
            if '' in amino_list:
                amino_list.remove('')
            mtype = ','.join(amino_list)
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"amino_list": mtype,
                                                                             'main_id':ObjectId(inserted_id),
                                                                             "file_path":predict_gff,
                                                                             }})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

    """

    def add_trna_predict_seq(self, inserted_id, specimen_id, fnn, struc):
        if not os.path.exists(fnn):
            self.bind_object.logger.info("不存在该文件！")
            return
        data_list = []
        gene_list = []
        fnnseq = {}
        strucseq = []
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
        with open(struc, "r") as lines:
            for l1 in lines:
                l2 = lines.next()
                l3 = lines.next()
                l4 = lines.next()
                l5 = lines.next()
                l6 = lines.next()
                if re.match("^\S+", l3):
                    l6 = l6 + lines.next()
                struc_one = l1 + l2 + l3 + l4 + l5 + l6
                strucseq.append(struc_one)
        for g in range(len(gene_list)):
            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id": gene_list[g],
                "seq_info": fnnseq[gene_list[g]][0],
                "seq": fnnseq[gene_list[g]][1],
                "struc": strucseq[g]
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["trna_predict_seq"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (fnn + "and" + struc, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % (fnn + "and" + struc))

    """