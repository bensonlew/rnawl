# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import pandas as pd
from Bio import SeqIO


class GenePredict(Base):
    def __init__(self, bind_object):
        super(GenePredict, self).__init__(bind_object)
        self._project_type = "bac_comparative"

    @report_check
    def add_gene_predict(self, gene_dir, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "基因预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "Gene_Predict_Origin",
            'gene_dir': gene_dir
        }
        collection = self.db["gene_predict"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_gene_cds_detail(self, inserted_id, specimen_id, cds_gff):
        """
        主要是用于cds详情表的调用
        :param inserted_id:
        :param specimen_id:
        :param cds_gff:
        :return:
        """
        with open(cds_gff, "r") as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:# 有标题
                line = line.strip().split("\t")
                data = {
                    "gene_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "gene": line[0],
                    "location": line[1].split("_ORF")[0],
                    "strand": line[4],
                    "start": int(line[2]),
                    "end": int(line[3]),
                    "gene_len": int(line[5]),
                    "prot_len": int(line[6]),
                    "init_codon": line[9],
                    "term_codon": line[10],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["gene_cds_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (cds_gff, e))

    @report_check
    def add_gene_rrna_detail(self, inserted_id, specimen_id, rrna_gff, fnn=None, fnn_path=None):
        """
        主要是用rrna详情表的调用
        如果预测结果中有16s必须要输入16s序列，否则会报错
        :param inserted_id:
        :param specimen_id:
        :param rrna_gff:
        :param fnn:
        :return:
        """
        rrna_seq ={}
        if fnn:
            for seq_record in SeqIO.parse(fnn, 'fasta'):
                seq_id = seq_record.id
                seq = seq_record.seq
                rrna_seq[seq_id] = str(seq)
        with open(rrna_gff, 'r') as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]: #去除表头
                line = line.strip().split("\t")
                location = line[1].split("_rRNA")[0]
                rrna_id = "rRNA" + str(line[1].split("_rRNA")[1])
                if line[5] == "-":
                    if int(line[2]) > int(line[3]):
                        start = int(line[2])
                        end = int(line[3])
                        gene_len = int(line[2]) - int(line[3]) + 1
                    else:
                        start = int(line[3])
                        end = int(line[2])
                        gene_len = int(line[3]) - int(line[2]) + 1
                else:
                    if int(line[3]) > int(line[2]):
                        start = int(line[2])
                        end = int(line[3])
                        gene_len = int(line[3]) - int(line[2]) +1
                    else:
                        start = int(line[3])
                        end = int(line[2])
                        gene_len = int(line[2]) - int(line[3]) +1
                if re.search(r'16S', line[7]):
                    type = "16S"
                elif re.search(r'5S', line[7]):
                    type = "5S"
                elif re.search(r'23S', line[7]):
                    type = "23S"
                else:
                    type = ""
                data = {
                    "gene_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "gene": line[0],
                    "location": location,
                    "strand": line[5],
                    "start": start,
                    "end": end,
                    "type": type,
                    "rrna_id": rrna_id,
                    "gene_len": gene_len,
                }
                if line[0] in rrna_seq.keys():
                    data["seq"] = rrna_seq[line[0]]
                if fnn_path:
                    data["s16_path"] = fnn_path
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["gene_rrna_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (rrna_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % rrna_gff)

    @report_check
    def add_gene_trna_detail(self, inserted_id, specimen_id, trna_gff, type, fnn=None):
        """
        主要是用于TRNA详情表的调用
        :param inserted_id:
        :param specimen_id:
        :param trna_gff:
        :return:
        """
        trna_seq ={}
        if fnn:
            for seq_record in SeqIO.parse(fnn, 'fasta'):
                seq_id = seq_record.id
                seq = seq_record.seq
                trna_seq[seq_id] = str(seq)
        with open(trna_gff, 'r') as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]: #去除表头
                line = line.strip().split("\t")
                location = line[1].split("_tRNA")[0]
                trna_id = "tRNA" + str(line[1].split("_tRNA")[1])
                if line[5] == "-":
                    if int(line[2]) > int(line[3]):
                        start = int(line[2])
                        end = int(line[3])
                        gene_len = int(line[2]) - int(line[3]) + 1
                    else:
                        start = int(line[3])
                        end = int(line[2])
                        gene_len = int(line[3]) - int(line[2]) + 1
                else:
                    if int(line[2]) < int(line[3]):
                        start = int(line[2])
                        end = int(line[3])
                        gene_len = int(line[3]) - int(line[2]) + 1
                    else:
                        start = int(line[3])
                        end = int(line[2])
                        gene_len = int(line[2]) - int(line[3]) + 1
                data = {
                    "gene_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "gene": line[0],
                    "location": location,
                    "start": start,
                    "end": end,
                    "type": line[5],
                    "trna_id": trna_id,
                    "gene_len": gene_len,
                    "ani_codon": line[6],
                }
                if type in ["seq", "gff"]:
                    data["score"] = line[8]
                else:
                    data["score"] = "-"
                if line[0] in trna_seq.keys():
                    data["seq"] = trna_seq[line[0]]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["gene_trna_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (trna_gff, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % trna_gff)












