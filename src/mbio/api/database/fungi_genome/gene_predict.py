# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class GenePredict(Base):
    def __init__(self, bind_object):
        super(GenePredict, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_gene_predict(self, file_path1, file_path2, file_path3, params=None):
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
            "file_path": [file_path1, file_path2, file_path3],
            "name": "CDSpredict_Origin",
            "suffix_nu": []
        }
        collection = self.db["gene_predict"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_gene_predict_specimen(self, inserted_id, specimen_id, sample_stat):
        data_list = []
        with open(sample_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "predict_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "num": line[1],
                    "total_base": line[2],
                    "aver_base": line[3],
                    "density": line[4],
                    "gc_ratio": line[5],
                    "gene_ratio": line[6],
                    "gene_base": line[7],
                    "gc_rate": line[8],
                    "inter_ratio": line[9],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["gene_predict_specimen"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (sample_stat, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % sample_stat)

    @report_check
    def add_gene_predict_detail(self, inserted_id, specimen_id, predict_gff):
        data_list = []
        suffix_nu = {}
        with open(predict_gff, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                scf = line[1].split("_")[0].capitalize()
                data = {
                    "predict_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "gene_id": line[0],
                    "location": scf,
                    "cds_id": line[1].split("_")[1],
                    "start": int(line[2]),
                    "end": int(line[3]),
                    "strand": line[4],
                    "nul_len": int(line[5]),
                    "prot_len": int(line[6]),
                    "init": line[9],
                    "termina": line[10],
                    "exon_num": line[11],
                    "exon" : line[12],
                    "intron": line[13]
                }
                if line[4] =='-'  and data['start'] < data['end']:
                    data['start'], data['end'] = data['end'], data['start']
                data_son = SON(data)
                data_list.append(data_son)
                gene_prefix = re.subn("[0-9]*$", "", line[0])[0]
                if gene_prefix not in suffix_nu.keys():
                    gene_suffix = re.subn(gene_prefix, "", line[0])[0]
                    suffix_nu[gene_prefix] = len(gene_suffix)
        try:
            collection = self.db["gene_predict_detail"]
            collection.insert_many(data_list)

            main_collection = self.db['gene_predict']
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})

        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            tmp = {}
            if '.' in specimen_id:
                specimen_id = specimen_id.replace('.','____')
            tmp[specimen_id] = suffix_nu
            main_collection = self.db["gene_predict"]
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$push": {"suffix_nu": tmp}}, False, True)
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

    @report_check
    def add_gene_predict_seq(self, inserted_id, specimen_id, fnn, faa):
        data_list = []
        gene_list = []
        fnnseq = {}
        faaseq = {}
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
        with open(faa, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(">"):
                    line = line.strip()
                    gene_id = (line.split(" ")[0]).split(">")[1]
                    faaseq[gene_id] = ""
                    # seq_info = line
                    # faaseq[gene_id].append(seq_info)
                else:
                    seq_fnn = line
                    faaseq[gene_id] = faaseq[gene_id] +  seq_fnn
        for gene in gene_list:
            sp_info = fnnseq[gene][0].split(' ')

            if sp_info[6] == '+' :
                seq_info_str = ' '.join([sp_info[0],sp_info[4],sp_info[5],sp_info[3],sp_info[1],sp_info[2]])
            else:
                seq_info_str = ' '.join([sp_info[0],sp_info[5],sp_info[4],sp_info[3],sp_info[2],sp_info[1]])

            data = {
                "predict_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id": gene,
                "seq_info": seq_info_str,
                "seq_fnn": fnnseq[gene][1],
                "seq_faa": faaseq[gene]
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["gene_predict_seq"]
            collection.insert_many(data_list)
        except Exception, e :
            self.bind_object.logger.info("导入%s结果表出错:%s" % (fnn + "and" + faa, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % (fnn + "and" + faa))

    @report_check
    def add_gene_predict_bar(self, inserted_id, specimen_id, distribute):
        step_data = {}
        with open(distribute, "r") as f:
            lines = f.readlines()
            for line in lines[1:-1]:
                line = line.strip().split('\t')
                step_data[line[0]] = int(line[1])
        data = {
            "predict_id": ObjectId(inserted_id),
            "specimen_id": specimen_id,
            "step_data": step_data,
        }
        try:
            collection = self.db["gene_predict_bar"]
            collection.insert_one(data)
        except Exception, e:
            self.bind_object.logger.error('导入%s信息出错：%s' % (distribute, e))
            self.bind_object.set_error('导入%s信息出错：%s' , variables=(distribute, e), code="52102101")
        else:
            self.bind_object.logger.info('导入%s信息成功！' % distribute)
