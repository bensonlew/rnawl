# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181214


from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class PredictCds(Base):
    def __init__(self, bind_object):
        super(PredictCds, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_predict_cds(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "基因预测主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "CDSpredict_Origin",
            "suffix_nu": [],
            "genome_id": genome_id,
        }
        collection = self.db["predict_cds"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_predict_cds_stat(self, inserted_id, genome_id, sample_stat):
        data_list = []
        with open(sample_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip('\r\n').split("\t")
                data = {
                    "predict_id": ObjectId(inserted_id),
                    "genome_id": genome_id,
                    "num": line[1],
                    "total_base": line[2],
                    "ave_base": line[3],
                    "density": line[4],
                    "gc_ratio": line[5],
                    "gene_ratio": line[6],
                    "inter_base": line[7],
                    "inter_gc_ratio": line[8],
                    "inter_ratio": line[9],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["predict_cds_stat"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (sample_stat, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % sample_stat)

    def add_predict_cds_detail(self, inserted_id, genome_id, predict_gff=None, fnn_path=None, faa_path=None, gff_path=None):
        data_list = []
        suffix_nu = {}
        with open(predict_gff, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip('\r\n').split("\t")
                data = {
                    "predict_id": ObjectId(inserted_id),
                    "genome_id": genome_id,
                    "gene_id": line[0],
                    "location": line[1].split("_")[0],
                    "cds_id": line[1].split("_")[1],
                    "start": int(line[2]),
                    "end": int(line[3]),
                    "strand": line[4],
                    "gene_len": int(line[5]),
                    "pro_len": int(line[6]),
                    "init": line[9],
                    "termina": line[10],
                }
                data_son = SON(data)
                data_list.append(data_son)
                gene_prefix = re.subn("[0-9]*$", "", line[0])[0]
                if gene_prefix not in suffix_nu.keys():
                    gene_suffix = re.subn(gene_prefix, "", line[0])[0]
                    suffix_nu[gene_prefix] = len(gene_suffix)
        try:
            collection = self.db["predict_cds_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (predict_gff, e))
        else:
            tmp = {}
            tmp[genome_id] = suffix_nu
            main_collection = self.db["predict_cds"]
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$push": {"suffix_nu": tmp}}, False, True)
            main_collection.update({'_id': ObjectId(inserted_id)},{'$set':{'main_id':ObjectId(inserted_id),
                                                                            "fnn_path": fnn_path,
                                                                            "faa_path": faa_path,
                                                                            "gff_path": gff_path}})
            self.bind_object.logger.info("导入%s结果表成功!" % predict_gff)

    @report_check
    def add_predict_cds_bar(self, inserted_id, genome_id, distribute):
        step_data = {}
        with open(distribute, "r") as f:
            lines = f.readlines()
            for line in lines[1:-1]:
                line = line.strip('\r\n').split('\t')
                step_data[line[0]] = int(line[1])
        data = {
            "predict_id": ObjectId(inserted_id),
            "genome_id": genome_id,
            "step_data": step_data,
        }
        try:
            collection = self.db["predict_cds_bar"]
            collection.insert_one(data)
        except Exception, e:
            self.bind_object.logger.error('导入%s信息出错：%s' % (distribute, e))
            self.bind_object.set_error('导入%s信息出错：%s' , variables=(distribute, e), code="51402701")
        else:
            self.bind_object.logger.info('导入%s信息成功！' % distribute)
