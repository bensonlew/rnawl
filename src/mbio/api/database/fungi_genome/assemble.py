# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.5.29

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Assemble(Base):
    def __init__(self, bind_object):
        super(Assemble, self).__init__(bind_object)
        self._project_type = "fungigenome"


    @report_check
    def add_assemble(self,analysis_type,seq_type,data_type,desc,ass_tool,kmer):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "analysis_type": analysis_type,
            "status": "faile",
            "desc": desc,
            "seq_type":seq_type,
            "data_type":data_type,
            "ass_tool": ass_tool,
            "kmer": kmer,
            "params":json.dumps({"type":"assemble"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["assemble"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_assemble_stat_uncomplete(self, assemble_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": assemble_id,
                    "specimen_id": specimen_id,
                    "scaf_no":line[0],
                    "scaf_base": line[1],
                    "scaf_no_large": line[2],
                    "scaf_base_large": line[3],
                    "scaf_len_max": line[4],
                    "scaf_n50": line[5],
                    "scaf_n90": line[6],
                    "gc_rate": line[7],
                    "n_rate": line[8],
                    "contig_no": line[9],
                    "contig_base": line[10],
                    "contig_no_large": line[11],
                    "contig_base_large": line[12],
                    "contig_len_max": line[13],
                    "contig_n50": line[14],
                    "contig_n90": line[15]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(assemble_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    def add_assemble_seq(self,assemble_id, file_path,specimen_id,type,seq_path):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[0:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": assemble_id,
                    "specimen_id": specimen_id,
                    "type": type,
                    "seq_id": line[0],
                    "len": float(line[1]),
                    "gc_rate": float(line[2]),
                    "base_stat": line[3],
                    "seq_path": seq_path
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_seq"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    @report_check
    def add_assemble_graphic(self, assemble_id, file_path,specimen_id,type):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": ObjectId(assemble_id),
                    "type": type,
                    "specimen_id": specimen_id,
                    "len": line[0],
                    "num": line[1],
                    "prop_cum": line[2],
                    "sort":int(line[3])
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_bar"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_assemble_assess_busco(self, assemble_id, file_path, specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": ObjectId(assemble_id),
                    "method": 'busco',
                    "specimen_id": specimen_id,
                    "comeplete": line[0],
                    "complete_dup": line[1],
                    "fragmented": line[2],
                    "missing": line[3]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_assess"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_assemble_assess_cegma(self, assemble_id, file_path, specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": ObjectId(assemble_id),
                    "method": 'cegma',
                    "specimen_id": specimen_id,
                    "comeplete": line[0],
                    "partial": line[1],
                    "total_orthol": line[2],
                    "average_orth": line[3]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_assess"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)