# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.3.5

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
        self._project_type = "bacgenome"


    @report_check
    def add_assemble(self,analysis_type,seq_type,data_type,desc, kmer, ass_tool):
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
            "software": ass_tool,
            "kmer": kmer,
            "params":json.dumps({"type":"assemble"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "version": "3.1"
        }
        collection = self.db["assemble"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_assemble_stat_uncomplete(self, assemble_id, file_path,specimen_id,file2=None):
        complete = ''
        if file2:
            with open(file2, "r") as f:
                lines = f.readlines()
                complete = lines[1].strip().split("\t")[0]
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
                    "contig_n90": line[15],
                    "completeness": complete,
                }
                if len(line) == 17:
                    data['depth'] = line[16]   #zouguanqing 20190327
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    def add_assemble_stat_complete(self,assemble_id, file_path,specimen_id,file2=None):
        complete = ''
        if file2:
            with open(file2, "r") as f:
                lines = f.readlines()
                complete =lines[1].strip().split("\t")[0]
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": assemble_id,
                    "specimen_id": specimen_id,
                    "chr_no": line[0],
                    "pla_no": line[1],
                    "g_size": line[2],
                    "gc_rate": line[3],
                    "completeness": complete,
                }
                if len(line) == 5:
                    data['depth'] = line[4]   #zouguanqing 20190327
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assemble_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assemble"]
            main_collection.update({"_id": ObjectId(assemble_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    def add_assemble_seq(self,assemble_id, file_path,specimen_id,type,seq_path, type_dict=None):
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
                    "seq_path": seq_path,
                }
                if type_dict:
                    list1 = type_dict[line[0]]
                    data["file_seqid"] = list1[0]
                    data["seq_type"] = list1[1]
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

    def add_assemble_complete_seq(self,assemble_id, file_path,specimen_id,seq_path, type_dict=None):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[0:]:
                line = line.strip().split("\t")
                list1 = type_dict[line[0]]
                data = {
                    "assemble_id": assemble_id,
                    "specimen_id": specimen_id,
                    "seq_id": line[0],
                    "type": line[0],
                    "len": float(line[1]),
                    "gc_rate": float(line[2]),
                    "base_stat": line[3],
                    "seq_path": seq_path,
                    "file_seqid": list1[0],
                    "seq_type": list1[1],
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
    def add_assemble_graphic(self, assemble_id, file_path,specimen_id,type,step):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assemble_id": ObjectId(assemble_id),
                    "type": type,
                    "step":step,
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


    def add_assemble_anno_pipline(self,main_id=None, specimen_id=None,analysis_type=None, len_file=None,):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "analysis_type": analysis_type,
                "status": "start",
                "desc": '序列长度统计',
                "data_type":'annotation',
                "params":json.dumps({"type":"annotation"}) ,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            collection = self.db["assemble"]
            inserted_id = collection.insert_one(insert_data).inserted_id
            collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
            return inserted_id
        else:
            data_list = []
            with open(len_file) as fr:
                fr.readline()
                for line in fr:
                    line = line.strip().split('\t')

                    data = {
                        "assemble_id": main_id,
                        "specimen_id": specimen_id,
                        "seq_id": line[0],
                        "type": line[0],
                        "len": float(line[1])
                    }
                    if len(line) > 2:
                        data['seq_sort_id'] = line[2]

                    data_son = SON(data)
                    data_list.append(data_son)
            try:
                collection = self.db["assemble_seq"]
                collection.insert_many(data_list)
                main_collection = self.db["assemble"]
                main_collection.update({"_id": main_id},
                                       {"$set": {"status": 'end'}})
            except Exception as e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (len_file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" %len_file)
