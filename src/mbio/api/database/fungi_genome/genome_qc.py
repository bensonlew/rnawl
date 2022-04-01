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


class GenomeQc(Base):
    def __init__(self, bind_object):
        super(GenomeQc, self).__init__(bind_object)
        self._project_type = "fungigenome"


    @report_check
    def add_datastat(self,analysis_type,desc,seq_type,data_type,qc_tool,min_len,min_qual):
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
            "qc_tool":qc_tool,
            "min_len": min_len,
            "min_qual": min_qual,
            "params":json.dumps({"type":"datastat_specimen"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["datastat"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_qc_stat_uncomplete(self, datastat_id, file_path,software):
        data_list = []
        if software == "fastp":
            software = "Fastp_v0.20.0"
        else:
            software = "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar,smrtanalysis"
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if int(line[1]) < 1000:
                    lib_type='PE'
                else:
                    lib_type = 'MP'
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "seq_type":lib_type,
                    "specimen_id": line[0],
                    "insert_size": line[1],
                    "read_len": line[2],
                    "raw_num": line[3],
                    "raw_base": line[4],
                    "raw_q20": line[5],
                    "raw_q30": line[6],
                    "clean_pair_num": line[7],
                    "clean_single_num": line[8],
                    "clean_base": line[9],
                    "clean_q20": line[10],
                    "clean_q30": line[11],
                    "software": software
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_qc"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_qc_stat_complete(self,datastat_id, file_path,lib_type,sample_name):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "seq_type": lib_type,
                    "specimen_id":sample_name ,
                    "clean_num_pac": line[0],
                    "clean_base": line[1],
                    "clean_max_pac": line[2],
                    "average_len": line[3]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_qc"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    @report_check
    def add_datastat_graphic(self, datastat_id, file_path,sample_name,seq_type,type,lib_type):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "seq_type": seq_type,
                    "type": type,
                    "specimen_id": sample_name,
                    "column": int(line[0]),
                    "min": line[2],
                    "max": line[3],
                    "q1": line[6],
                    "err_qual": 10**(float(line[5])/(-10))*100,
                    "median": line[7],
                    "q3": line[8],
                    "a": line[12],
                    "c": line[13],
                    "g": line[14],
                    "t": line[15],
                    "n": line[16],
                    "lib_type":lib_type,
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_graphic"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_pacbio_graphic(self, datastat_id, file_path,sample_name,analysis_type,data_type):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if analysis_type == '':
                    analysis_type = 'clean'
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "specimen_id": sample_name,
                    "categories": line[0],
                    "data": line[1],
                    "line_data": line[2],
                    "analysis_type":analysis_type,
                    "data_type":data_type
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_pacbio_graphic"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_specimen(self, datastat_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "init_name": line[0],
                    "anay_name": line[0],
                    "lib": line[2],
                    "raw_files": line[1]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_specimen"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(datastat_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_uncomplete_gene(self, datastat_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "specimen_id": line[0],
                    "gene_prefix": line[3],
                    "origin_prefix": line[3],
                    "assmble_file": line[1],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_gene"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_complete_gene(self, datastat_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                seq_type = ''
                if re.search(r'gene', line[1]):
                    seq_type = 'Chromosome'
                elif re.search(r'p', line[1]):
                    name = line[3].split('_')
                    se = name[0]
                    if len(se) == 1:
                        seq_type = se[0].upper() + 'lasmid'
                    else:
                        seq_type = se[0].upper() + 'lasmid' + se[1].upper()
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "specimen_id": line[0],
                    "genome_type": line[2],
                    "gene_prefix": line[3],
                    "origin_prefix": line[3],
                    "assmble_file": line[1],
                    "seq_type": seq_type,
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_gene"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_uncomplete_summary(self, datastat_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "specimen_id": line[0],
                    "g_size": line[1],
                    "scaf_no": line[2],
                    "gc_rate": line[3],
                    "cds_no": line[4],
                    "trna_no": line[6],
                    "rrna_no": line[5],
                    "gene_kegg": line[7],
                    "gene_cog": line[8]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_summary"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_complete_summary(self, datastat_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "specimen_id": line[0],
                    "g_size": line[1],
                    "chr_no": line[2],
                    "plas_no": line[3],
                    "gc_rate": line[4],
                    "cds_no": line[5],
                    "trna_no": line[7],
                    "rrna_no": line[6],
                    "gene_kegg": line[8],
                    "gene_cog": line[9]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_summary"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_datastat_tree_detail(self, datastat_id, file_path, specimen_id, marker, type):
        data_list = []
        with open(file_path, "r") as f:
            line = f.readline()
            #line = re.sub('[)]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):', "):", line)
        data = {
            "datastat_id": ObjectId(datastat_id),
            "marker": marker,
            "value": line,
            "type": type,
            "specimen_id": specimen_id,
        }
        data_son = SON(data)
        data_list.append(data_son)
        try:
            collection = self.db["datastat_tree"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)
