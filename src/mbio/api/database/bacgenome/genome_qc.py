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
import os




class GenomeQc(Base):
    def __init__(self, bind_object):
        super(GenomeQc, self).__init__(bind_object)
        self._project_type = "bacgenome"


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
            "version": "3.1"
        }
        collection = self.db["datastat"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
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
                    "software":software
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

    def add_qc_stat_complete(self,datastat_id, file_path, lib_type, sample):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "seq_type": lib_type,
                    "specimen_id": sample,
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
            #zouguanqing 20180810 >>>
            if seq_type == 'raw' :
                clean_path = re.sub('raw_fastxstat$', 'clean_fastxstat', file_path)
                if not os.path.exists(clean_path):
                    raise Exception,"没有{}文件".format(clean_path)
                clean_file = open(clean_path)
                num = len(clean_file.readlines())
                lines = lines[:num]
            ###<<<
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
                data = {
                    "datastat_id": ObjectId(datastat_id),
                    "analysis_type": analysis_type,
                    "data_type": data_type,
                    "specimen_id": sample_name,
                    "categories": line[0],
                    "data": line[1],
                    "line_data": line[2]
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
    def add_datastat_specimen(self, datastat_id, file_path, genome_path=None, gff_path=None):
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
                if genome_path:
                    data["genome_path"] = genome_path
                if gff_path:
                    data["gff_path"] = gff_path
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["datastat_specimen"]
            collection.insert_many(data_list)
            main_collection = self.db["datastat"]
            main_collection.update({"_id": ObjectId(datastat_id)},
                                   {"$set": {"status": 'end'}})
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
                if re.search(r'^gene', line[3]):
                    seq_type = 'Chromosome'
                elif re.search(r'_gene', line[3]):
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
                    #"gene_kegg": line[7],  # 细菌升级v2 版去掉 zouguanqing 20190402
                    #"gene_cog": line[8]
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
                line = line.rstrip('\n').split("\t")
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
                    #"gene_kegg": line[8],  # 细菌升级v2 版去掉 zouguanqing 20190402
                    #"gene_cog": line[9]
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
            ##guanqing.zou 20180810 >>>
            pat=re.compile('[^(,]*\|[^:,)]*')
            g=pat.findall(line)
            for i in g:
                spi=i.split('|',1)
                new=spi[1]+'|'+spi[0]
                line=re.sub(i.replace('|','\|'),new,line)
            ##<<<

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

    def add_datastat_blast(self,blast_file, datastat_id,type):   #zouguanqing 20190327
        if not isinstance(datastat_id, ObjectId):
            datastat_id = ObjectId(datastat_id)
        insert_data = []
        if type in ['hgene']:
            head = ['sample', 'hit', 'identity', 'coverage', 'evalue', 'score']
            head_index = [0, 1, 2, 14, 10, 11]
            index = dict(zip(head, head_index))
            with open(blast_file) as f:
                for line in f:
                    tmp = {}
                    sp = line.strip().split('\t')
                    for k in index.keys():
                        if k == 'coverage':
                            tmp[k] = round(float(sp[index[k]]) * 100, 2)
                        elif k == 'hit':
                            dd = sp[index[k]].split("|")
                            if type == 'hgene':
                                tmp['hit'] = dd[0]
                                tmp['hit_spe'] = dd[1]
                            else:
                                tmp['hit'] = dd[0].split("_rRNA")[0]
                                tmp['hit_spe'] = dd[1]
                        else:
                            tmp[k] = sp[index[k]]
                    tmp['datastat_id'] = datastat_id
                    tmp['type'] = type
                    insert_data.append(tmp)
        elif type in ['16s']:
            with open(blast_file) as f:
                for line in f:
                    tmp = {}
                    sp = line.strip().split('\t')
                    num =abs(int(sp[7])-int(sp[6])+1)/float(sp[3])*100
                    tmp['sample'] = sp[0]
                    tmp['hit'] = sp[1].split("|")[0].split("_rRNA")[0]
                    tmp['hit_spe'] = sp[1].split("|")[1]
                    tmp['identity'] = sp[2]
                    tmp['coverage'] = num
                    tmp['evalue'] = sp[10]
                    tmp['score'] = sp[11]
                    tmp['type'] = type
                    tmp['datastat_id'] = datastat_id
                    insert_data.append(tmp)
        try:
            collection = self.db["datastat_blast"]
            collection.insert_many(insert_data)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (blast_file, e))
        else:
             self.bind_object.logger.info("导入%s结果表成功!" % blast_file)

    def add_gene(self,datastat_id, sample, gene_pre=None, map=None):
        if not isinstance(datastat_id,ObjectId):
            datastat_id = ObjectId(datastat_id)
        insert_list = []
        if map:
            for k in map.keys():
                genome_type = re.sub('\d*$','',map[k])
                insert = {
                    "datastat_id":datastat_id,
                    "specimen_id" : sample,
                    "origin_prefix" : gene_pre[sample][map[k]],
                    "gene_prefix" :  gene_pre[sample][map[k]],
                    "genome_type" : genome_type,
                    "seq_type" : map[k]
                }
                insert_list.append(insert)
        else:
            insert = {
                "datastat_id":datastat_id,
                "specimen_id" : sample,
                "origin_prefix" : gene_pre[sample],
                "gene_prefix" :  gene_pre[sample]
            }
            insert_list.append(insert)
        try:
            collection = self.db['datastat_gene']
            collection.insert_many(insert_list)

        except Exception, e:
            self.bind_object.logger.info('导入datastat_gene 信息失败')

    def add_sample(self,datastat_id,sample):
        if not isinstance(datastat_id,ObjectId):
            datastat_id = ObjectId(datastat_id)
        insert = {
            "datastat_id":datastat_id,
            "anay_name" : sample,
            "init_name" : sample,
        }
        try:
            collection = self.db['datastat_specimen']
            collection.insert_one(insert)
        except Exception, e:
            self.bind_object.logger.info('datastat_specimen 信息失败')
                    

