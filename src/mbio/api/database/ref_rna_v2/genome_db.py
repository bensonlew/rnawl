# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import sys
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from pymongo import MongoClient


class GenomeDb(ApiBase):
    def __init__(self, bind_object):
        super(GenomeDb, self).__init__(bind_object)
        #self._db_name = Config().MONGODB + '_ref_rna'

    def add_genome_db(self, json_file):
        '''
        添加数据库物种信息
        '''
        collection = self.db["sg_genome_db"]
        data_list = []
        with open(json_file, 'r') as f:
            genome_json = json.loads(f.read())
            for key, value in genome_json.items():
                data_list.append(SON(value))

        if data_list:
            try:
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入注释统计信息失败:%s" % blast_path)

    def get_new_genome_id(self, genome_id):
        collection = self.db["sg_genome_db"]
        try:
            genome_dict = collection.find_one({"genome_id":genome_id})
        except:
            raise Exception("查找数据库失败")

        if genome_dict:
            genome_list = collection.find({})
            genome_list = sorted(genome_list, key = lambda genome:int(genome["genome_id"][2:]))
            genome_id_new = genome_list[-1]["genome_id"]
            genome_id_new = genome_id_new[0:2] + str(int(genome_id_new[2:]) + 1).zfill(4)
            print "genome id %s is exists us %s"%(genome_id, genome_id_new)
            return genome_id_new
        else:
            return genome_id

    def add_genome_db_from_v1(self, genome_id, json_file):
        '''
        根据v1json添加数据库物种信息
        '''
        collection = self.db["sg_genome_db"]
        data_list = []
        genome_id = self.get_new_genome_id(genome_id)
        with open(json_file, 'r') as f:
            genome_json = json.loads(f.read())
            for key, value in genome_json.items():
                collection = self.db["sg_genome_db"]
                data_list.append(SON(self.convert_v12v2(value, genome_id)))
                genome_id = genome_id[0:2] + str(int(genome_id[2:]) + 1).zfill(4)

        if data_list:
            try:
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入注释统计信息失败:%s" % blast_path)

    def convert_v12v2(self, single_dict, genome_id):
        '''
        根据v1json推测v2
        '''
        single_dict["genome_id"] = genome_id

        if "name" in single_dict:
            pass
        else:
            single_dict["name"] = single_dict["dna_fa"].split("/")[1]

        single_dict["organism_name"] = single_dict["name"]

        # 推测基因组版本号
        if single_dict.has_key("genome_version") and single_dict["genome_version"]:
            single_dict['assembly'] = single_dict["genome_version"]
        elif single_dict['assembly'] != "":
            pass
        elif single_dict['ensemble_release'] != "" or single_dict['accession'] != "":
            single_dict['assembly'] = single_dict['ensemble_release'] + "_" + single_dict['accession']
        else:
            single_dict['assembly'] = "unknown"

        # 推测注释版本号
        if single_dict.has_key("annot_version") and single_dict["annot_version"]:
            single_dict['annot_version'] = single_dict["annot_version"]
        elif single_dict['ensemble_release'] != "":
            single_dict['annot_version'] = single_dict['ensemble_release']
        else:
            single_dict['annot_version'] = "i-sanger"

        single_dict['anno_path_v2'] = single_dict['anno_path'] + "_v2"
        single_dict['g2t2p'] = single_dict['anno_path_v2'] + "/g2t2p"
        single_dict['organism_name'] = single_dict['name']
        single_dict['genome_id'] = genome_id
        single_dict['transcript'] = os.path.dirname(single_dict['dna_fa']) + "/transcript.fa"
        return single_dict


    def add_single_genome_db(self, genome_id, json_file, genome_version=None, annot_version=None):
        '''
        添加数据库物种信息
        '''
        collection = self.db["sg_genome_db"]
        has_genome_id = collection.find_one({"genome_id":genome_id})
        if has_genome_id:
            raise Exception("genome_id {}已存在请换一个")
        data_list = []
        with open(json_file, 'r') as f:
            genome_json = json.loads(f.read())
        genome_json["genome_id"] = genome_id

        if "name" in genome_json:
            pass
        else:
            genome_json["name"] = genome_json["dna_fa"].split("/")[1]

        genome_json["organism_name"] = genome_json["name"]

        if "g2t2p" in genome_json:
            pass
        else:
            genome_json["g2t2p"] = genome_json["anno_path_v2"] + "/g2t2p"

        if "transcript" in genome_json:
            pass
        else:
            genome_json["transcript"] = os.path.dirname(genome_json['dna_fa']) + "/transcript.fa"


        # 推测基因组版本号
        if genome_version:
            genome_json['assembly'] = genome_version
        elif genome_json['assembly'] != "":
            pass
        elif genome_json['ensemble_release'] != "" or genome_json['accession'] != "":
            genome_json['assembly'] = genome_json['ensemble_release'] + "_" + genome_json['accession']
        else:
            genome_json['assembly'] = "unknown"

        # 推测注释版本号
        if annot_version:
            genome_json['annot_version'] = annot_version
        elif genome_json['ensemble_release'] != "":
            genome_json['annot_version'] = genome_json['ensemble_release']
        else:
            genome_json['annot_version'] = "i-sanger"
        print genome_json

        try:
            collection.insert_one(SON(genome_json))
        except Exception, e:
            raise Exception("导入数据库失败")

    def get_genome_dict_by_task_id(self, task_id):
        "根据taskid获取基因组属性字典"
        if task_id:
            collection = self.db["sg_task"]
            try:
                genome_id = collection.find_one({"task_id":task_id})['genome_id']
                return self.get_genome_dict_by_genomeid(genome_id)
            except Exception:
                self.bind_object.logger.info("sg_task找不到task_id{}".format(task_id))
                return None

        else:
            self.bind_object.logger.info("task_id为空")
            return None

    def get_genome_dict_by_genomeid(self, genome_id):
        "根据genome_id获取genome属性字典"
        collection = self.db["sg_genome_db"]
        try:
            genome_dict = collection.find_one({"genome_id":genome_id})
        except Exception:
            self.bind_object.logger.info("genome_db找不到genome_id{}".format(genome_id))
        return genome_dict

    def updata_from_isanger(self):
        "根据genome_id获取genome属性字典"
        isanger_client = MongoClient("mongodb://rna:y6a3t5n1w0y7@10.100.1.10/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
        i_db = isanger_client.sanger_ref_rna_v2
        i_collection = i_db["sg_genome_db"]
        collection = self.db["sg_genome_db"]
        try:
            genome_dicts = i_collection.find()
            for genome_dict in genome_dicts:
                collection.insert_one(SON(genome_dict))
        except Exception:
            raise Exception("更新失败")

    def updata_to_isanger(self, genome_id):
        "根据genome_id 将一个物种从 tsg数据库 更新到isanger"
        collection = self.db["sg_genome_db"]
        genome_dict = self.get_genome_dict_by_genomeid(genome_id)

        isanger_client = MongoClient("mongodb://rna:y6a3t5n1w0y7@10.100.1.10/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
        i_db = isanger_client.sanger_ref_rna_v2
        i_collection = i_db["sg_genome_db"]
        collection = self.db["sg_genome_db"]
        try:
            i_collection.insert_one(SON(genome_dict))
        except Exception:
            raise Exception("更新失败")

    def updata_to_smallrna(self, genome_id):
        "根据genome_id 将一个物种从 tsg数据库 更新到small_rna"
        collection = self.db["sg_genome_db"]
        genome_dict = self.get_genome_dict_by_genomeid(genome_id)
        small_client = self._config.get_mongo_client(mtype="small_rna")
        small_db = small_client[self._config.get_mongo_dbname("small_rna")]
        small_collection = small_db["sg_genome_db"]
        try:
            small_collection.insert_one(SON(genome_dict))
        except Exception:
            raise Exception("更新失败")

    def update_status(self):
        "根据genome_id 将一个物种从 tsg数据库 更新到isanger"
        collection = self.db["sg_genome_db"]
        genome_dicts = collection.find()
        for genome_dict in genome_dicts:
            # print genome_dict
            insert_dict={"status":'true', "secret":'false', "recommend":'false'}
            for k in ["status", "secret", "recommend"]:
                if k in genome_dict:
                    insert_dict[k] = genome_dict[k]
            self.update_db_record(
                "sg_genome_db",
                record_id=genome_dict["_id"],
                insert_dict=insert_dict
            )

    def update_status2(self):
        "根据genome_id 将一个物种从 tsg数据库 更新到isanger"
        collection = self.db["sg_genome_db"]
        genome_dicts = collection.find({"secret":True})
        for genome_dict in genome_dicts:
            # print genome_dict
            insert_dict={"secret":'true'}
            print genome_dict
            self.update_db_record(
                "sg_genome_db",
                record_id=genome_dict["_id"],
                insert_dict=insert_dict
            )


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("python genome_db.py add_v1 GM0010 test.json" + "\n" +
                        "python genome_db.py add_v2 GM0010 test.json genome_version annot_version\n" +
                        "python genome_db.py from_isanger\n" +
                        "python genome_db.py to_small_rna GM0010\n" +
                        "python genome_db.py to_isanger GM0010\n")
    if sys.argv[1] == "add_v2":
        if len(sys.argv) == 4:
            GenomeDb(None).add_single_genome_db(sys.argv[2], sys.argv[3])
        elif len(sys.argv) == 5:
            GenomeDb(None).add_single_genome_db(sys.argv[2], sys.argv[3], sys.argv[4])
        elif len(sys.argv) > 5:
            GenomeDb(None).add_single_genome_db(sys.argv[2], sys.argv[4], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "add_v1":
        if len(sys.argv) >= 4:
            GenomeDb(None).add_genome_db_from_v1(sys.argv[2], sys.argv[3])
        else:
            pass
    elif sys.argv[1] == "from_isanger":
        GenomeDb(None).updata_from_isanger()
    elif sys.argv[1] == "to_isanger":
        GenomeDb(None).updata_to_isanger(sys.argv[2])
    elif sys.argv[1] == "to_smallrna":
        GenomeDb(None).updata_to_smallrna(sys.argv[2])
    elif sys.argv[1] == "update_status":
        api = GenomeDb(None)
        # print api._config
        if len(sys.argv) >= 3:
            api._config.DBVersion = sys.argv[2]
        api.update_status()
    elif sys.argv[1] == "update_status2":
        api = GenomeDb(None)
        # print api._config
        if len(sys.argv) >= 3:
            api._config.DBVersion = sys.argv[2]
        api.update_status2()

    else:
        raise Exception("不可识别该版本{}".format(sys.argv[1]))
