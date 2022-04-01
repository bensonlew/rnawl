# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

from __future__ import print_function

import json
import logging
import os
import random
import types
from collections import OrderedDict

from biocluster.config import Config
from bson import SON
from bson.objectid import ObjectId

from mainapp.models.workflow import Workflow
# from .core.base import Base


class Wgcna(object):
    def __init__(self, project_type="whole_transcriptome"):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

    def get_bam_path(self, sample_name, task_id):
        collection = self.db["specimen"]
        main_info = collection.find_one({'old_name': sample_name, 'task_id': task_id, "about_qc": "after"})
        bam_path = main_info["bam_path"]
        return bam_path

    def get_new_id(self, task_id, otu_id=None):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
        """
        if otu_id:
            new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
        else:
            # id_ = '%f' % time.time()
            # ids = str(id_).strip().split(".")
            # new_id = "{}_{}_{}".format(task_id, ids[0][5:], ids[1])  #改成时间来命名workflow id
            new_id = "{}_{}_{}".format(task_id, random.randint(1000, 10000), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id

    def get_main_info(self, main_id, collection_name, task_id):
        """
        根据主表id获得整个主表/记录的信息
        :param main_id:
        :param collection_name:
        :return: 整个主表信息
        """
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        main_info = collection.find_one({'main_id': main_id, 'task_id': task_id})
        return main_info

    def get_power(self, prepare_id):
        if isinstance(prepare_id, types.StringTypes):
            main_id = ObjectId(prepare_id)
        elif isinstance(prepare_id, ObjectId):
            main_id = prepare_id
        else:
            raise Exception("prepare_id参数必须为字符串或者ObjectId类型!")
        collection = self.db["wgcna_prepare_detail"]
        prepare_info = collection.find_one({'prepare_id': main_id})
        return prepare_info

    def get_main_info_by_record(self, main_table, **kwargs):
        """
        主表字段查询信息
        :param task_id:
        :return:
        """
        collection = self.db[main_table]
        result = collection.find_one(kwargs)
        return result

    def get_geneset_name(self, geneset_id):
        if ',' not in geneset_id:
            geneset_id = geneset_id
            conn = self.db['geneset']
            geneset_records = conn.find_one({"main_id": ObjectId(geneset_id)})
            geneset_name = geneset_records['name']
        else:
            geneset_id_list = geneset_id.split(',')
            geneset_name_list = list()
            for i in geneset_id_list:
                conn = self.db['geneset']
                geneset_records = conn.find_one({"main_id": ObjectId(i)})
                geneset_name = geneset_records['name']
                geneset_name_list.append(geneset_name)
            geneset_name = "|".join(geneset_name_list)
        return geneset_name

    def insert_main_table(self, collection_name, data):
        """
        新建主表或往主表插入一条记录，并返回该条记录的_id
        :param collection_name:
        :param data:
        :return: 插入记录的id
        """
        main_id = self.db[collection_name].insert_one(SON(data)).inserted_id
        conn = self.db[collection_name]
        task_id = conn.find_one({'_id': main_id})['task_id']
        conn.update({'_id': main_id, "task_id": task_id}, {"$set": {'main_id': main_id}}, upsert=True)
        return main_id

    def find_table_id(self, table_name, **kwargs):
        """
        根据提供的关键字信息查询某表中符合条件的记录，返回该条记录的id
        :param table_name: table name
        :param kwargs:
        :return: 返回符合条件的第一个记录的_id
        """
        conn = self.db[table_name]
        return conn.find(kwargs, {"_id": 1})[0]["_id"]

    def get_control_id(self, main_id):
        collection = self.db['specimen_group_compare']
        try:
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            elif isinstance(main_id, ObjectId):
                main_id = main_id
            else:
                raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
            data = collection.find_one({"_id": main_id})
            if "compare_names" in data.keys():
                compare_names = json.loads(data['compare_names'])
                return compare_names
            else:
                print("{}没有分组的对照信息！".format(str(main_id)))

        except Exception:
            print("{}不存在".format(str(main_id)))

    def insert_geneset_info(self, geneset_id, col_name, col_id):
        collection = self.db['geneset']
        geneset_list = []
        if not isinstance(geneset_id, types.StringTypes):
            raise Exception("输入geneset_id参数必须为字符串类型!")
        geneset_list.extend([ObjectId(x) for x in geneset_id.split(",")])
        if isinstance(col_id, types.StringTypes):
            col_id = ObjectId(col_id)
        elif isinstance(col_id, ObjectId):
            pass
        else:
            raise Exception("输入col_id参数必须为字符串或者ObjectId类型!")
        try:
            for geneset_id in geneset_list:
                result = collection.find_one({"_id": geneset_id})
                result["is_use"] = 1
                collection.update({"_id": geneset_id}, {"$set": result})
        except Exception:
            print("没有找到geneset_id:{}".format(geneset_id))
        try:
            collection = self.db[col_name]
            collection.find_one({"_id": col_id})
        except:
            "没有找到col_id:{} in {}".format(col_id, col_name)
        for geneset_id in geneset_list:
            opts = {"geneset_id": geneset_id, "col_name": col_name, "col_id": col_id}
            collection = self.db["geneset_info"]
            collection.insert_one(opts)
        return True

    def update_group_compare_is_use(self, task_id, main_id):
        collection = self.db['specimen_group_compare']
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            pass
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        result = collection.find_one({"main_id": main_id, "task_id": task_id})
        result["is_use"] = 1
        collection.update({"main_id": main_id, "task_id": task_id}, {"$set": result})
        return True

    def update_group_is_use(self, task_id, main_id, library="long"):
        collection = self.db['specimen_group']
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            pass
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        result = collection.find_one({"main_id": main_id, "task_id": task_id})
        result["is_use"] = 1
        collection.update({"main_id": main_id, "task_id": task_id}, {"$set": result})
        return True

    def get_task_info(self, task_id):
        collection = self.db['task']
        document = collection.find_one({"task_id": task_id})
        return document

    def get_bamlist_path(self, task_id):
        """
        根据task_id到sg_task获得bamlist相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result['bamlist']

    def get_bed_path(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result['bedpath']

    def get_unigenefa_path(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result['unigene_fa']

    def delete_geneset(self, geneset_id, task_id):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            pass
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['geneset_info']
        status_collection = self.db['sg_status']
        results = collection.find({"geneset_id": geneset_id})
        for result in results:
            if "table_name" in result:
                col_name = result["table_name"]
            else:
                col_name = result["col_name"]
            if "table_id" in result:
                col_id = result["table_id"]
            else:
                col_id = result["col_id"]
            print(col_id)
            col = self.db[col_name]
            print(col_name)
            try:
                col_result = col.find_one({"_id": col_id})
                if col_result["task_id"] == task_id:
                    col_result["params"] = ""
                    col.update({"_id": col_id}, {"$set": col_result})
            except:
                print("不能找到对应id {} in {}".format(col_id, col_name))
            try:
                col_result = status_collection.find_one({"table_id": col_id})
                col_result["status"] = "deleted"
                status_collection.update({"table_id": col_id}, {"$set": col_result})
            except:
                print("不能找到对应table_id {} in {}".format(col_id, col_result))
        collection = self.db["geneset"]
        collection_detail = self.db["geneset_detail"]
        result = collection.find_one({"main_id": geneset_id, "task_id": task_id})
        if result:
            collection_detail.delete_many({"geneset_id": result["_id"]})
            collection.remove({"main_id": geneset_id, "task_id": task_id})
        return True

    def delete_main_table(self, collection, task_id):
        table = self.db[collection]
        status = self.db["sg_status"]
        result = table.find_one({"task_id": task_id})
        if result:
            table.remove({"task_id": task_id})
            status.remove({"table_id": result["_id"], "task_id": task_id})

    def insert_seq(self, mongo_data):
        """
        用于前端查询序列信息. 每个task将只存一条序列信息
        :param mongo_data:
        :return:
        """
        task_id = mongo_data["task_id"]
        collection = self.db['query_seq']
        result = collection.find_one({"task_id": task_id})
        if result:
            collection.update({"task_id": task_id}, {"$set": mongo_data})
        else:
            collection.insert_one(mongo_data)

    def insert_rmatsmodel(self, mongo_data):
        """
        用于交互可变剪切模式图分析. 每个task将只存一条分析记录
        :param mongo_data:
        :return:
        """
        task_id = mongo_data["task_id"]
        collection = self.db['splicing_rmats_model']
        result = collection.find_one({"task_id": task_id})
        if result:
            collection.update({"task_id": task_id}, {"$set": mongo_data})
        else:
            collection.insert_one(mongo_data)

    def insert_godag(self, mongo_data):
        """
        用于GO富集邮箱无环图分析. 每个task将只存一条分析记录
        :param mongo_data:
        :return:
        """
        task_id = mongo_data["task_id"]
        collection = self.db['geneset_go_dag']
        result = collection.find_one({"task_id": task_id})
        if result:
            collection.update({"task_id": task_id}, {"$set": mongo_data})
        else:
            result = collection.insert_one(mongo_data)
        return result['_id']

    def insert_none_table(self, collection):
        return self.db[collection].insert_one({}).inserted_id

    def get_snp_info(self, task_id, is_report):
        collection = self.db['snp']
        result = collection.find_one({"task_id": task_id, "is_report": is_report})
        return result

    def get_callvcf_path(self, task_id, is_report):
        """
        这个也用来获取snp主表的sample_names字段的值，获取他的个数
        根据task_id到sg_snp获得call_vcf相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['snp']
        result = collection.find_one({"task_id": task_id, "is_report": is_report})
        return result['call_vcf_path'], len(result['sample_names'])

    def get_ssr_info(self, task_id):
        collection = self.db['ssr']
        result = collection.find_one({"task_id": task_id})
        return result

    def update_status_failed(self, collection, doc_id):
        """
        改特定_id主表的status状态从start为failed，主要用于特殊投递任务失败

        params collection: 主表collection名称
        params doc_id: 主表_id
        """
        self.db[collection].update_one({'_id': ObjectId(doc_id), "status": "start"}, {"$set": {'status': 'failed'}})

    def get_libtype(self, task_id):
        collection = self.db['exp']
        result = collection.find_one({"task_id": task_id})
        if not result or ('libtype' not in result):
            return 'fr'
        return result['libtype']

    def get_mean_read_len(self, task_id):
        collection = self.db['specimen']
        result = collection.find_one({"task_id": task_id, "about_qc": "after"})
        if not result:
            return 149
        mean_len = float(result['total_bases']) / float(result["total_reads"])
        return round(mean_len)

    def get_group_dict(self, main_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['specimen_group']
        result = collection.find_one({"_id": main_id})
        group_names = result["category_names"]
        sample_list = [x.keys() for x in result["specimen_names"]]
        group_dict = OrderedDict(zip(group_names, sample_list))
        return group_dict

    def get_count_path(self, main_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['exp']
        result = collection.find_one({"_id": main_id})
        if result["exp_type"] == "FPKM":
            result1 = collection.find_one(
                {"task_id": result["task_id"], "method": "RSEM", "exp_level": result["exp_level"]})
            return result1['count_file']
        else:
            return result['count_file']

    def get_exp_params_info(self, main_id, task_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['exp']
        result = collection.find_one({'main_id': main_id, 'task_id': task_id})
        return result

    def get_exp_id(self, level, task_id):
        # if isinstance(main_id, types.StringTypes):
        #     main_id = ObjectId(main_id)
        # elif isinstance(main_id, ObjectId):
        #     main_id = main_id
        # else:
        #     raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['exp']
        result = collection.find_one({'level': level, 'task_id': task_id})
        return result

    def get_genesetkeggenrich_params_info(self, main_id, task_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['geneset_kegg_enrich']
        result = collection.find_one({'main_id': main_id, 'task_id': task_id})
        return result

    def get_table_info_by_task_id(self, task_id, table_name="sg_task"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id})
        return result

    def get_json(self):
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    def get_annotation_stat_info(self, task_id, table_name="annotation_stat"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id})
        return result

    def get_des_type(self, task_id, table_name="annotation_stat"):
        collection_task = self.db['task']
        result_task = collection_task.find_one({"task_id": task_id})
        if "genome_id" in result_task.keys():
            project_type = 'ref_rna_v2'
            db = Config().get_mongo_client(mtype=project_type, dydb_forbid=True)[Config().get_mongo_dbname(project_type, dydb_forbid=True)]
            genome_id = result_task["genome_id"]
            col = db["sg_genome_db"]
            genome_info = col.find_one({"genome_id": genome_id})
            db_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish"
            des = os.path.join(db_path, genome_info["bio_mart_annot"])
            des_type = genome_info["biomart_gene_annotype"]
            ref_loc = genome_info["dna_fa"]
            ref_fa = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/" + ref_loc
            return des, des_type, ref_fa
        else:
            collection = self.db[table_name]
            result = collection.find_one({"task_id": task_id, "type": "origin"})
            species_name = result['species_name']
            self.json_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
            self.json_dict = self.get_json()
            des = os.path.join(os.path.split(self.json_path)[0], self.json_dict[species_name]["bio_mart_annot"])
            des_type = self.json_dict[species_name]["biomart_gene_annotype"]
            return des, des_type

    def get_pep(self, task_id, table_name="annotation_stat"):
        collection_task = self.db['sg_task']
        result_task = collection_task.find_one({"task_id": task_id})
        if "genome_id" in result_task.keys():
            genome_id = result_task["genome_id"]
            col = self.db["genome_db"]
            genome_info = col.find_one({"genome_id": genome_id})
            db_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish"
            pep = os.path.join(db_path, genome_info["pep"])
            return pep
        else:
            collection = self.db[table_name]
            result = collection.find_one({"task_id": task_id, "type": "origin"})
            species_name = result['species_name']
            self.json_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
            self.json_dict = self.get_json()
            pep = os.path.join(os.path.split(self.json_path)[0], self.json_dict[species_name]["pep"])
            return pep

    # -------------判断富集表里是否有新基因--------------
    def consistent_in_enrich_and_diff(self, main_id, diff_id, typee="KEGG", table_name="geneset_kegg_enrich"):
        if typee == "KEGG":
            table_name = "geneset_kegg_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        else:
            table_name = "geneset_go_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        if result:
            geneset_id = json.loads(result['params'])['geneset_id']
        else:
            raise Exception("没有找到指定的富集主表：{}".format(main_id))
        geneset_list = self.db['geneset_detail'].find_one({"geneset_id": ObjectId(geneset_id)})['seq_list']
        has_new = False
        if geneset_list:
            for gene in geneset_list:
                if gene.startswith('MSTR') or gene.startswith('TCON') or gene.startswith('XLOC'):
                    has_new = True
        diff_main = self.db['diff'].find_one({"main_id": ObjectId(diff_id)})
        diff_type = json.loads(diff_main['params'])['kind']
        if diff_type == 'ref':
            if has_new:
                return False
            else:
                return True
        else:
            return True

    def update_db_record(self, table_name, record_id, **kwargs):
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception('type of record_id must be ObjectId')
        conn = self.db[table_name]
        conn.update({'_id': record_id}, {'$set': kwargs}, upsert=True)

    def get_genome_info_by_genome_id(self, genome_id):
        """
        根据genomei_id 获取基因组信息
        :param task_id:
        :return:
        """
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        collection = db['sg_genome_db']
        result = collection.find_one({
            "genome_id": genome_id
        })
        return result

    def get_samples_by_task_id(self, task_id):
        collection = self.db['specimen']
        results = collection.find({"task_id": task_id, "library": "small"})
        if type(results) == dict():
            results = [results]
        else:
            results = list(results)
        # print("*** result is ")
        # print(results)
        sample_list = [result["old_name"] for result in results]
        return sample_list

    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i + 3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                self.logger.warn('fail to insert records into table {} -> ({})'.format(table_name, e))
            else:
                self.logger.warn('fail to insert record into table {} -> ({})'.format(table_name, e))
        else:
            if record_num >= 2:
                self.logger.info('succeed in inserting records into table {}'.format(table_name))
            else:
                self.logger.info('succeed in inserting record into table {}'.format(table_name))
            return table_id


if __name__ == "__main__":
    print('no test Now')
