# -*- coding: utf-8 -*-
from __future__ import print_function
from bson.objectid import ObjectId
import types
from bson import SON
from mainapp.models.workflow import Workflow
# from .core.base import Base
import random
import json
from collections import OrderedDict
import os
from biocluster.config import Config
__author__ = 'gdq'


class Wgcna(object):
    def __init__(self, project_type="medical_transcriptome"):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


    def get_bam_path(self, sample_name, task_id):
        collection = self.db["sg_specimen"]
        main_info = collection.find_one({'old_name': sample_name, 'task_id':task_id, "about_qc": "after"})
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

    def get_enrich_info_by_geneset(self, main_id, collection_name, task_id):
        """
        根据主表id获得整个主表/记录的信息
        :param main_id:
        :param collection_name:
        :return: 整个主表信息
        """
        # if isinstance(main_id, types.StringTypes):
        #     main_id = ObjectId(main_id)
        # elif isinstance(main_id, ObjectId):
        #     main_id = main_id
        # else:
        #     raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        main_info = collection.find_one({'geneset_id': main_id, 'task_id': task_id,"status":"end"})
        return main_info

    def get_power(self, prepare_id):
        if isinstance(prepare_id, types.StringTypes):
            main_id = ObjectId(prepare_id)
        elif isinstance(prepare_id, ObjectId):
            main_id = prepare_id
        else:
            raise Exception("prepare_id参数必须为字符串或者ObjectId类型!")
        collection = self.db["sg_wgcna_prepare_detail"]
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
        collection = self.db['sg_specimen_group_compare']
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
        collection = self.db['sg_geneset']
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
            collection = self.db["col_name"]
            collection.find_one({"_id": col_id})
        except:
            "没有找到col_id:{} in {}".format(col_id, col_name)
        for geneset_id in geneset_list:
            opts = {"geneset_id": geneset_id, "col_name": col_name, "col_id": col_id}
            collection = self.db["sg_geneset_info"]
            collection.insert_one(opts)
        return True

    def update_group_compare_is_use(self, task_id, main_id):
        collection = self.db['sg_specimen_group_compare']
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

    def update_group_is_use(self, task_id, main_id):
        collection = self.db['sg_specimen_group']
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

    def get_genesetname(self,task_id=None,geneset_id = None):
        '''
        根据task_id和geneset_id获取基因集名称
        '''
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            pass
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_geneset']
        result = collection.find_one({"task_id": task_id,"main_id":geneset_id})
        geneset_name = result["name"]
        return geneset_name


    def get_task_info(self, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def get_rawgroup_id(self,task_id):
        """
        根据task_id到sg_specimen_group获得原始group表信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_specimen_group']
        result = collection.find_one({"task_id": task_id})
        return result["main_id"]

    def get_align_method(self, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_mapping']
        result = collection.find_one({"task_id": task_id})
        align_method = result["method"]
        return align_method


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

    def delete_difftable(self,diff_id,task_id):
        if isinstance(diff_id, types.StringTypes):
            diff_id = ObjectId(diff_id)
            print("fansile")
        elif isinstance(diff_id, ObjectId):
            pass
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        print(diff_id)
        status_collection = self.db['sg_status']
        try:
            col_result = status_collection.find_one({"table_id": diff_id,"task_id":task_id})
            col_result["status"] = "deleted"
            status_collection.update({"table_id": diff_id}, {"$set": col_result})
        except:
            print("不能找到对应table_id {} in {}".format(diff_id, col_result))
        col = self.db["sg_diff"]
        try:
            col_result = col.find_one({"main_id": diff_id,"task_id":task_id })
            print("qishizhaodaole")
            if col_result["task_id"] == task_id:
                col_result["params"] = ""
                col.update({"main_id": diff_id,"task_id":task_id}, {"$set": col_result})
        except:
            print("不能找到对应id {} in {}".format(diff_id, "sg_diff"))
        if col_result:
            geneset_ids = col_result["genesets"]
            for geneset_id in  geneset_ids:
                self.delete_geneset(geneset_ids[geneset_id],task_id)
        return True



    def delete_geneset(self, geneset_id, task_id):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            pass
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_geneset_info']
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
        collection = self.db["sg_geneset"]
        collection_detail = self.db["sg_geneset_detail"]
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

    def delete_genome_db_table(self, genome_id):
        table = self.db["sg_genome_db"]
        result = table.find_one({"genome_id": genome_id})
        if result:
            table.remove({"genome_id": genome_id})

    def insert_seq(self, mongo_data):
        """
        用于前端查询序列信息. 每个task将只存一条序列信息
        :param mongo_data:
        :return:
        """
        task_id = mongo_data["task_id"]
        collection = self.db['sg_query_seq']
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
        collection = self.db['sg_splicing_rmats_model']
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
        collection = self.db['sg_geneset_go_dag']
        result = collection.find_one({"task_id": task_id})
        if result:
            collection.update({"task_id": task_id}, {"$set": mongo_data})
        else:
            result = collection.insert_one(mongo_data)
        return result['_id']

    def insert_none_table(self, collection):
        return self.db[collection].insert_one({}).inserted_id

    def get_snp_info(self, task_id, is_report):
        collection = self.db['sg_snp']
        result = collection.find_one({"task_id": task_id, "is_report": is_report})
        return result

    def get_callvcf_path(self, task_id, is_report):
        """
        这个也用来获取snp主表的sample_names字段的值，获取他的个数
        根据task_id到sg_snp获得call_vcf相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_snp']
        result = collection.find_one({"task_id": task_id, "is_report": is_report})
        return result['call_vcf_path'], len(result['sample_names'])

    def get_ssr_info(self, task_id):
        collection = self.db['sg_ssr']
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
        collection = self.db['sg_exp']
        result = collection.find_one({"task_id": task_id})
        if not result or ('libtype' not in result):
            return 'fr'
        return result['libtype']

    def get_mean_read_len(self, task_id):
        collection = self.db['sg_specimen']
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
        collection = self.db['sg_specimen_group']
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

        collection = self.db['sg_exp']
        result = collection.find_one({"_id": main_id})
        if 'count_file' in result and result['count_file'].startswith("s3"):
            # 判断有无定量文件
            count_file = result['count_file']
        else:
            result1 = collection.find_one({"task_id": result["task_id"], "method" : "RSEM", "level": result["level"]})
            count_file = result1['count_file']
        return count_file

        # try:
        #     return result['count_file']
        # except:
        #     raise Exception("Can't find count_file!")

    def get_exp_params_info(self, main_id, task_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_exp']
        result = collection.find_one({'main_id': main_id, 'task_id': task_id})
        return result

    def get_genesetkeggenrich_params_info(self,main_id,task_id):
        if isinstance(main_id,types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id,ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection=self.db['sg_geneset_kegg_enrich']
        result = collection.find_one({'main_id':main_id,'task_id':task_id})
        return result

    def get_table_info_by_task_id(self, task_id, table_name="sg_task"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id})
        return result

    def get_json(self):
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    def get_annotation_stat_info(self, task_id, table_name="sg_annotation_stat"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        return result

    # def get_des_type(self, task_id, table_name="sg_annotation_stat"):
    #     collection_task = self.db['sg_task']
    #     result_task = collection_task.find_one({"task_id": task_id})
    #     if "genome_id" in result_task.keys():
    #         genome_id = result_task["genome_id"]
    #         col = self.db["sg_genome_db"]
    #         genome_info = col.find_one({"genome_id": genome_id})
    #         db_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish"
    #         des = os.path.join(db_path, genome_info["bio_mart_annot"])
    #         des_type = genome_info["biomart_gene_annotype"]
    #         return des, des_type
    #     else:
    #         collection = self.db[table_name]
    #         result = collection.find_one({"task_id": task_id, "type": "origin"})
    #         species_name = result['species_name']
    #         self.json_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
    #         self.json_dict = self.get_json()
    #         des = os.path.join(os.path.split(self.json_path)[0], self.json_dict[species_name]["bio_mart_annot"])
    #         des_type = self.json_dict[species_name]["biomart_gene_annotype"]
    #         return des, des_type

    def get_des_type(self, task_id, table_name="annotation_stat"):
        collection_task = self.db['sg_task']
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

    def get_pep(self, task_id, table_name="sg_annotation_stat"):
        collection_task = self.db['sg_task']
        result_task = collection_task.find_one({"task_id": task_id})
        if "genome_id" in result_task.keys():
            genome_id = result_task["genome_id"]
            col = self.db["sg_genome_db"]
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

    #-------------判断富集表里是否有新基因--------------
    def consistent_in_enrich_and_diff(self, main_id, diff_id, typee="KEGG", table_name="sg_geneset_kegg_enrich"):
        if typee == "KEGG":
            table_name = "sg_geneset_kegg_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "DISGENET":
            table_name = "sg_geneset_disgenet_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "DO":
            table_name = "sg_geneset_do_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "REACTOME":
            table_name = "sg_geneset_reactome_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        else:
            table_name = "sg_geneset_go_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        if result:
            geneset_id=json.loads(result['params'])['geneset_id']
        else:
            raise Exception("{}没有找到指定的富集主表：{}".format(table_name,main_id))
        geneset_list = self.db['sg_geneset_detail'].find_one({"geneset_id": ObjectId(geneset_id)})['seq_list']
        has_new=False
        if geneset_list:
            for gene in geneset_list:
                if gene.startswith('MSTR') or gene.startswith('TCON') or gene.startswith('XLOC'):
                    has_new = True
        diff_main = self.db['sg_diff'].find_one({"main_id": ObjectId(diff_id)})
        diff_type = json.loads(diff_main['params'])['kind']
        if diff_type == 'ref':
            if has_new:
                return False
            else:
                return True
        else:
            return True

    # def get_new_id(self, task_id, otu_id=None):
    #     """
    #     根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
    #     """
    #     if otu_id:
    #         new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
    #     else:
    #         # id_ = '%f' % time.time()
    #         # ids = str(id_).strip().split(".")
    #         # new_id = "{}_{}_{}".format(task_id, ids[0][5:], ids[1])  #改成时间来命名workflow id
    #         new_id = "{}_{}_{}".format(task_id, random.randint(1000, 10000), random.randint(1, 10000))
    #     workflow_module = Workflow()
    #     workflow_data = workflow_module.get_by_workflow_id(new_id)
    #     if len(workflow_data) > 0:
    #         return self.get_new_id(task_id, otu_id)
    #     return new_id

    def consistent_in_enrich_and_diff_diff(self, main_id, diff_id, typee="KEGG", table_name="sg_geneset_kegg_enrich"):
        if typee == "KEGG":
            table_name = "sg_diff_geneset_kegg_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "DISGENET":
            table_name = "sg_diff_geneset_disgenet_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "DO":
            table_name = "sg_diff_geneset_do_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        elif typee == "REACTOME":
            table_name = "sg_diff_geneset_reactome_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        else:
            table_name = "sg_diff_geneset_go_enrich"
            collection = self.db[table_name]
            result = collection.find_one({"main_id": ObjectId(main_id)})
        if result:
            geneset_id=json.loads(result['params'])['geneset_id']
        else:
            raise Exception("{}没有找到指定的富集主表：{}".format(table_name,main_id))
        geneset_list = self.db['sg_geneset_detail'].find_one({"geneset_id": ObjectId(geneset_id)})['seq_list']
        has_new=False
        if geneset_list:
            for gene in geneset_list:
                if gene.startswith('MSTR') or gene.startswith('TCON') or gene.startswith('XLOC'):
                    has_new = True
        diff_main = self.db['sg_diff'].find_one({"main_id": ObjectId(diff_id)})
        diff_type = json.loads(diff_main['params'])['kind']
        if diff_type == 'ref':
            if has_new:
                return False
            else:
                return True
        else:
            return True

    def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
        if record_id is not None:
            if isinstance(record_id, types.StringTypes):
                record_id = ObjectId(record_id)
            elif isinstance(record_id, ObjectId):
               record_id = record_id
            else:
                self.bind_object.set_error('type of main id must be String or ObjectId', code="53701107")
        conn = self.db[table_name]
        if query_dict:
            if record_id is not None:
                query_dict.update({'_id': record_id})
        else:
            if record_id is not None:
                query_dict = {'_id': record_id}
            else:
                self.bind_object.set_error('query dict must be provided while record id is None', code="53701108")
        if insert_dict:
            kwargs.update(insert_dict)
        conn.update(query_dict, {'$set': kwargs}, upsert=True)

    def delete_search_records(self,table_name,task_id,remain_number):
        if table_name == "sg_asprofile_search":
            collection = self.db[table_name]
            num = collection.find({"task_id": task_id,"status":"end"}).count()
            if num == remain_number:
                deleted_record = collection.find_one({"task_id": task_id,"status":"end"})
                deleted_main_id = deleted_record["main_id"]
                collection.remove({"main_id": deleted_main_id, "task_id": task_id})
                detail_collection = self.db["sg_asprofile_search_detail"]
                detail_collection.delete_many({"asprofile_search_id": deleted_main_id})
        if table_name == "sg_snp_search":
            collection = self.db[table_name]
            num = collection.find({"task_id": task_id,"status":"end"}).count()
            if num == remain_number:
                deleted_record = collection.find_one({"task_id": task_id,"status":"end"})
                deleted_main_id = deleted_record["main_id"]
                collection.remove({"main_id": deleted_main_id, "task_id": task_id})
                detail_collection = self.db["sg_snp_search_detail"]
                detail_collection.delete_many({"snp_search_id": deleted_main_id})


if __name__ == "__main__":
    print('no test Now')
