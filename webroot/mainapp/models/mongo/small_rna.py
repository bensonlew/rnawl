# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

from __future__ import print_function
from bson.objectid import ObjectId
import types
from bson import SON
from mainapp.models.workflow import Workflow
from .core.base import Base
import random
import json
from collections import OrderedDict
from biocluster.config import Config
import datetime

'''
to_file是把数据库里面的内容导出到一个文件，接口里面的那个to_file赋值后面的括号的参数就是文件的名字，这个文件作为tool的输入
small_rna里面定义了很多函数，可以直接获取mongodb里面的数据，获得具体的参数值，如文件路径信息，task_id信息,根据这里面的函数的需求
'''


class SmallRna(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(SmallRna, self).__init__(bind_object=self._bind_object)
        self._project_type = 'small_rna'

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
            new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id

    def get_genesets(self, task_id, geneset_type, is_use=None):
        """
        获取task_id下全部的基因集名称
        """
        geneset_names = list()
        collection = self.db["sg_geneset"]
        if is_use:
            results = collection.find({'task_id': task_id, "status": "end", "type": geneset_type, "is_use": 1})
        else:
            results = collection.find({'task_id': task_id, "status": "end", "type": geneset_type})
        if results:
            for result in results:
                geneset_names.append(result["name"])
        return geneset_names

    def get_geneset_batch_flag(self, task_id):
        collection = self.db["sg_geneset_batch"]
        geneset_batches = collection.find({"task_id": task_id})
        flag_list = list()
        if geneset_batches:
            for geneset_batch in geneset_batches:
                if 'flag' in geneset_batch:
                    print(geneset_batch)
                    flag_list.append(geneset_batch)
            if flag_list:
                flags = sorted(flag_list, key=lambda geneset_batch:int(geneset_batch["flag"]))
                print(flags[-1]["flag"])
                flag = int(flags[-1]["flag"]) + 1
            else:
                flag = 1
        else:
            flag = 1
        return flag

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
        Insert a new document in specific collection, build new one if that collection does not exist.
        :param collection_name: name of the collection (str)
        :param data: main_info of the inserted document with essential key "task_id" (dict)
        :return: _id of the inserted document (ObjectId)
        """
        main_id = self.db[collection_name].insert_one(SON(data)).inserted_id
        conn = self.db[collection_name]
        task_id = conn.find_one({"_id": main_id})["task_id"]
        conn.update({"_id": main_id, "task_id": task_id}, {"$set": {"main_id": main_id}}, upsert=True)
        return main_id

    def delete_main_table(self, collection_name, main_id):
        """
        根据id删除主表
        """
        conn = self.db[collection_name]
        conn.delete_many({"_id": main_id})


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
            data = collection.find_one({"_id":main_id})
            if "compare_names" in data.keys():
                compare_names = json.loads(data['compare_names'])
                return compare_names
            else:
                print("{}没有分组的对照信息！".format(str(main_id)))

        except Exception:
            print("{}不存在".format(str(main_id)))

    def insert_geneset_info(self, geneset_ids, col_name, col_id):
        collection = self.db["sg_geneset"]
        geneset_list = list()
        if not isinstance(geneset_ids, types.StringTypes):
            raise Exception("geneset_ids must be str!")
        geneset_list.extend([ObjectId(x) for x in geneset_ids.split(",")])
        if isinstance(col_id, types.StringTypes):
            col_id = ObjectId(col_id)
        elif isinstance(col_id, ObjectId):
            pass
        else:
            raise Exception("col_id must be str or ObjectId!")
        try:
            for geneset_id in geneset_list:
                result = collection.find_one({"_id": geneset_id})
                result["is_use"] = 1
                collection.update({"_id": geneset_id}, {"$set": result})
        except:
            raise Exception("Can not find _id: {} in sg_geneset".format(geneset_id))
        try:
            collection = self.db[col_name]
            collection.find_one({"_id": col_id})
        except:
            raise Exception("Can not find _id: {} in {}".format(col_id, col_name))
        for geneset_id in geneset_list:
            opts = {
                "geneset_id": geneset_id,
                "col_name": col_name,
                "col_id": col_id
            }
            collection = self.db["sg_geneset_info"]
            collection.insert_one(opts)
        return True

    def get_task_info(self, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def get_genome_info(self, organism_name, assembly, annot_version):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        collection = db['sg_genome_db']
        result = collection.find_one({
            "organism_name": organism_name,
            "assembly": assembly,
            "annot_version": annot_version
        })
        return result

    def get_genome_info_new(self, genome_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        collection = db['sg_genome_db']
        result = collection.find_one({
            "genome_id": genome_id,
        })
        return result

    def get_genome_info_by_id(self, genome_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        collection = db['sg_genome_db']
        result = collection.find_one({
            "genome_id": genome_id,
        })
        return result

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
            raise Exception("geneset_id must be str or ObjectId!")
        collection = self.db["sg_geneset_info"]
        status_collection = self.db['sg_status']
        results = collection.find({"geneset_id":geneset_id})
        for result in results:
            col_name = result["col_name"]
            col_id = result["col_id"]
            print(col_id)
            col = self.db[col_name]
            print(col_name)
            try:
                col_result = col.find_one({"_id":col_id})
                if col_result["task_id"] == task_id:
                    col_result["params"] = ""
                    col.update({"_id":col_id}, {"$set":col_result})
            except:
                print("Can not find _id: {} in {}".format(col_id, col_name))
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

    def check_assest_for_demo(self):
        collection = self.db['sg_task']
        nums = collection.count({"task_id": {"$regex": "refrna_demo_mouse.*"}})
        if nums:
            if nums <= 5:
                return False
            else:
                return True

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

    def update_main_start(self, collection, main_id, mongo_data):
        """
        改特定_id主表的status状态从start为failed，主要用于特殊投递任务失败

        params collection: 主表collection名称
        params doc_id: 主表_id
        """
        self.db[collection].update_one({'_id': ObjectId(main_id)},
                                       {"$set": mongo_data})



    def update_status_failed(self, collection, doc_id):
        """
        改特定_id主表的status状态从start为failed，主要用于特殊投递任务失败

        params collection: 主表collection名称
        params doc_id: 主表_id
        """
        self.db[collection].update_one({'_id': ObjectId(doc_id), "status": "start"},
                                       {"$set": {'status': 'failed'}})

    def get_libtype(self, task_id):
        collection = self.db['sg_exp']
        result = collection.find_one({"task_id": task_id})
        if 'libtype' not in result:
            return None
        return result['libtype']

    def get_mean_read_len(self, task_id):
        collection = self.db['sg_specimen']
        result = collection.find_one({"task_id": task_id, "about_qc": "after"})
        mean_len = float(result['total_bases'])/ float(result["total_reads"])
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

    def get_samples_by_task_id(self, task_id):
        collection = self.db['sg_specimen']
        results = collection.find({"task_id": task_id, "about_qc": "before"})
        if type(results) == dict():
            results = [results]
        else:
            results = list(results)
        # print("*** result is ")
        # print(results)
        sample_list = [result["old_name"] for result in results]
        return sample_list

    def get_count_path(self, main_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_exp']
        result = collection.find_one({"_id": main_id})
        return result['count_file']

    def get_exp_params_info(self, main_id, task_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id must be str or ObjectId!")
        conn = self.db['sg_exp']
        result = conn.find_one({'main_id': main_id, 'task_id': task_id})
        return result

    def get_exp_info(self, task_id):
        "get small rna exp matrix"

        conn = self.db['sg_exp']
        result = conn.find_one({"exp_type" : "tpm", 'task_id': task_id})
        return result

    def get_target_params_info(self, target_id):
        if isinstance(target_id, types.StringTypes):
            target_id = ObjectId(target_id)
        elif isinstance(target_id, ObjectId):
            target_id = target_id
        else:
            raise Exception("target_id must be str or ObjectId!")
        conn = self.db['sg_target']
        result = conn.find_one({'_id': target_id})
        return result

    def get_exp_tpm_params_info(self, task_id):
        conn = self.db['sg_exp']
        result = conn.find_one({'exp_type': "tpm", 'task_id': task_id})
        return result

    def get_count_exp_id(self, main_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception('main_id must be str or ObjectId!')
        conn = self.db['sg_exp']
        tpm_doc = conn.find_one({'main_id': main_id})
        count_doc = conn.find_one({'project_sn': tpm_doc['project_sn'],
                                   'task_id': tpm_doc['task_id'],
                                   'exp_type': 'count'})
        return str(count_doc['main_id'])

    def get_annotation_stat_info(self, task_id, table_name="sg_annotation_stat"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        return result

    def insert_seq(self, mongo_data):
        """
        用于前端查询序列信息，每个task将只存一条序列信息。
        """
        task_id = mongo_data['task_id']
        collection = self.db['sg_query_seq']
        result = collection.find_one({'task_id': task_id})
        if result:
            collection.update({'task_id': task_id}, {'$set': mongo_data})
        else:
            collection.insert_one(mongo_data)

    def update_db_record(self, table_name, record_id, **kwargs):
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception('type of record_id must be ObjectId')
        conn = self.db[table_name]
        conn.update({'_id': record_id}, {'$set': kwargs}, upsert=True)

    def update_db_record_new(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
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

if __name__ == "__main__":
    print('no test Now')
