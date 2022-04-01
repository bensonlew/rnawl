# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
from bson import SON
# from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
import random
import json


class Metagenomic(Meta):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Metagenomic, self).__init__(self._bind_object)
        self._project_type = "metagenomic"
        # self.db_name = Config().MONGODB + '_metagenomic'

    def get_geneset_info(self, geneset_id):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            geneset_id = geneset_id
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['geneset']
        print collection
        result = collection.find_one({"_id": geneset_id})
        return result

    def get_sample_name(self, geneset_id):
        result = self.get_geneset_info(geneset_id)
        sample_str = result['specimen']
        sample_list = sample_str.split(',')
        return sample_list

    def get_anno_path(self, geneset_id, database, task_id, nr_method="best_hit"):  ## add by shaohua.yuan
        # from geneset_id get anno path
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            geneset_id = geneset_id
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        anno_main_table = "anno_" + database
        collection = self.db[anno_main_table]
        if database == "cazy":
            origin = "CAZy_Origin"
        else:
            origin = database.upper() + "_Origin"
        if database in ["cog","kegg","cazy","ardb","card","vfdb"]:
            result = collection.find_one({"geneset_id": geneset_id,"task_id": task_id,"name":origin})
        elif database in ["nr"]:
            if nr_method == "best_hit":
                origin = origin
            elif  nr_method == "lca":
                origin = database.upper() + "_Origin_LCA"
            elif  nr_method == "de_unclassified":
                origin = database.upper() + "_Origin_Deunclassified"
            result = collection.find_one({"geneset_id": geneset_id,"task_id": task_id,"name":origin})
        elif database == "probio":
            result1 = collection.find_one({"geneset_id": geneset_id,"task_id": task_id, "is_origin":1})
            result = collection.find_one({"geneset_id": geneset_id,"task_id": task_id, "is_origin":1, "nr_method": nr_method})
            if result1 and not result:
                result = "gene_anno_info"
        else:
            result = collection.find_one({"geneset_id": geneset_id,"task_id": task_id, "is_origin" : 1})
        if not result:
            raise Exception('anno_{}没有相应geneset_id:{}对应的信息 client{}!'.format(database,geneset_id, self.db))
        return result

    def find_origin_genesetID(self, collection, task_id):  # add by shaohua.yuan
        # 根据task_id查找原始基因集
        collection = self.db[collection]
        result = collection.find_one({"$and": [{"task_id": task_id}, {"type": 1}]})
        if result:
            geneset_id = result["_id"]
        else:
            raise Exception('没有相应geneset表！')
        return geneset_id

    def from_id_get_result(self,collection,main_id, main_name="_id", condition=None): # add by shaohua.yuan
        if collection == "anno_mvirdb":
            collection = "anno_mvir"
        collection = self.db[collection]
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        if condition:
            condition.update({main_name: main_id})
            result = collection.find_one(condition)
        else:
            result = collection.find_one({"_id": main_id})
        if not result:
            raise Exception('{}没有相应的id{}！'.format(collection,main_id))
        return result


    def export_geneset_table(self, geneset_id, method):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            geneset_id = geneset_id
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['geneset']
        result = collection.find_one({"_id": geneset_id})
        if not result:
            raise Exception('基因集id没有找到对应的表信息')
        task_id = result['task_id']
        gene_list = result['gene_list_length']
        if not gene_list:
            raise Exception('找不到gene_list文件')
        basic_main = collection.find_one({'task_id': task_id, 'type': 1})
        if basic_main.has_key(method):
            geneset_abu = basic_main[method]
        elif method == "ppm":
            geneset_abu = basic_main["reads_num"] # 通过该文件计算
        if not geneset_abu:
            if method == "ppm":
                geneset_abu = basic_main["reads_num"] # 通过该文件计算
            else:
                raise Exception('找不到基因集%s丰度表' % method)
        geneset_abu = self.use_s3(geneset_abu)
        gene_list = self.use_s3(gene_list)
        return geneset_abu,gene_list

    def use_s3(self, path):
        if path.startswith("rere"):
            return "s3://" + path  # 防止tsg_31796任务路径出错
        else:
            return path

    def get_anno_info(self,anno_id,anno_type):
        if isinstance(anno_id, types.StringTypes):
            anno_id = ObjectId(anno_id)
        elif isinstance(anno_id, ObjectId):
            anno_id = anno_id
        else:
            raise Exception("输入anno_id参数必须为字符串或者ObjectId类型!")
        if anno_type == "mvirdb":
            anno_type = "mvir"
        anno_db_name = "anno_" + anno_type
        collection = self.db[anno_db_name]
        result = collection.find_one({"_id": anno_id})
        return result

    def get_nr_info(self, nr_method, task_id):
        method_map_name = {
            "best_hit": "NR_Origin",
            "de_unclassified": "NR_Origin_Deunclassified",
            "lca": "NR_Origin_LCA"
        }
        collection = self.db["anno_nr"]
        result = collection.find_one({"name": method_map_name[nr_method], "task_id": task_id})
        return result

    def get_abu_file(self,params,task_id):
        collection = self.db['abund_table_path']
        basic_main = collection.find_one({'task_id': task_id, 'params': params})
        if not basic_main:
            return 0
        else:
            return basic_main
            # if basic_main['status'] == 'start':
            #     return 1

    def get_group_name(self, group_id, second_group=''):
        """
        根据分组方案id获取分组方案名字
        :param group_id: 分组方案id
        :return: 分组方案名字
        """
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, types.StringTypes):
                group_id = ObjectId(group_id)
        else:
            raise Exception("group_detail必须为ObjectId对象或其对应的字符串!")
        collection = self.db['specimen_group']
        result = collection.find_one({'_id': group_id})
        gname = result['group_name']
        if second_group != '"null"':
            gname = gname + ',' + 'second_group'
        return gname

    def find_origin_geneset_info(self, task_id):  # add by zouguanqing
        # 根据task_id查找原始基因集
        collection = self.db['geneset']
        result = collection.find_one({"$and": [{"task_id": task_id}, {"type": 1}]})
        if result:
            return result
        else:
            raise Exception('没有相应geneset表!')

    def find_anno_main_table(self,collection,task_id):   # add by zouguanqing
        collection1 = self.db[collection]
        result = collection1.find_one({"$and": [{"task_id": task_id}, {"is_origin": 1}]})
        if result:
            #if result['status'] == 'hide':
            return result
            #else:
            #    return False
        else:
            return False

    def find_anno_nr_table(self,task_id,name):   # add by zouguanqing
        collection1 = self.db['anno_nr']
        result = collection1.find_one({"$and": [{"task_id": task_id}, {"name": name}]})  # 'NR_Origin'
        if result:
            #if result['status'] == 'end':
            return result
            #else:
            #    return False
        else:
            return False

    def update_ppm_path(self, task_id, path): #add by shaohua.yuan
        collection = self.db['geneset']
        result = collection.find_one({"task_id": task_id})
        if not result:
            raise Exception('task_id {} not in geneset mongo!')
        geneset_id = result["_id"]
        collection.update_one({'_id': geneset_id}, {'$set': {'ppm': path}})

    def common_find_one(self,db_name, condition):   #zouguanqing
        collection = self.db[db_name]
        result = collection.find_one(condition)
        if not result:
            return None
        else:
            return result

    def common_update_one(self, db_name, table_id, update_value):   #zouguanqing
        collection = self.db[db_name]
        collection.update_one({'_id':table_id},{"$set":update_value})

    def rm_main(self, task_id, db_name, main_id):
        """
        根据task_id、主表名称和主表id删除对应的主表
        :param task_id:
        :param db_name:
        :param main_id:
        :return:
        author:qingchen.zhang @20190927
        """
        if isinstance(main_id, types.StringTypes):
            main_table_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_table_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        demo_list = []
        sg_task = self.db['sg_task'].find({})
        new_task_id = task_id.rstrip('_1')
        for demo in sg_task:
            demo_id = demo['demo_id']
            is_demo = demo['is_demo']
            new_task = demo['task_id']
            if is_demo == 1:
                demo_list.append(new_task)
                if demo_id not in demo_list:
                    demo_list.append(demo_id)
        print len(demo_list)
        if new_task_id not in demo_list:
            main_db = self.db[db_name]
            results = main_db.find({'_id': main_table_id})
            if results:
                main_db.delete_many({'_id': main_table_id})
            else:
                raise Exception("根据main_id未找到对应的结果!")
        else:
            print("task任务为demo不能删除数据")

    def rm_main_detail(self, task_id, db_name, main_id, table_id):
        """
        从table_detail中删除对应主表的数据，在更新主表的时候删除原来已经存在的数据
        逻辑是：先去判断sg_task表中的demo_id，如果没有这个demo_id，那么去执行删除任务，否则的话就执行跳过
        :param task_id:任务id
        :param db_name:对应详情表名称
        :param main_id:详情表对应主表的main_id
        :param table_id:详情表对应主表main_id的名称
        :return:
        author:qingchen.zhang@20190927
        """
        if isinstance(main_id, types.StringTypes):
            main_table_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_table_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        demo_list = []
        sg_task = self.db['sg_task'].find({})
        new_task_id = task_id.rstrip('_1')
        for demo in sg_task:
            demo_id = demo['demo_id']
            is_demo = demo['is_demo']
            new_task = demo['task_id']
            if is_demo == 1:
                demo_list.append(new_task)
                if demo_id not in demo_list:
                    demo_list.append(demo_id)
        print (demo_list)
        print new_task_id
        if new_task_id not in demo_list:
            db_detail = self.db[db_name]
            results = db_detail.find({table_id: main_table_id})
            if results:
                db_detail.delete_many({table_id: main_table_id})
            else:
                raise Exception("根据main_id未找到对应的结果!")
            print("删除数据完成")
        else:
            print("task任务为demo不能删除数据")

    def get_nr_version(self, task_id):
        """
        根据输入task_id返回nr版本号，用于判断不同的nr版本
        :param task_id:
        :return:
        """
        nr_db = self.db['anno_nr']
        nr_version = ''
        result = nr_db.find_one({"task_id": task_id, "name": "NR_Origin"})
        if result:
            if "settled_params" in result:
                nr_version = json.loads(result['settled_params'])['version']
            else:
                nr_version = ''
        return nr_version

    def find_version_from_task(self, task_id):
        """
        根据task_id任务查询工作流版本号
        qingchen.zhang 依此来兼容新老项目
        :param task_id:
        :return:
        """
        db = self.db["sg_task"]
        rets = db.find_one({"task_id": task_id})
        if rets:
            if "database" in rets:
                database = json.loads(rets['database'])
                if "update_time" in database:
                    update_time = database['update_time']
                    if update_time in ['202011']:
                        version = "kegg_v94.2"
                        return version
            else:
                version = ""
                return version
        else:
            return 'not find'

    def get_new_env_units(self, units, env_id):
        env_info = self.get_mongo_common('env', ObjectId(env_id))
        units = json.loads(units)
        if 'env_true_name' in env_info:
            true_env = env_info['env_true_name']
        else:
            true_env = None
        if true_env:
            env_names = env_info['env_names'].split(',')
            env_name_dict = dict(zip(true_env, env_names))
            for i in range(len(units)):
                n = units[i]
                units[i] = env_name_dict[n]
        return ','.join(units)
