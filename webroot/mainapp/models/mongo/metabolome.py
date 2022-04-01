# -*- coding: utf-8 -*-
# from mainapp.config.db import get_mongo_client

from bson.objectid import ObjectId
import types
from bson import SON
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
import random, os, re
import json
from boto.s3.bucket import Bucket
#from biocluster.api.file.lib.transfer import TransferManager
from biocluster.file import exists, list_dir


class Metabolome(Meta):
    '''
    获取mongo字段
    '''
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Metabolome, self).__init__(self._bind_object)
        self._project_type = "metabolome"

    def get_metab_table(self, metab_table_id, task_id, type=None):
        metab_table_id = self.check_id(metab_table_id)
        # collection = self.db['metab_table']
        collection = self.db['exp']
        result = collection.find_one({"_id": metab_table_id})
        print result
        if not result:
            raise Exception('{}没有相应的id{}！'.format(collection, metab_table_id))
        if type == "pos":
            table = os.path.join(result['table_path'].split(',')[0], "metab_abund.txt")
        elif type == "neg":
            table = os.path.join(result['table_path'].split(',')[1], "metab_abund.txt")
        elif type == "mix":
            table = os.path.join(result['table_path'].split(',')[2], "metab_abund.txt")
        elif not type:
            table = os.path.join(result['table_path'], "metab_abund.txt")
        else:
            raise Exception('metab_table类型type参数错误')
        if not "/mnt/ilustre/" in table:
            sanger_table = "/mnt/ilustre/data/" + table
            tsanger_table = "/mnt/ilustre/tsanger-data/" + table
            if exists(table):
                table = table
            elif exists(tsanger_table):
                table = tsanger_table
            elif exists(sanger_table):
                table = sanger_table
            else:
                raise Exception('metab_table path not exists-{}'.format(table))
        else:
            if not exists(table):
                raise Exception('metab_table path not exists-{}'.format(table))
        return table

    def get_metab_desc(self, metab_table_id, task_id, type=None):
        metab_table_id = self.check_id(metab_table_id)
        # collection = self.db['metab_table']
        collection = self.db['exp']
        result = collection.find_one({"_id": metab_table_id})
        if not result:
            raise Exception('{}没有相应的id{}！'.format(collection, metab_table_id))
        if type == "pos":
            table = os.path.join(result['table_path'].split(',')[0], "metab_desc.txt")
        elif type == "neg":
            table = os.path.join(result['table_path'].split(',')[1], "metab_desc.txt")
        elif type == "mix":
            table = os.path.join(result['table_path'].split(',')[2], "metab_desc.txt")
        elif not type:
            table = os.path.join(result['table_path'], "metab_desc.txt")
        else:
            raise Exception('metab_table类型type参数错误')
        if not "/mnt/ilustre/" in table:
            sanger_table = "/mnt/ilustre/data/" + table
            tsanger_table = "/mnt/ilustre/tsanger-data/" + table
            if exists(table):
                table = table
            elif exists(sanger_table):
                table = sanger_table
            elif exists(tsanger_table):
                table = tsanger_table
            else:
                raise Exception('metab_table path not exists-{}'.format(table))
        else:
            if not exists(table):
                raise Exception('metab_table path not exists-{}'.format(table))
        return table

    def get_metab_info(self, metab_table_id, task_id):
        metab_table_id = self.check_id(metab_table_id)
        # collection = self.db['metab_table']
        collection = self.db['exp']
        result = collection.find_one({"_id": metab_table_id, "task_id": task_id})
        #result = collection.find_one({"_id": metab_table_id})
        if not result:
            result = collection.find_one({"main_id": metab_table_id, "task_id": task_id})
        if not result:
            raise Exception('{}没有相应的id:{},{}！'.format(collection, metab_table_id, task_id))
        table_name = result['name']
        project_sn = result['project_sn']
        #task_id = result['task_id']
        project_type = result['type']
        return project_sn, project_type, table_name

    def get_project_info(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        if not result:
            raise Exception('{}没有响应的id:{}!'.format(collection, task_id))
        project_sn = result['project_sn']
        member_id = result['member_id']
        project_type = result['type']
        return project_sn, project_type, member_id

    def get_metab_path(self, task_id):
        collection = self.db['exp']
        result = collection.find_one({"task_id": task_id, "is_raw": 1})
        if not result:
            raise Exception('没有查询到相应任务的原始表：{}'.format(task_id))
        path = result['table_path']
        new_list = []
        if not "/mnt/ilustre/" in path:
            tmp_paths = path.split(",")
            print tmp_paths[0] + "metab_abund.txt"
            tmp_path = tmp_paths[0].rstrip("/") + "/metab_abund.txt"
            if exists(tmp_path):
                path_prefix = ""
            elif exists("/mnt/ilustre/data/" + tmp_path):
                path_prefix = "/mnt/ilustre/data/"
            elif exists("/mnt/ilustre/tsanger-data/" + tmp_path):
                path_prefix = "/mnt/ilustre/tsanger-data/"
            else:
                raise Exception('metab_table路径参数错误-{}'.format(path))
            for each in tmp_paths:
                each = each.rstrip("/") + "/"
                tmp_each = path_prefix + each
                new_list.append(tmp_each)
            path = ",".join(new_list)
        return path

    def get_diff_dir(self, diff_table_id, task_id):
        diff_table_id = self.check_id(diff_table_id)
        # collection = self.db['metab_table']
        collection = self.db['exp_diff']
        #result = collection.find_one({"_id": diff_table_id, "task_id": task_id})
        result = collection.find_one({"_id": diff_table_id})
        #result = collection.find_one({"_id": metab_table_id})
        if not result:
            raise Exception('{}没有相应的id:{},{}！'.format(collection, diff_table_id, task_id))
        diff_dir = result['diff_dir'].rstrip("/") + "/"
        my_bucket = Config().get_project_region_bucket(project_type="metabolome")
        #if my_bucket in diff_dir or "/mnt/ilustre/" in diff_dir:
        if my_bucket in diff_dir or "/mnt/ilustre/" in diff_dir or "s3nb://metabolome/" in diff_dir:
            diff_dir = diff_dir
        else:
            sanger_dir = "/mnt/ilustre/data/" + diff_dir
            tsanger_dir = "/mnt/ilustre/tsanger-data/" + diff_dir
            if list_dir(sanger_dir):
                diff_dir = sanger_dir
            elif list_dir(tsanger_dir):
                diff_dir = tsanger_dir
        return diff_dir

    def get_work_type(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        if not result:
            raise Exception('{}没有相应的任务id{}'.format(collection, task_id))
        return result['type']

    def metab_table_main_id(self, metab_table_id, task_id):
        metab_table_id = self.check_id(metab_table_id)
        # collection = self.db['metab_table']
        collection = self.db['exp']
        result = collection.find_one({"_id": metab_table_id, "task_id": task_id})
        #result = collection.find_one({"_id": metab_table_id})
        if not result:
            raise Exception('{}没有相应的id:{},{}！'.format(collection, metab_table_id, task_id))
        metab_table_main_id = result['main_id']
        return metab_table_main_id

    def check_id(self, myid):
        if isinstance(myid, types.StringTypes):
            myid = ObjectId(myid)
        elif isinstance(myid, ObjectId):
            myid = myid
        else:
            raise Exception("输入参数{}必须为字符串或者ObjectId类型!".format(myid))
        return myid

    def insert_main_id(self, collection, obj_id):
        return self.db[collection].update({"_id": ObjectId(obj_id)}, {'$set': {"main_id": ObjectId(obj_id)}})

    def insert_set_info(self, set_id, col_name, col_id):
        collection = self.db['metab_set']
        set_list = []
        if not isinstance(set_id, types.StringTypes):
            raise Exception("输入代谢集set_id参数必须为字符串类型！")
        set_list.extend([ObjectId(x) for x in set_id.split(",")])
        if isinstance(col_id, types.StringTypes):
            col_id = ObjectId(col_id)
        elif isinstance(col_id, ObjectId):
            pass
        else:
            raise Exception("输入col_id参数必须为字符串或者ObjectId类型")
        metab_set_id = "metab_set_id"
        try:
            for metab_set_id in set_list:
                result = collection.find_one({"_id": metab_set_id})
                result["is_use"] = 1
                collection.update({"_id": metab_set_id}, {"$set": result})
        except Exception:
            print("没有找到set_id:{}".format(metab_set_id))
        try:
            collection = self.db["col_name"]
            collection.find_one({"_id": col_id})
        except:
            "没有找到col_id: {} in {}".format(col_id, col_name)
        for metab_set_id in set_list:
            opts = {
                "set_id": metab_set_id,
                "col_name": col_name,
                "col_id": col_id
            }
            collection = self.db["metab_set_info"]
            collection.insert_one(opts)
        return True

    def delete_set(self, set_id):
        if isinstance(set_id, types.StringTypes):
            set_id = ObjectId(set_id)
        elif isinstance(set_id, ObjectId):
            pass
        else:
            raise Exception("输入set_id参数必须为字符串或者ObjectId类型")
        #collection = self.db["metab_set_info"]
        #results = collection.find({"set_id": set_id})
        #for result in results:
        #    col_name = result["col_name"]
        #    col_id = result["col_id"]
        metabset_r = self.db['metab_set'].find_one({"_id":set_id})
        if 'not_delete' in metabset_r:
            return True
        if metabset_r:
            task_id = metabset_r['task_id']
            status_r = self.db['sg_status'].find({"task_id":task_id,"metabset_id":{"$exists":1}})
            for status in status_r:
                if str(set_id) in  status['metabset_id'].split(','):

                    col_name = status['type_name']
                    col_id = status['table_id']
                    print "delete set record: col_name=%s, col_id=%s" % (col_name, col_id)
                    col = self.db[col_name]
                    try:
                        col_result = col.find_one({"_id": col_id})
                        col_result["params"] = ""
                        col.update({"_id": col_id}, {"$set": col_result})
                    except:
                        print("表{}:不能找到对应id: {}".format(col_name, col_id))

                    #更新任务状态表
                    self.db['sg_status'].update({"_id":status["_id"]},{"$set":{"status":"delete"}})

            #删除metab_set
            self.db["metab_set"].remove({"_id": set_id})
            #删除代谢集的detail表
            self.db['metab_set_detail'].remove({"set_id": set_id})

        return True

    '''
    def insert_main_table(self, collection, data):
        return self.db[collection].insert_one(SON(data)).inserted_id
    '''

    def get_hmdb_anno(self, task_id):
        collection = self.db['anno_hmdb']
        result = collection.find_one({"task_id": task_id, "name":"AnnoHmdb_Origin","status" :"end"})
        if not result:
            #raise Exception('{}没有相应的任务id{}'.format(collection, task_id))
            info = ""
        else:
            info = result['anno_hmdb']
        return info

    def get_kegg_species(self, task_id):
        collection = self.db["anno_keggp"]
        result = collection.find_one({"name" : "AnnoKeggp_Origin", "task_id": task_id})
        if not result:
            print '{}没有相应的id:{}！'.format(collection, task_id)
            return None
            #raise Exception('{}没有相应的id:{}！'.format(collection, task_id))
        params = result["params"]
        params_dict = json.loads(params)
        species = params_dict["organism"]
        return species

    def conmon_find_one(self,table_name,search):
        collection = self.db[table_name]
        r = collection.find_one(search)
        if r:
            return r
        else:
            return None

    def update_mongo(self,db_name,search, change):
        db = self.db[db_name]
        rets = db.find(search)
        if rets:
            for ret in rets:
                db.update({"_id":ret["_id"]},{"$set":change})
                return str(ret["_id"])
        else:
            return 'not find'

    def find_version_from_task(self, task_id):
        """
        根据task_id任务查询工作流版本号
        :param task_id:
        :return:
        """
        db = self.db["sg_task"]
        rets = db.find_one({"task_id": task_id})
        if rets:
            if "database" in rets:
                database = json.loads(rets['database'])
                version = database['kegg']
            else:
                version = ""
            return version
        else:
            return 'not find'
    
    # add by zhaoyuzhuo 20210927
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

    def get_main_info_by_record(self, main_table, **kwargs):
        """
        主表字段查询信息
        :param task_id:
        :return:
        """
        collection = self.db[main_table]
        result = collection.find_one(kwargs)
        return result
    
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

    # add by zhaoyuzhuo ,date:20211208
    def get_metab_desc1(self, task_id, type=None):
        collection = self.db['exp']
        result = collection.find_one({"task_id": task_id, "name" : "raw"})
        if not result:
            raise Exception('{}没有相应的task_id{}！'.format(collection, task_id))
        if type == "pos":
            table = os.path.join(result['table_path'].split(',')[0], "metab_desc.txt")
        elif type == "neg":
            table = os.path.join(result['table_path'].split(',')[1], "metab_desc.txt")
        elif type == "mix":
            table = os.path.join(result['table_path'].split(',')[2], "metab_desc.txt")
        elif not type:
            table = os.path.join(result['table_path'], "metab_desc.txt")
        else:
            raise Exception('metab_table类型type参数错误')
        if not "/mnt/ilustre/" in table:
            sanger_table = "/mnt/ilustre/data/" + table
            tsanger_table = "/mnt/ilustre/tsanger-data/" + table
            if exists(table):
                table = table
            elif exists(sanger_table):
                table = sanger_table
            elif exists(tsanger_table):
                table = tsanger_table
            else:
                raise Exception('metab_table path not exists-{}'.format(table))
        else:
            if not exists(table):
                raise Exception('metab_table path not exists-{}'.format(table))
        return table

    # add by zhaoyuzhuo ,date:20220126
    def get_scale_type(self, metab_table_id, task_id):
        metab_table_id = self.check_id(metab_table_id)
        collection = self.db['exp']
        result = collection.find_one({"_id": metab_table_id, "task_id": task_id})
        params = json.loads(result["params"])
        scale = params["scale"]
        if not result:
            raise Exception('{}没有相应的id{}！'.format(collection, metab_table_id))
        return scale
