# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
from bson import SON
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
from mainapp.models.mongo.core.base import Base
import random
import json

class RefRna(Meta):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(RefRna, self).__init__(bind_object=self._bind_object)
        self._project_type='ref_rna'
        #self.db_name = Config().MONGODB + '_ref_rna'
        #super(RefRna, self).__init__(db=self.db_name)

    def get_new_id(self, task_id, main_id):
        new_id = "%s_%s_%s" % (task_id, main_id[-4:], random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, main_id)
        return new_id

    def get_main_info(self, main_id, collection_name):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        main_info = collection.find_one({'_id': main_id})
        return main_info

    def insert_main_table(self, collection_name, data):
        return self.db[collection_name].insert_one(SON(data)).inserted_id
    
    def get_express_id(self, task_id, _type, express_method): #add by khl 
         """
         暂时实现的功能是根据表达量的软件如RSEM和表达量水平FPKM 获取表达量的主表(tsanger_ref_rna["sg_express"]) id
         :params _type: 表达量水平fpkm/tpm
         :params query_type: gene or transcript
         :params express_method: 表达量方法 featurecounts/rsem
         """
         #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
         collection=self.db["sg_express"]
         db=collection.find({"task_id":task_id})
         new_id = []
         for d in db:
                 _id = d["_id"]
                 #print d
                 params=json.loads(d['params'])
                 """
                 print params
                 print params['type']
                 print params['express_method']
                 print _type
                 print express_method
                 """
                 if _type == params['type'] and express_method == params['express_method'].lower():
                       return _id
                 else:
                       pass
         if not new_id:   #没有找到对应的表达量信息
             return False
         
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
    
    def get_class_code_id(self, task_id):
        """
        获得class_code的主表信息，为差异分析得到sequence_id的gene_name
        """
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = self.db["sg_express_class_code"]
        try:
            db=collection.find_one({"task_id":task_id})
            return db["_id"]
        except Exception:
            print "没有找到task_id: {}对应的class_code_id".format(task_id)
    
    def get_control_id(self,main_id):
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
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
                print compare_names
                return compare_names
            else:
                print "{}没有分组的对照信息！".format(str(main_id))
    
        except Exception:
            print "{}不存在".format(str(main_id))

    def insert_geneset_info(self, geneset_id, col_name, col_id):
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = self.db['sg_geneset']
        geneset_list = []
        if not isinstance(geneset_id, types.StringTypes):
            # geneset_id = ObjectId(geneset_id)
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
                result = collection.find_one({"_id":geneset_id})
                result["is_use"] = 1
                collection.update({"_id":geneset_id},{"$set":result})
        except Exception:
            print "没有找到geneset_id:{}".format(geneset_id)
        try:
            collection = self.db["col_name"]
            collection.find_one({"_id":col_id})
        except:
            "没有找到col_id:{} in {}".format(col_id, col_name)
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
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result


    def delete_geneset(self, geneset_id):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(data.geneset_id, ObjectId):
            pass
        else:
            raise Exception("输入geneset_id参数必须为字符串或者ObjectId类型!")
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = self.db['sg_geneset_info']
        results = collection.find({"geneset_id":geneset_id})
        for result in results:
            col_name = result["col_name"]
            col_id = result["col_id"]
            print col_id
            col = self.db[col_name]
            print col_name
            try:
                col_result = col.find_one({"_id":col_id})
                col_result["params"] = ""
                col.update({"_id":col_id}, {"$set":col_result})
            except:
                print "不能找到对应id{} in {}".format(col_id, col_name)
        collection = self.db["sg_geneset"]
        result = collection.find_one({"_id":geneset_id})
        if result:
            collection.remove({"_id":geneset_id})
        return True

    def check_assest_for_demo(self):
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = self.db['sg_task']
        nums = collection.count({"task_id":{"$regex": "refrna_demo_mouse.*"}})
        if nums:
            if nums <= 5:
                return False
            else:
                return True

    def insert_seq(self, mongo_data):
        task_id = mongo_data["task_id"]
        #mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = self.db['sg_query_seq']
        result = collection.find_one({"task_id": task_id})
        if result:
            collection.update({"task_id": task_id}, {"$set": mongo_data})
        else:
            collection.insert_one(mongo_data)

    def get_table_info_by_task_id(self, task_id, table_name="sg_task"):
        collection = self.db[table_name]
        result = collection.find_one({"task_id": task_id})
        return result


if __name__ == "__main__":
    data=RefRna()
    d = data.get_express_id("tsg_2000","fpkm","featurecounts")
    #print d
    data.get_control_id("5924f2a77f8b9a201d8b4567")

