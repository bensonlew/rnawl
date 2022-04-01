# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from .core.base import Base
from bson.objectid import ObjectId
import types
from bson import SON
import re
import datetime
import json
from collections import OrderedDict
from types import StringTypes


class Metaasv(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Metaasv, self).__init__(self._bind_object)
        self._project_type = "metaasv"
        self.level = {
            9: "asv", 8: "s__", 7: "g__", 6: "f__", 5: "o__",
            4: "c__", 3: "p__", 2: "k__", 1: "d__"
        }


    def get_otu_table_info(self, asv_id):

        if isinstance(asv_id, types.StringTypes):
            asv_id = ObjectId(asv_id)
        elif isinstance(asv_id, ObjectId):
            asv_id = asv_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['asv']
        result = collection.find_one({"_id": asv_id})
        return result

    def get_task_info(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def insert_main_table(self, collection, data):
        return self.db[collection].insert_one(SON(data)).inserted_id

    def insert_none_table(self, collection):
        return self.db[collection].insert_one({}).inserted_id

    def insert_main_table_new(self, collection, obj_id, data):
        # return self.db[collection].find_one_and_update({"_id": ObjectId(obj_id)}, {'$set': data}, upsert=True)
        return self.db[collection].update({"_id": ObjectId(obj_id)}, {'$set': data}, upsert=True)

    def update_status_failed(self, collection, doc_id):
        """
        改特定_id主表的status状态从start为failed，主要用于特殊投递任务失败

        params collection: 主表collection名称
        params doc_id: 主表_id
        """
        self.db[collection].update_one({'_id': ObjectId(doc_id), "status": "start"}, {"$set": {'status': 'failed'}})

    def update_workflow_id(self, collection, main_id, workflow_id):
        """
        """
        self.db.workflowid2analysisid.insert_one({'workflow_id': workflow_id,
                                                  "main_id": ObjectId(main_id),
                                                  'collection': collection,
                                                  "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})

    def sampleIdToName(self, task_id, sampleIds):
        """
        将一个用逗号隔开的样本ID的集合转换成样本名，返回一个用逗号隔开的样本名的集合
        """
        myIds = re.split("\s*,\s*", sampleIds)
        collection = self.db["specimen"]
        main_result = collection.find_one({"task_id": task_id})
        main_id = main_result["_id"]
        collection_detail = self.db["specimen_detail"]
        mySampleNames = list()
        for id_ in myIds:
            if id_ == "":
                raise Exception("存在空的sample_id")
            if not isinstance(id_, ObjectId):
                if isinstance(id_, types.StringTypes):
                    id_ = ObjectId(id_)
                else:
                    raise Exception("样本id必须为ObjectId对象或者其对应的字符串！")
            result = collection_detail.find_one({"specimen_id": main_id,"_id": id_})
            if not result:
                raise Exception("无法根据传入的_id:{}在sg_speciem表里找到相应的记录".format(str(id_)))
            mySampleNames.append(result["specimen"])
        mySamples = ",".join(mySampleNames)
        return mySamples

    def get_main_info(self, main_id, collection_name):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        express_info = collection.find_one({'_id': main_id})
        return express_info

    def get_mongo_common(self,db_name, search):
        collection = self.db[db_name]
        result = collection.find_one(search)
        if result:
            return result
        else:
            return False

    def filter_json_sort(filter_detail):
        filters = json.loads(filter_detail)
        temp = []
        for i in filters:
            if isinstance(i, dict):
                temp.append(OrderedDict(sorted(i.items(), key=lambda t: t[0])))
        return temp

    def get_group_name(self, group_id, lefse=False, second_group=''):
        """
        根据分组方案id获取分组方案名字
        :param group_id: 分组方案id
        :return: 分组方案名字
        """
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception("group_detail必须为ObjectId对象或其对应的字符串!")
        collection = self.db['specimen_group']
        result = collection.find_one({'_id': group_id})
        gname = result['group_name']
        if lefse and second_group:
            gname = gname + ',' + 'second_group'
        return gname

    def get_otu_sample_name(self, otu_id):
        """
        获取otu表样本名字
        """
        if isinstance(otu_id, StringTypes):
            otu_id = ObjectId(otu_id)
        elif isinstance(otu_id, ObjectId):
            otu_id = otu_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        collection_2 = self.db['asv_specimen']
        result_2 = collection_2.find({"asv_id": otu_id})
        collection_3 = self.db['specimen_detail']
        sample_name = []
        for i in result_2:
            specimen_id = ObjectId(i['specimen_id'])
            result_3 = collection_3.find_one({"_id": specimen_id})
            sample_name.append(result_3["specimen"])
        return sample_name

    def get_diversity_table_info(self, alpha_diversity_id):

        if isinstance(alpha_diversity_id, types.StringTypes):
            asv_id = ObjectId(alpha_diversity_id)
        elif isinstance(alpha_diversity_id, ObjectId):
            asv_id = alpha_diversity_id
        else:
            raise Exception("输入alpha_diversity_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['alpha_diversity']
        result = collection.find_one({"_id": asv_id})
        return result

    def convert_json(self, filter_json):
        filter_json = json.loads(filter_json)
        temp = []
        for key in filter_json.keys():
            key_value = filter_json[key]
            key_dict = {}
            if key not in temp:
                key_dict[key] = key_value

                temp.append(key_dict)

    def get_otu_table_path_from_asv_id(self, asv_id):
        """
        根据asv_id从MongoDB中获取相对丰度表
        :param asv_id:
        :return:
        """
        if isinstance(asv_id, types.StringTypes):
            asv_id = ObjectId(asv_id)
        elif isinstance(asv_id, ObjectId):
            asv_id = asv_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['asv']
        result = collection.find_one({"_id": asv_id})
        if result:
            asv_path = result["asv_path"]
            return asv_path

    def get_asv_id_from_table(self, table_name, table_id):
        """
        主要是根据table_name和table_id获取asv_id
        :param table_name: 主表名称
        :param table_id: 主表id
        :return:
        """
        main_collection = self.db[table_name]
        main_id = ObjectId(table_id)
        try:
            result = main_collection.find_one({"_id": main_id})
            if result:
                asv_id = result["asv_id"]
                return asv_id
        except:
            raise Exception("未从主表{}中查到对应的信息{}".format(table_name, table_id))

    def export_otu_table_by_level(self, otuId, targetPath, level=9):
        """
        用于asv表的实时接口--生成MongoDB的数据
        :param otuId:
        :param targetPath:
        :param level:
        :return:
        """
        self._try_write(targetPath)
        level = int(level)
        print "正在导出级别为{}:{}的ASV表格:{}".format(level, self.level[level], targetPath)
        collection = self.db['asv_specimen']
        results = collection.find({"asv_id": ObjectId(otuId)})
        if not results.count():
            raise Exception("asv_id: {}在asv_specimen表中未找到！".format(otuId))
        samples = list()
        for result in results:
            if "specimen_id" not in result:
                raise Exception("asv_id:{}错误，请使用新导入的asv表的id".format(otuId))
            sp_id = result['specimen_id']
            my_collection = self.db['specimen_detail']
            my_result = my_collection.find_one({"_id": sp_id})
            if not my_result:
                raise Exception("意外错误，样本id:{}在specimen_detail表里未找到".format(sp_id))
            samples.append(my_result["specimen"])
        if "total" not in samples:
            samples.append("total")
        if "percent" not in samples:
            samples.append("percent")
        collection = self.db['asv_detail_level']
        name_dic = dict()
        ##这里查询的表格变了，主要是因为页面展示需要percent和total两列信息，这样直接根据工作流计算的结果在此基础上做累加的算法，其他不做改变
        results = collection.find({"asv_id": ObjectId(otuId), "level_id": 9})
        if not results.count():
            raise Exception("asv_id: {}在asv_detail_level表中未找到！".format(otuId))
        for col in results:## 将各水平相同的level水平进行合并和样本合并
            tmp = level + 1
            new_classify_name = self._create_classify_name(col, tmp)
            if new_classify_name not in name_dic:
                name_dic[new_classify_name] = dict()
                for sp in samples:
                    name_dic[new_classify_name][sp] = float(col[sp])
            else:
                for sp in samples:
                    name_dic[new_classify_name][sp] += float(col[sp])
        with open(targetPath, "wb") as f:
            f.write("ASV ID\t%s\n" % "\t".join(samples))
            for k in name_dic.iterkeys():
                line = k
                for s in samples:
                    line += "\t" + str(name_dic[k][s])
                line += "\n"
                f.write(line)
        return targetPath

    def _create_classify_name(self, col, tmp):
        for i in range(1, 10):
            if self.level[i] not in col:
                print(col)
                raise Exception("Mongo数据库中的taxonomy信息不完整, 缺少{}".format(self.level[i]))
        new_col = list()
        for i in range(1, tmp):
            new_col.append(col[self.level[i]])
        return "; ".join(new_col)

    def _try_write(self, targetPath):
        try:
            with open(targetPath, "wb"):
                pass
        except Exception as e:
            raise Exception(e)
        return True

    def add_otu_detail(self, otu_path, from_otu_table, level):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            raise Exception("无法根据传入的_id:{}在asv表里找到相应的记录".format(str(from_otu_table)))
        project_sn = result['project_sn']
        task_id = result['task_id']
        covered_level = list()
        if "select_id" in result:
            covered_level = result["select_id"]
            covered_level = json.loads(covered_level)
            if level not in covered_level:
                covered_level.append(int(level))
        else:
            covered_level.append(9)
            if level not in covered_level:
                covered_level.append(int(level))
        covered_level = list(set(covered_level))
        covered_level.sort()
        result["select_id"] = json.dumps(covered_level)
        collection.update({"_id": from_otu_table}, {"$set": result}, upsert=False)
        insert_data = list()
        with open(otu_path, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split('\t', head)
            new_head = head[1:]
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                otu_detail['asv_id'] = from_otu_table
                otu_detail['project_sn'] = project_sn
                otu_detail['task_id'] = task_id
                otu_detail["level_id"] = int(level)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = int(float(sample_num[i]))
                otu_detail["total"] = int(float(line[-2]))
                otu_detail["percent"] = float(line[-1])
                insert_data.append(otu_detail)
        try:
            collection = self.db['asv_detail_level']
            collection.insert_many(insert_data)
        except Exception as e:
            raise Exception("导入asv_detail_level表格失败：{}".format(e))
        else:
            print "导入asv_detail_level表格成功"

    def find_asv_name(self, task_id):
        """

        :param task_id:
        :return:
        """
        collection = self.db["asv_set"]
        asv_name_list = []
        results = collection.find({"task_id": task_id})
        if results:
            for result in results:
                asv_name = result["name"]
                if asv_name not in asv_name_list:
                    asv_name_list.append(asv_name)
        return asv_name_list

    def get_bugbase_info(self, bugbase_id):

        if isinstance(bugbase_id, types.StringTypes):
            bugbase_id = ObjectId(bugbase_id)
        elif isinstance(bugbase_id, ObjectId):
            bugbase_id = bugbase_id
        else:
            raise Exception("输入bugbase_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['bugbase']
        result = collection.find_one({"_id": bugbase_id})
        return result

    def get_faprotax_info(self, faprotax_id):

        if isinstance(faprotax_id, types.StringTypes):
            faprotax_id = ObjectId(faprotax_id)
        elif isinstance(faprotax_id, ObjectId):
            faprotax_id = faprotax_id
        else:
            raise Exception("输入faprotax_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['faprotax']
        result = collection.find_one({"_id": faprotax_id})
        return result