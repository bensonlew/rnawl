# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from mbio.packages.metagenomic.common import get_old_mongo
import functools

# client = Config().get_mongo_client(mtype="metagenomic")
# db = client[Config().get_mongo_dbname("metagenomic")]

client = None
db = None


def get_db(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        global client, db
        client = Config().get_mongo_client(mtype="metagenomic")
        db = client[Config().get_mongo_dbname("metagenomic")]
        print("id_covert db:{}".format(db))
        return func(*args, **kw)
    return wrapper


def get_mongo(db_version=None):
    if db_version == 0:
        client, db = get_old_mongo()
    else:
        client = Config().get_mongo_client(mtype="metagenomic")
        db = client[Config().get_mongo_dbname("metagenomic")]
    return client, db


class Test(object):
    def __init__(self):
        self.db = Config().get_mongo_client(mtype="metagenomic")[Config().get_mongo_dbname("metagenomic")]

    def runTest(self):
        raise Exception("self.db:{}".format(self.db))


# @get_db
def name2id(data, type=None, db_version=None):
    if type not in ["group", "raw", "task"]:
        raise Exception("type 类型错误：%s，必须为group/raw/task" % str(type))
    specimen_dic = {}
    if type == 'group':
        group, group_detail = get_group_collection(data, db_version)
        for group in group_detail:
            for id in group.keys():
                # specimen_dic[group[id]] = _get_objectid(id)
                specimen_dic[group[id]] = id
    elif type == 'raw':
        specimen_dic = get_raw_collection(data, db_version)[2]
    elif type == 'task':
        # print "### name2id ###"
        # print(get_raw_collection_from_task(data, db_version=db_version))
        # print "### name2id ###"
        specimen_dic = get_raw_collection_from_task(data, db_version=db_version)[2]
    return specimen_dic


# @get_db
def id2name(data, type=None, db_version=None):
    if type not in ["group", "raw", "task"]:
        raise Exception("type类型错误： %s, 必须为group/raw/task" % str(type))
    specimen_dic = {}
    if type == 'group':
        group, group_detail = get_group_collection(data, db_version)
        for group in group_detail:
            for id in group.keys():
                # specimen_dic[_get_objectid(id)] = group[id]
                specimen_dic[id] = group[id]
    elif type == 'raw':
        specimen_dic = get_raw_collection(data, db_version=db_version)[3]
    elif type == 'task':
        specimen_dic = get_raw_collection_from_task(data, db_version=db_version)[3]
    return specimen_dic


# @get_db
def get_raw_collection_from_task(task_id, check=True, db_version=None):
    db = get_mongo(db_version)[1]
    print("get_raw_collection_from_task db:{}, task_id: {}".format(db, task_id))
    raw_stat = db['data_stat']
    if check:   # task id 类似这种形式的sanger_160815_24551 将查不到raw data表，需对task id做处理
        sp = task_id.split('_')
        if len(sp) >2:
            task_id = '_'.join([sp[0], sp[1]])
    raw_stat_find = raw_stat.find_one({"task_id": task_id, "type": "raw"})
    # print raw_stat_find[u"_id"]
    # print(get_raw_collection(raw_stat_find["_id"], db_version))
    return get_raw_collection(raw_stat_find["_id"], db_version)


# @get_db
def get_raw_collection(data, db_version=None):  # data_stat id
    specimen_list = []
    specimen_id_list = []
    specimen_dic = {}
    specimen_id_dic = {}
    db = get_mongo(db_version)[1]
    data = _get_objectid(data)
    raw = db['data_stat_detail']
    raw_finds = raw.find({"data_stat_id": data})
    if not raw_finds or raw_finds.count() == 0:
        raise Exception('data_stat原始序列id:%s在data_stat_detail表中没有找到' % str(data))
    for find in raw_finds:
        specimen_list.append(find["specimen_name"])
        if "origin_id" in find.keys():
            specimen_id_list.append(str(find["origin_id"]))
            specimen_id_dic[str(find["origin_id"])] = find["specimen_name"]
            specimen_dic[find["specimen_name"]] = str(find["origin_id"])
        else:
            specimen_id_list.append(str(find["_id"]))
            specimen_id_dic[str(find["_id"])] = find["specimen_name"]
            specimen_dic[find["specimen_name"]] = str(find["_id"])
    # print("get_raw_collection return: 1:{}, 2{}, 3:{}, 4:{}".format(specimen_list, specimen_id_list, specimen_dic, specimen_id_dic))
    return specimen_list, specimen_id_list, specimen_dic, specimen_id_dic


# @get_db
def get_group_collection(data, db_version=None):  # data is specimen_group _id
    db = get_mongo(db_version)[1]
    group_table = db['specimen_group']
    data = _get_objectid(data)
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    if not group_schema or len(group_schema) == 0:
        raise Exception("无法根据传入的group_id:{}在specimen_group表里找到相应的记录".format(data))
    # return  group_schema.keys()
    return group_schema["category_names"], group_schema["specimen_names"]
    # return group_schema.category_names, group_schema.specimen_names


def _get_objectid(data):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
        else:
            try:
                data = ObjectId(data)
            except:
                raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
    return data


def test():
    # data1, data2 = get_group_collection("59f9888ec81e103f648b4570")
    # data = name2id("whole_api", type='task')
    # data = name2id(ObjectId("59fc4ad3a4e1af6371847a27"), type='raw')
    # data = name2id("59fc4ad3a4e1af6371847a27", type="raw")
    # data = name2id(ObjectId("59f9888ec81e103f648b4570"), type="group")
    # data = name2id("59f9888ec81e103f648b4570", type="group")
    data = id2name("whole_api", type='task')
    # data = id2name(ObjectId("59fc4ad3a4e1af6371847a27"), type='raw')
    # data = id2name("59f9888ec81e103f648b4570", type="group")
    # out = data.keys()
    # print type(out)
    # print out
    print data
    print data[ObjectId('59f7e488a4e1af151a15d787')]
    # print data["LD31"]
    # print data


if __name__ == "__main__":
    test()
