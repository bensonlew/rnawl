# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
import json
from collections import OrderedDict
from biocluster.config import Config


class Base(object):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        self._config = Config()
        self._db = None
        self._ref_db = None
        self._project_type = self.__check_project_type()

    def __check_project_type(self):
        module_name = self.__module__
        mlist = module_name.split(".")
        mlist.pop()
        m = re.match(r'mainapp.models.mongo.(?:instant|submit).([\w_]+)', ".".join(mlist))
        if m:
            return m.group(1)
        else:
            return None

    @property
    def _client(self):
        return self._config.get_mongo_client(mtype=self._project_type)

    @property
    def _ref_client(self):
        return self._config.get_mongo_client(mtype=self._project_type, ref=True)

    @property
    def bind_object(self):
        return self._bind_object

    @property
    def db(self):
        if self._db is None:
            self._db = self._client[self._config.get_mongo_dbname(self._project_type)]
        return self._db

    @property
    def ref_db(self):
        if self._ref_db is None:
            self._ref_db = self._ref_client[self._config.get_mongo_dbname(self._project_type, True)]
        return self._ref_db

    @staticmethod
    def SortDict(param):
        if not isinstance(param, dict):
            raise Exception("传入的param不是一个字典")
        new_param = OrderedDict(sorted(param.items()))
        params = json.dumps(new_param)
        params = re.sub(':\s+', ':', params)
        params = re.sub(',\s+', ',', params)
        return params

    # def addSgStatus(self, objId, dbName, desc=None):
    #     collection = self.db[dbName]
    #     tableName = collection.find_one({"_id": objId})["name"]
    #     collection = self.db["sg_status"]
    #     insertData = {
    #         "table_id": objId,
    #         "table_name": tableName,
    #         "task_id": self.bind_object.taskId,
    #         "type_name": dbName,
    #         "status": "end",
    #         "is_new": "new",
    #         "desc": desc,
    #         "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    #         "params": self._params,
    #         "submit_location": self.bind_object.data.submit_location
    #     }
    #     collection.insert_one(insertData)


# class ApiManager(object):
#     """
#     api对象， 用于管理进行mongo数据库导入的各个api
#     """
#     def __init__(self, bind_object):
#         self._bind_object = bind_object
#         self._api_dict = dict()
#
#     def __getattr__(self, name, position="webroot"):
#         if name not in self._api_dict.keys():
#             self._api_dict[name] = self._get_api(name, position)
#         return self._api_dict[name]
#
#     def api(self, name, position="webroot"):
#         if name not in self._api_dict.keys():
#             self._api_dict[name] = self._get_api(name, position)
#         return self._api_dict[name]
#
#     def _get_api(self, name, position):
#         """
#         当position值为webroot的时候
#         用名字获取api对象, 名字应该从mainapp/models/mongo/instant 下一级开始， 比如名字为
#         meta.pan_core, 则模块的路径为mainapp/models/mongo/instant/meta/pan_core, 类名为PanCoreMongo
#         当position值为mbio的时候
#         用名字获取api对象, 名字应该从mbio/api/database 下一级开始， 比如名字为
#         meta.pan_core, 则模块的路径为mbio/api/database/pan_core, 类名为PanCore
#         """
#         if position == "webroot":
#             moduleName = "mainapp.models.mongo.instant." + name
#             name = re.split("\.", name)
#             className = name.pop(-1)
#             l = className.split("_")
#             l.append("Mongo")
#             l = [el.capitalize() for el in l]
#             className = "".join(l)
#             imp = importlib.import_module(moduleName)
#             obj = getattr(imp, className)(self._bind_object)
#             return obj
#         elif position == "mbio":
#             moduleName = "mbio.api.database.{}".format(name)
#             name = re.split("\.", name)[-1]
#             l = l = className.split("_")
#             l = [el.capitalize() for el in l]
#             className = "".join(l)
#             imp = importlib.import_module(name)
#             obj = getattr(imp, className)(self._bind_object)
#             return obj
#         else:
#             raise Exception("position的值应该为webroot或者mbio中的一个")
