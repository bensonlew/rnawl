# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from biocluster.config import Config
import importlib
import functools
from biocluster.core.function import get_clsname_form_path
import re
# import datetime
# from types import StringTypes
import types
import sys

PY3 = sys.version_info[0] == 3


class Base(object):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        self._config = Config()
        self._db = None
        self._ref_db = None
        self.manager = None
        self._project_type = self.__check_project_type()

    def __check_project_type(self):
        module_name = self.__module__
        mlist = module_name.split(".")
        mlist.pop()
        m = re.match(r'mbio.api.database.([\w_]+)', ".".join(mlist))
        if m:
            return m.group(1)
        else:
            return None

    def __del__(self):
        self._client.close()

    @property
    def project_type(self):
        return self._project_type

    @project_type.setter
    def project_type(self, value):
        self._project_type = value

    def set_project_type_as_workflow(self):
        if hasattr(self._config, "PROJECT_TYPE"):
            self._project_type = self._config.PROJECT_TYPE
        else:
            raise Exception("Workflow没有设置PROJECT_TYPE信息!")

    @property
    def _client(self):
        if self._config.DBVersion:
            return self._config.get_mongo_client(mtype=self._project_type, db_version=self._config.DBVersion)
        else:
            return self._config.get_mongo_client(mtype=self._project_type)

    @property
    def _ref_client(self):
        if self._config.DBVersion:
            return self._config.get_mongo_client(mtype=self._project_type, ref=True, db_version=self._config.DBVersion)
        else:
            return self._config.get_mongo_client(mtype=self._project_type, ref=True)

    @property
    def bind_object(self):
        return self._bind_object

    @property
    def db(self):
        print("testaa--------------------------")
        print(self._db)
        if self._db is None:
            if self._config.DBVersion:
                self._db = self._client[self._config.get_mongo_dbname(self._project_type,
                                                                      db_version=self._config.DBVersion)]
            else:
                self._db = self._client[self._config.get_mongo_dbname(self._project_type)]
        return self._db

    @property
    def ref_db(self):
        if self._ref_db is None:
            if self._config.DBVersion:
                self._ref_db = self._ref_client[self._config.get_mongo_dbname(self._project_type, True,
                                                                              db_version=self._config.DBVersion)]
            else:
                self._ref_db = self._ref_client[self._config.get_mongo_dbname(self._project_type, True)]
        return self._ref_db

    @property
    def path(self):
        class_name = str(type(self))
        m = re.match(r"<class\s\'(.*)\'>", class_name)
        class_name = m.group(1)
        paths = class_name.split(".")
        paths.pop()
        paths.pop(0)
        paths.pop(0)
        paths.pop(0)
        return ".".join(paths)

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


class ApiManager(object):
    """
    管理API对象
    """
    def __init__(self, bind_object, play_mod=False, debug=False):
        self._bind_object = bind_object
        self._api_dict = {}
        self._call_record = []
        self._play_mode = play_mod
        self.debug = debug

    def __getattr__(self, name):
        if name not in self._api_dict.keys():
            self._api_dict[name] = self._get_api(name)
        return self._api_dict[name]

    def api(self, name):
        if name not in self._api_dict.keys():
            self._api_dict[name] = self._get_api(name)
        return self._api_dict[name]

    def del_api(self, name):
        if name in self._api_dict.keys():
            del self._api_dict[name]

    @property
    def call_record(self):
        return self._call_record

    @property
    def play_mode(self):
        return self._play_mode

    def _get_api(self, name):
        """
        获取api对象

        :param name:
        :return:
        """
        class_name = get_clsname_form_path(name, tp="")
        module = importlib.import_module("mbio.api.database.%s" % name.lower())
        lib_obj = getattr(module, class_name)(self._bind_object)
        lib_obj.manager = self
        return lib_obj

    def add_call_record(self, api_name, func_name, args, kwargs):
        rc = CallRecord(api_name, func_name, args, kwargs, self)
        index = len(self._call_record)
        self._call_record.append(rc)
        return index

    def get_call_records_list(self):
        r_list = []
        for r in self._call_record:
            r_list.append(r.get_record())
        return r_list

    def load_call_records_list(self, record_list):
        self._call_record = []
        for r in record_list:
            cr = CallRecord(r["api"], r["func"], r["args"], r["kwargs"], self)
            self._call_record.append(cr)

    def play(self):
        i = 0
        for cr in self._call_record:
            if self.debug:
                self._bind_object.logger.debug("running #%s\t%s.%s(%s,%s) ..." %
                                               (i, cr.api_name, cr.func_name, ",".join(["%s" % a for a in cr.args]),
                                                ",".join(["%s=%s" % (key, value) for key, value in cr.kwargs.items()])))
            cr.run()
            i += 1


class CallRecord(object):
    """
    调用API截停记录
    """
    def __init__(self, api_name, func_name, args, kwargs, api_manager):
        self.api_manager = api_manager
        self.api_name = api_name
        self.func_name = func_name
        self.args = args
        self.kwargs = kwargs
        self._return_value = None

    @property
    def return_value(self):
        return self._return_value

    def run(self):
        p = re.compile(r"^\$(\d+)\$$")
        index = 0
        for key in self.args:
            if PY3:
                check = isinstance(key, str)
            else:
                check = isinstance(key, types.StringTypes)
            if check:
                match = p.match(key)
                if match:
                    call_index = int(match.group(1))
                    self.args[index] = self.api_manager.call_record[call_index].return_value
            index += 1
        for k, v in self.kwargs.items():
            if PY3:
                check = isinstance(v, str)
            else:
                check = isinstance(v, types.StringTypes)
            if check:
                match = p.match(v)
                if match:
                    call_index = int(match.group(1))
                    self.kwargs[k] = self.api_manager.call_record[call_index].return_value
        api = self.api_manager.api(self.api_name)
        func = getattr(api, self.func_name)
        self._return_value = func(*self.args, **self.kwargs)

    def get_record(self):
        record = {
            "api": self.api_name,
            "func": self.func_name,
            "args": self.args,
            "kwargs": self.kwargs
        }
        return record


def report_check(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        if args[0].manager.play_mode:
            return f(*args, **kwargs)
        if args[0].bind_object.IMPORT_REPORT_DATA is not True:
            return False
        else:
            if args[0].bind_object.IMPORT_REPORT_AFTER_END is True:
                list_args = list(args)
                api_obj = list_args.pop(0)
                index = args[0].manager.add_call_record(api_obj.path, f.__name__, list_args, kwargs)
                return "$%s$" % index
            else:
                return f(*args, **kwargs)
    return wrapper
