# -*- coding: utf-8 -*-

"""跨mongo版本获取老项目的mongo连接, 基于biocluster.config 复制修改"""

import os
from pymongo import MongoClient
import grpc
from biocluster.proto import workflow_guide_pb2_grpc
from biocluster.proto import public_pb2, config_pb2, config_pb2_grpc
import traceback
import sys
import gevent
import web


class Config(object):

    def __init__(self):
        self._init_config_times = 0
        self.wfm_port = 7321
        self._mongodb_info = {}
        self._mongo_type_client = {}
        self._grpc_times = 0
        self.current_mode = "workflow"
        self.ntm_port = 7322
        self.DBVersion = None
        self.ProjectID = ""

    @property
    def mongo_client(self):
        mdbclienttype = self.get_mongo_client_type()
        return self._mongo_type_client[mdbclienttype]

    @property
    def biodb_mongo_client(self):
        mdbclienttype = self.get_mongo_client_type(ref=True)
        return self._mongo_type_client[mdbclienttype]

    def get_mongo_client_type(self, mtype="default", ref=False, db_version=0):
        db_projectsn = ""
        dydb = 0
        try:
            data = web.input()
            if hasattr(data, "db_version"):
                db_version = data.db_version
            if hasattr(data, "project_sn"):
                db_projectsn = data.db_projectsn
            if hasattr(data, "dydb"):
                dydb = data.dydb
        except:
            pass
        if self.DBVersion:
            db_version = self.DBVersion
        if dydb == 0:
            db_projectsn = ""
        if self.ProjectID !="":
            db_projectsn = self.ProjectID
        try:
            db_version = int(db_version)
        except:
            pass
        mdbinfotype= "{}_{}".format(mtype, db_version)
        if mdbinfotype not in self._mongodb_info.keys():
            self._get_project_db_config(mtype=mtype,db_version=db_version,project_id=db_projectsn)

        if ref:
            dburi = self._mongodb_info[mdbinfotype]["refuri"]
            mdbclienttype = "%s_ref_%s" % (mtype, db_version)
        else:
            dburi = self._mongodb_info[mdbinfotype]["uri"]
            mdbclienttype = "%s_%s" % (mtype, db_version)
        if mdbclienttype not in self._mongo_type_client.keys():
            self._mongo_type_client[mdbclienttype] = MongoClient(dburi, connect=False)

        return mdbclienttype


    def get_mongo_client(self, mtype="default", ref=False, db_version=0, task_id=""):
        """
        mdbinfotype：mongodbname info type  每个项目的参考库与数据库的uri与dbname的信息
        mdbclienttype: mongodb client type  每个项目的client连接，包含参考库与数据库
        :param mtype: 项目类型
        :param ref:  是否为参考库
        :param db_version: 数据库的版本，默认为0，可以设置1
        :return:
        """
        mdbclienttype = self.get_mongo_client_type(mtype,ref,db_version)
        return self._mongo_type_client[mdbclienttype]

    def get_mongo_dbname(self, mtype="default", ref=False, db_version=0, task_id=""):
        """
        获取数据库的名字
        :param mtype:
        :param ref:
        :param db_version:
        :return:
        """
        try:
            data = web.input()
            if hasattr(data, "db_version"):
                db_version = data.db_version
        except:
            pass
        if self.DBVersion:
            db_version = self.DBVersion
        try:
            db_version = int(db_version)
        except:
            pass
        mdbinfotype = "{}_{}".format(mtype, db_version)
        if mdbinfotype not in self._mongodb_info.keys():
            self.get_mongo_client_type(mtype,ref,db_version)
        if ref:
            return self._mongodb_info[mdbinfotype]["refdbname"]
        else:
            return self._mongodb_info[mdbinfotype]["dbname"]

    def _get_project_db_config(self, mtype="default", db_version=0,project_id=""):
        self._grpc_times += 1
        env_dist = os.environ
        port = self.wfm_port
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            self.current_mode = "tool"
        if self.current_mode == "tool":
            port = self.ntm_port
        print ("_get_project_db_config_port:{}".format(port))
        try:
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = config_pb2_grpc.ConfigServerStub(channel)
                response = stub.GetMongoDB(config_pb2.DBType(
                    type=mtype,
                    version=str(db_version),
                    projectid = project_id,
                ))
                cfg = {
                    "uri": response.uri,
                    "dbname": response.dbname,
                    "refuri": response.refuri,
                    "refdbname": response.refdbname,
                }
                mdbinfotype = "{}_{}".format(mtype, db_version)
                self._mongodb_info[mdbinfotype] = cfg
                # print "self._mongodb_info:{}".format(self._mongodb_info)
                self._grpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._grpc_times > 6:
                print("获取mongodb配置发生错误超过6次,退出运行: %s" % e)
                sys.stdout.flush()
                raise e
            elif self._grpc_times > 3:
                print("获取mongodb配置发生错误超过3次,换current_mode为tool重试: %s" % e)
                sys.stdout.flush()
                self.current_mode = "tool"
                gevent.sleep(10)
                self._get_project_db_config(mtype, db_version,project_id)
            else:
                print("获取mongodb配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self._get_project_db_config(mtype, db_version,project_id)
