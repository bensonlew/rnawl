# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""系统配置，获取ip端口、日志、数据库"""

import ConfigParser
import socket
# import random
import os
from .core.singleton import singleton
# import struct
# import platform
import re
# import importlib
# import web
# import subprocess
from pymongo import MongoClient
# import subprocess
# from IPy import IP
import boto.s3.connection
import boto
# import threading
import grpc
from .proto import workflow_guide_pb2_grpc
from .proto import public_pb2, config_pb2, config_pb2_grpc
import traceback
import sys
import gevent
import web

# web.config.debug = False


@singleton
class ConfigTrans(object):

    def __init__(self):
        self.rcf = ConfigParser.RawConfigParser()
        self.rcf.read(os.path.dirname(os.path.realpath(__file__)) + "/main.conf")
        self.PROJECT_TYPE = ""
        self.JOB_PLATFORM = ""
        self.JOB_QUEUE = ""
        self.WORK_DIR = self.rcf.get("Command", "work_dir")
        self.SOFTWARE_DIR = self.rcf.get("Command", "software_dir")
        self._workdir = self.rcf.get("Command", "work_dir")
        self.wpm_user = self.rcf.get("Command", "user")
        self.SCRIPT_DIR = ""
        self.PACKAGE_DIR = ""
        self._init_config_times = 0
        self.wfm_port = 7321
        self._mongodb_info = {}
        self._mongo_type_client = {}
        self._webauth_info = {}
        self.current_workflow_id = ""
        self._rgw_conn = {}
        self._rgw_account = {}
        self._grpc_times = 0
        self._project_region_bucket = {}
        self._path_region_bucket = {}
        self.current_tool_id = ""
        self.current_mode = "workflow"
        self.ntm_port = 7322
        self.WFM_INSTANT_TIMEOUT = 600  # 设置即时任务通过grpc请求发起任务运行超时时间
        self.WFM_SUBMIT_TIMEOUT = 300   # 设置投递任务通过grpc请求发起任务响应超时
        self.RGW_ENABLE = True  # 是否将文件传输到对象存储上去
        self.get_workdir_time = 0
        self.DBVersion = None
        # self.init_config(project_type="", workflow_id="", grpc_port=self.ntm_port)

    def init_config(self, project_type, workflow_id, wfm_grpc_port=7321, instant=False, ntm_grpc_port=7322):
        self._init_config_times += 1
        self.PROJECT_TYPE = project_type
        self.current_workflow_id = workflow_id
        if wfm_grpc_port:
            self.wfm_port = wfm_grpc_port
        else:
            self.wfm_port = 7321
        self.ntm_port = ntm_grpc_port if ntm_grpc_port else 7322
        try:
            with grpc.insecure_channel('localhost:%s' % wfm_grpc_port) as channel:
                stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
                response = stub.GetRunInfo(public_pb2.Workflow(instant=instant, workflow_id=workflow_id,
                                                               process_id=int(os.getpid())))
                self.JOB_PLATFORM = response.platform
                self.JOB_QUEUE = response.defaut_queue
                # self.WORK_DIR = "/mnt/ilustre/users/sanger-dev/bc2/workspace/"
                self.WORK_DIR = response.workspace
                # self.SOFTWARE_DIR = response.software_dir
                self.SCRIPT_DIR = response.script_dir
                self.PACKAGE_DIR = response.package_dir
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._init_config_times > 6:
                print("获取配置发生错误超过6次,退出运行: %s" % e)
                sys.stdout.flush()
                raise e
            elif self._init_config_times > 3:
                print("获取配置发生错误超过3次,换grpc_port重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self.init_config(project_type, workflow_id, 7321, instant, 7322)
            # if self._init_config_times > 3:
            #     print("获取配置发生错误超过3次,退出运行: %s" % e)
            #     raise e
            else:
                print("获取配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self.init_config(project_type, workflow_id, wfm_grpc_port, instant, ntm_grpc_port)

    def init_current_tool_id(self, tool_id):
        self.current_tool_id = tool_id

    @property
    def mongo_client(self):
        db_version = 0
        if self.DBVersion:
            db_version = self.DBVersion
        else:
            try:
                data = web.input()
                if hasattr(data, "db_version"):
                    db_version = data.db_version
            except:
                pass
        try:
            db_version = int(db_version)
        except:
            pass
        print "config mongo_client dbversion:{}".format(db_version)
        mdbclienttype = "{}_{}".format("default", db_version)
        mdbinfotype = "{}_{}".format("default", db_version)
        self._get_project_db_config(db_version=db_version)
        dburi = self._mongodb_info[mdbinfotype]["uri"]
        if mdbclienttype not in self._mongo_type_client.keys():
            self._mongo_type_client[mdbclienttype] = MongoClient(dburi, connect=False)
        return self._mongo_type_client[mdbclienttype]

    @property
    def biodb_mongo_client(self):
        db_version = 0
        if self.DBVersion:
            db_version = self.DBVersion
        else:
            try:
                data = web.input()
                if hasattr(data, "db_version"):
                    print "webinput dbversion:{}".format(data.db_version)
                    db_version = data.db_version
            except:
                pass
        try:
            db_version = int(db_version)
        except:
            pass
        print "config biodb_mongo_client dbversion:{}".format(db_version)
        mdbclienttype = "{}_ref_{}".format("default", db_version)
        mdbinfotype = "{}_{}".format("default", db_version)
        self._get_project_db_config(db_version=db_version)
        dburi = self._mongodb_info[mdbinfotype]["refuri"]
        if mdbclienttype not in self._mongo_type_client.keys():
            self._mongo_type_client[mdbclienttype] = MongoClient(dburi, connect=False)
        return self._mongo_type_client[mdbclienttype]

    def get_mongo_client(self, mtype="default", ref=False, db_version=0):
        """
        mdbinfotype：mongodbname info type  每个项目的参考库与数据库的uri与dbname的信息
        mdbclienttype: mongodb client type  每个项目的client连接，包含参考库与数据库
        :param mtype: 项目类型
        :param ref:  是否为参考库
        :param db_version: 数据库的版本，默认为0，可以设置1
        :return:
        """
        # if mtype:
        #     if ref:
        #         type_key = "%s_ref" % mtype
        #     else:
        #         type_key = mtype
        #     if type_key not in self._mongo_type_client.keys():
        #         if self.rcf.has_option("MONGO", "%s_uri" % type_key):
        #             self._mongo_type_client[type_key] = MongoClient(self.rcf.get("MONGO", "%s_uri" % type_key),
        #                                                             connect=False)
        #         else:
        #             if ref:
        #                 if self.rcf.has_option("MONGO", "%s_uri" % mtype):
        #                     self._mongo_type_client[type_key] = self.get_mongo_client(mtype=mtype)
        #                 else:
        #                     self._mongo_type_client[type_key] = self.biodb_mongo_client
        #             else:
        #                 self._mongo_type_client[type_key] = self.mongo_client
        #     return self._mongo_type_client[type_key]
        # else:
        #     if not ref:
        #         return self.mongo_client
        #     else:
        #         return self.biodb_mongo_client
        print "DBVersion is:{}".format(self.DBVersion)
        if self.DBVersion:
            print "this way"
            db_version = self.DBVersion
        else:
            try:
                data = web.input()
                if hasattr(data, "db_version"):
                    print "webinput dbversion:{}".format(data.db_version)
                    db_version = data.db_version
            except:
                pass
        try:
            db_version = int(db_version)
        except:
            pass
        print "config dbversion:{}".format(db_version)
        # mongodb info type 全局变量中存储mongo的数据库信息
        mdbinfotype= "{}_{}".format(mtype, db_version)
        if mdbinfotype not in self._mongodb_info.keys():
            self._get_project_db_config(mtype, db_version)
        if ref:
            dburi = self._mongodb_info[mdbinfotype]["refuri"]
            mdbclienttype = "%s_ref_%s" % (mtype, db_version)
        else:
            dburi = self._mongodb_info[mdbinfotype]["uri"]
            mdbclienttype = "%s_%s" % (mtype, db_version)
        if mdbclienttype not in self._mongo_type_client.keys():
            self._mongo_type_client[mdbclienttype] = MongoClient(dburi, connect=False)

        return self._mongo_type_client[mdbclienttype]

    def get_mongo_dbname(self, mtype="default", ref=False, db_version=0):
        """
        获取数据库的名字
        :param mtype:
        :param ref:
        :param db_version:
        :return:
        """
        if self.DBVersion:
            db_version = self.DBVersion
        else:
            try:
                data = web.input()
                if hasattr(data, "db_version"):
                    db_version = data.db_version
            except:
                pass
        try:
            db_version = int(db_version)
        except:
            pass
        mdbinfotype = "{}_{}".format(mtype, db_version)
        if mdbinfotype not in self._mongodb_info.keys():
            self._get_project_db_config(mtype, db_version)

        if ref:
            return self._mongodb_info[mdbinfotype]["refdbname"]
        else:
            return self._mongodb_info[mdbinfotype]["dbname"]

    def _get_project_db_config(self, mtype="default", db_version=0,port=""):
        self._grpc_times += 1
        # env_dist = os.environ
        # if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
        #     self.current_mode = "tool"
        # if self.current_mode == "tool":
        #     port = self.ntm_port
        # else:
        #     port = self.wfm_port
        if port == "":
            port = self.wfm_port
        print "_get_project_db_config_port:{}".format(port)
        try:
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = config_pb2_grpc.ConfigServerStub(channel)
                response = stub.GetMongoDB(config_pb2.DBType(
                    type=mtype,
                    version=str(db_version),
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
                self._get_project_db_config(mtype, db_version,self.ntm_port)
            else:
                print("获取mongodb配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self._get_project_db_config(mtype, db_version)

    def get_webauth(self, ptype = "tool_lab"):
        self._get_webauth_config(ptype)
        return self._webauth_info[ptype]

    def _get_webauth_config(self, ptype="tool_lab",port=""):
        self._grpc_times += 1
        # env_dist = os.environ
        # if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
        #     self.current_mode = "tool"
        # if self.current_mode == "tool":
        #     port = self.ntm_port
        # else:
        #     port = self.wfm_port
        if port == "":
            port = self.wfm_port
        print "get_webauth_config_port:{}".format(port)
        hostname = socket.gethostname()
        ip = socket.gethostbyname(hostname)
        # print("获取webauth:hostname %s, ip %s,port %s",hostname,ip,port)
        print "获取webauth:hostname {}, ip {},port {}".format(hostname,ip,port)
        try:
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = config_pb2_grpc.ConfigServerStub(channel)
                response = stub.GetWebAuth(config_pb2.ProjectType(
                    type=ptype,
                ))
                # authkey: 458b97de4c0bb5bf416c8cea208309ed
                # client: client03
                # binds_id: 5e8c0a091b1800007c006b1a
                # interface_id: 1348
                # env_name: offline
                cfg = {
                    "authkey": response.key,
                    "client": response.client,
                    "binds_id": response.bindsid,
                    "env_name": response.envname,
                    "interface_id":response.interfaceid,
                    "url":response.url,
                }
                self._webauth_info[ptype] = cfg
                self._grpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._grpc_times > 6:
                print("获取webauth配置发生错误超过6次,退出运行: %s" % e)
                sys.stdout.flush()
                raise e
            elif self._grpc_times > 3:
                print("获取webauth配置发生错误超过3次,换current_mode为tool重试: %s" % e)
                sys.stdout.flush()
                self.current_mode = "tool"
                gevent.sleep(10)
                self._get_webauth_config(ptype,self.ntm_port)
            else:
                print("获取webauth配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self._get_webauth_config(ptype)


    def get_rgw_conn(self, region="s3", bucket="sanger", new=False):
        key = "%s_%s" % (region, bucket)
        if key not in self._rgw_account.keys():
            self._grpc_times += 1
            try:
                env_dist = os.environ
                if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
                    self.current_mode = "tool"
                if self.current_mode == "tool":
                    port = self.ntm_port
                else:
                    port = self.wfm_port
                with grpc.insecure_channel('localhost:%s' % port) as channel:
                    stub = config_pb2_grpc.ConfigServerStub(channel)
                    response = stub.GetRgwAccount(config_pb2.BucketInfo(
                        region=region,
                        bucket=bucket,
                    ))
                    self._rgw_account[key] = response
                    self._grpc_times = 0
            except Exception as e:
                exstr = traceback.format_exc()
                print(exstr)
                sys.stdout.flush()
                if self._grpc_times > 3:
                    print("获取rgw账户发生错误超过3次,退出运行: %s" % e)
                    sys.stdout.flush()
                    raise e
                else:
                    print("获取rgw账户发生错误10秒后重试: %s" % e)
                    sys.stdout.flush()
                    gevent.sleep(10)
                    self.get_rgw_conn(region, bucket)
        if new:
            conn = boto.connect_s3(
                aws_access_key_id=self._rgw_account[key].accesskey,
                aws_secret_access_key=self._rgw_account[key].secretkey,
                host=self._rgw_account[key].host,
                port=self._rgw_account[key].port,
                is_secure=self._rgw_account[key].issecure,
                # uncomment if you are not using ssl
                calling_format=boto.s3.connection.OrdinaryCallingFormat(),
            )
            self._rgw_conn[key] = conn
            return conn
        else:
            if key not in self._rgw_conn.keys():
                self._rgw_conn[key] = boto.connect_s3(
                    aws_access_key_id=self._rgw_account[key].accesskey,
                    aws_secret_access_key=self._rgw_account[key].secretkey,
                    host=self._rgw_account[key].host,
                    port=self._rgw_account[key].port,
                    is_secure=self._rgw_account[key].issecure,
                    # uncomment if you are not using ssl
                    calling_format=boto.s3.connection.OrdinaryCallingFormat(),
                )
            return self._rgw_conn[key]

    # def get_rgw_download_chunk_size(self):
    #     size = self.rcf.get("RGW", "download_chunk_size")
    #     return get_from_friendly_size(size)
    #
    # def get_http_download_chunk_size(self):
    #     size = self.rcf.get("HTTP", "download_chunk_size")
    #     return get_from_friendly_size(size)
    #
    # def get_rgw_upload_chunk_size(self):
    #     size = self.rcf.get("RGW", "upload_chunk_size")
    #     return get_from_friendly_size(size)
    #
    # def get_rgw_min_size_to_split(self):
    #     size = self.rcf.get("RGW", "min_size_to_split")
    #     return get_from_friendly_size(size)
    #
    # def get_http_min_size_to_split(self):
    #     size = self.rcf.get("HTTP", "min_size_to_split")
    #     return get_from_friendly_size(size)
    #
    # def get_rgw_max_threads(self):
    #     return int(self.rcf.get("RGW", "max_threads"))

    def get_project_region_bucket(self, project_type="default"):
        # if self.rcf.has_option("RGW", "%s_bucket" % project_type):
        #     return self.rcf.get("RGW", "%s_bucket" % project_type)
        # else:
        #     return self.rcf.get("RGW", "default_bucket")
        if project_type in self._project_region_bucket.keys():
            return self._project_region_bucket[project_type]
        self._grpc_times += 1
        try:
            env_dist = os.environ
            if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
                self.current_mode = "tool"
            if self.current_mode == "tool":
                port = self.ntm_port
            else:
                port = self.wfm_port
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = config_pb2_grpc.ConfigServerStub(channel)
                response = stub.GetProjectBucket(config_pb2.ProjectType(
                    type=project_type,
                ))
                self._project_region_bucket[project_type] = response.url
                self._grpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._grpc_times > 3:
                print("获取project bucket配置发生错误超过3次,退出运行: %s" % e)
                sys.stdout.flush()
                raise e
            else:
                print("获取project bucket配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self.get_project_region_bucket(project_type)
        return self._project_region_bucket[project_type]

    def get_bucket_from_path(self, path):
        # for k, v in self.rcf.items("PROJECT_MODULE"):
        #     if v:
        #         p_list = re.split(r'\s*\|\s*', v)
        #         for p in p_list:
        #             m = re.match(r"^regexp\((.+)\)$", p)
        #             if m:
        #                 exp = m.group(1)
        #                 m1 = re.match(exp, path)
        #                 if m1:
        #                     return self.get_project_region_bucket(k)
        #             else:
        #                 if p == path:
        #                     return self.get_project_region_bucket(k)
        # return self.get_project_region_bucket()
        if path in self._path_region_bucket.keys():
            return self._path_region_bucket[path]
        self._grpc_times += 1
        try:
            env_dist = os.environ
            if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
                self.current_mode = "tool"
            if self.current_mode == "tool":
                port = self.ntm_port
            else:
                port = self.wfm_port
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = config_pb2_grpc.ConfigServerStub(channel)
                response = stub.GetBucketFromPath(config_pb2.Path(
                    path=path,
                ))
                print "url:{}".format(response)
                self._path_region_bucket[path] = response.url
                self._grpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._grpc_times > 3:
                print("获取path bucket配置发生错误超过3次,退出运行: %s" % e)
                sys.stdout.flush()
                raise e
            else:
                print("获取path bucket配置发生错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(10)
                self.get_bucket_from_path(path)
        return self._path_region_bucket[path]

    # @property
    # def LISTEN_IP(self):
    #     if self._listen_ip is None:
    #         self._listen_ip = self.get_listen_ip()
    #     return self._listen_ip
    #
    # def get_listen_ip(self):
    #     """
    #     获取配置文件中IP列表与本机匹配的IP作为本机监听地址
    #     """
    #     # def getip(ethname):
    #     #     s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    #     #     fcn = importlib.import_module("fcntl")
    #     #     return socket.inet_ntoa(fcn.ioctl(s.fileno(), 0X8915, struct.pack("256s", ethname[:15]))[20:24])
    #     # set_ipl_ist = self.rcf.get("Network", "ip_list")
    #     ip_ranges = self.rcf.get("Network", "ip_range")
    #     range_list = re.split('\s*,\s*', ip_ranges)
    #     ip_range_lists = []
    #     for rg in range_list:
    #         ip_range_lists.append(IP(rg))
    #     if 'Windows' in platform.system():
    #         local_ip_list = socket.gethostbyname_ex(socket.gethostname())
    #         for lip in local_ip_list[2]:
    #             for sip in ip_range_lists:
    #                 if lip in sip:
    #                     return lip
    #     if platform.system() == 'Linux' or platform.system() == 'Darwin':
    #         # return getip("eth1")
    #         ipstr = '([0-9]{1,3}\.){3}[0-9]{1,3}'
    #         ipconfig_process = subprocess.Popen("/sbin/ifconfig", stdout=subprocess.PIPE)
    #         output = ipconfig_process.stdout.read()
    #         if platform.linux_distribution()[1].startswith("7."):
    #             ip_pattern = re.compile('(inet %s)' % ipstr)
    #         else:
    #         # if platform.system() == "Linux":
    #             ip_pattern = re.compile('(inet addr:%s)' % ipstr)
    #         pattern = re.compile(ipstr)
    #         for ipaddr in re.finditer(ip_pattern, str(output)):
    #             ip = pattern.search(ipaddr.group())
    #             # print ip.group()
    #             for sip in ip_range_lists:
    #                 if ip.group() in sip:
    #                     return ip.group()
    #     raise Exception("内网Network网段设置错误,没有找到IP属于网段%s!" % ip_ranges)
    #
    # @property
    # def LISTEN_PORT(self):
    #     if self._listen_port:
    #         return self._listen_port
    #     else:
    #         return self.get_listen_port()

    # def get_listen_port(self):
    #     # writer = 'yuguo'
    #     """
    #     获取配置文件中start port，end_port ,
    #     在其之间随机生成一个端口，返回一个未被占用的端口
    #     """
    #     start_port = self.rcf.get("Network", "start_port")
    #     end_port = self.rcf.get("Network", "end_port")
    #     ip = self.LISTEN_IP
    #     lpt = random.randint(int(start_port), int(end_port)+1)
    #     # nm = nmap.PortScanner()
    #     # while 1:
    #     #     x = random.randint(int(start_port), int(end_port)+1)
    #     #     # 使用nmap检测端口状态
    #     #     nm.scan(ip, str(x))
    #     #     s = nm[ip]['tcp'][x]['state']
    #     #     if s == 'closed':
    #     #         return x
    #     #         break
    #     # while 1:
    #     #     try:
    #     #         ss = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    #     #         ss.connect((ip, int(lpt)))
    #     #         ss.shutdown(2)
    #     #         # lpt is opened
    #     #         ss.close()
    #     #     except:
    #     #         # lpt is down
    #     #         return lpt
    #     #         break
    #
    #     try:
    #         ss = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    #         ss.connect((ip, int(lpt)))
    #         ss.shutdown(2)
    #         ss.close()
    #     except socket.error:
    #         return lpt
    #     else:
    #         return self.get_listen_port()

    # def get_db(self):
    #     if not self._db_client:
    #         if self.DB_TYPE == "mysql":
    #             # self._db_client = web.database(dbn=self.DB_TYPE, host=self.DB_HOST,
    #             # db=self.DB_NAME, user=self.DB_USER, passwd=self.DB_PASSWD, port=int(self.DB_PORT))
    #             import MySQLdb
    #             self._db_client = MySQLdb.connect(host=self.DB_HOST, user=self.DB_USER, passwd=self.DB_PASSWD,
    #                                               db=self.DB_NAME, port=int(self.DB_PORT))
    #     return self._db_client
    #

    # @staticmethod
    def get_netdata_config(self, type_name):
        print ("work_dir:%s" % self.WORK_DIR)
        # type_list = re.split(r"\s*,\s*", self.rcf.get("NETDATA", "types"))
        # if type_name not in type_list:
        #     raise Exception("Unkown netdata %s" % type_name)
        # options = self.rcf.options("NETDATA")
        # type_dict = {}
        # for opt in options:
        #     if re.match("^" + type_name, opt):
        #         type_dict[opt] = self.rcf.get("NETDATA", opt)
        # return type_dict
        libs = {
            "sanger_type": "httpcache",
            "sanger_path": "http://bcl.i-sanger.com/data/",

            "tsanger_type": "httpcache",
            "tsanger_path": "http://bcl.tsanger.com/data/",

            "tsg_type": "httpcache",
            "tsg_path": "http://bcl.tsanger.com/data/",

            "i-sanger_type": "httpcache",
            "i-sanger_path": "http://bcl.i-sanger.com/data/",

            "s3_type": "s3cache",
            "s3_cache_dir": os.path.join(self.WORK_DIR, "s3cache"),

            "http_type": "httpcache",
            "http_cache_dir": os.path.join(self.WORK_DIR, "httpcache"),
            # "http_cache_dir": "/mnt/lustre/users/sanger/workspace/httpcache",
        }
        type_dict = {}
        for opt in libs.keys():
            if re.match("^" + type_name, opt):
                type_dict[opt] = libs[opt]
        return type_dict

    @staticmethod
    def get_netdata_lib(type_name):
        libs = {
            "sanger": "httpcache",
            "tsanger": "httpcache",
            "tsg": "httpcache",
            "i-sanger": "httpcache",
            "s3": "s3cache",
            "http": "httpcache",
            "filelist": "filelistcache",
            "fileapi": "fileapicache"
        }
        if type_name in libs.keys():
            return libs[type_name]

    # def get_api_type(self, client):
    #     if self.rcf.has_option("API", client):
    #         return self.rcf.get("API", client)
    #     else:
    #         return None
    def get_api_type(self, client):
        api_dict = {"client01": "sanger", "client02": "split_data", "client03": "tsanger"}
        if client in api_dict.keys():
            return api_dict[client]
        else:
            return None

    # def get_use_api_clients(self):
    #     return self.rcf.options("API")
    #
    # @property
    # def UPDATE_EXCLUDE_API(self):
    #     if self._update_exclude_api is not None:
    #         return self._update_exclude_api
    #     if self.rcf.has_option("API_UPDATE", "exclude_api"):
    #         exclude_api = self.rcf.get("API_UPDATE", "exclude_api")
    #         if exclude_api.strip() != "":
    #             self._update_exclude_api = re.split(r"\s*,\s*", exclude_api)
    #             return self._update_exclude_api
    #         else:
    #             self._update_exclude_api = []
    #             return self._update_exclude_api
    #     else:
    #         self._update_exclude_api = []
    #         return self._update_exclude_api
    #
    # def __del__(self):
    #     if self._mongo_client:
    #         self.mongo_client.close()
    #     if self._biodb_mongo_client:
    #         self._biodb_mongo_client.close()
    #

    def get_file_type(self, file_path):
        type_name = "local"
        m = re.match(r"^http://.*|^https://.*", file_path)
        if m:
            type_name = "http"
        else:
            m1 = re.match(r"^([\w\-]+)://.*", file_path)
            if m1:
                type_name = "s3"
            else:
                m2 = re.match(r"^([\w\-]+):/*(.*)$", file_path)
                if m2:
                    type_name = self.get_netdata_lib(m2.group(1))
                    if type_name == "httpcache":
                        type_name = "http"
                else:
                    for k, v in self.get_convert_list().items():
                        if file_path.startswith(k):
                            type_name = "http"
        return type_name

    @staticmethod
    def get_convert_list():
        # convert_list = re.split(r"\s*\|+\s*", self.rcf.get("HTTP", "path_to_convert"))
        convert_list = {
            "/mnt/ilustre/data/": "http://bcl.i-sanger.com/data/",
            "/mnt/ilustre/tsanger-data/": "http://bcl.tsanger.com/data/",
        }
        # for c in convert_list:
        #     m = re.match(r"^(.*)\((.*)\)$", c)
        #     if m:
        #         prefix = m.group(1)
        #         to_prefix = m.group(2)
        #         self._http_convert_list[prefix] = to_prefix
        return convert_list

    def convert_path_to_http(self, path):
        m = re.match(r"^http://.*|^https://.*", path)
        if m:
            return path
        else:
            m1 = re.match(r"^([\w\-]+):/?(.*)$", path)
            if m1:
                head = m1.group(1)
                key = m1.group(2)
                return os.path.join(self.get_netdata_config(head)["%s_path" % head], key)
            else:
                for k, v in self.get_convert_list().items():
                    if path.startswith(k):
                        return path.replace(k, v)
                return path

    def convert_real_path(self, path):
        type_name = self.get_file_type(path)
        if type_name == "s3" or type_name == "local":
            return path
        elif type_name == "http":
            return self.convert_path_to_http(path)
        else:
            m = re.match(r"^([\w\-]+):/*(.*)$", self._remote_path)
            if m:
                file_type = m.group(1)
                file_path = m.group(2)
                return os.path.join(self.get_netdata_config(file_type)["%s_path" % file_type], file_path)
        return path

    def get_work_dir(self):
        """
        该函数是rna那边要用的
        这部分是在接口中调用，因为接口中调用这一步的时候，还没有初始化init——config，所以这单独处理
        :param :
        :return:
        """
        self.init_config()
        return self._workdir
