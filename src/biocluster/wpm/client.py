# -*- coding: utf-8 -*-
# __author__ = 'HD'
import grpc
import json
import gevent
from biocluster.proto import webapi_pb2, webapi_pb2_grpc
from biocluster.config import Config
# from .db import Postgresql
import sys


def worker_client():
    return Run()


class Run(object):
    """
    该函数功能是提交workflow的任务到workflow队列，用于发起工作流+发送工作流的暂停，取消，取消暂停的消息。
    提交工作流的时候json_data与当前流程中的json.data是一致的
    提交工作流的运行状态的时候，json_data:{"id": workflow_id, "msg": "stop"}
    """
    def __init__(self):
        self.json_data = {}
        self.instant = False
        self.client = "client03"
        self.config = Config()
        self.timeout = 600
        # self._grpc_times = 0
        self._cmdreq_times = 0
        self._filedeletereq_times = 0

    def add_task(self, jsondata):
        self.json_data = jsondata
        if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
            print "add workflow %s format error!" % self.json_data
            return {"success": False, "info": "json格式错误"}
        wfm_port = self.config.wfm_port
        if "instant" in self.json_data.keys():
            self.instant = self.json_data["instant"]
        if "wfm_port" in self.json_data.keys():
            wfm_port = self.json_data['wfm_port']
        print "wfm_port:{}".format(wfm_port)
        if "client" in self.json_data.keys():
            self.client = self.json_data['client']
        # if "DBVersion" not in self.json_data.keys():  # 默认DBVersion是1
        #     self.json_data["DBVersion"] = 1
        self.timeout = self.config.WFM_INSTANT_TIMEOUT if self.instant else self.config.WFM_SUBMIT_TIMEOUT
        self.json_data["WPM"] = True
        try:
            with grpc.insecure_channel('localhost:%s' % wfm_port) as channel:
                stub = webapi_pb2_grpc.WebApiStub(channel)
                response = stub.Submit(webapi_pb2.Task(
                    id=self.json_data["id"],
                    client=self.client,
                    json=json.dumps(self.json_data),
                    instant=self.instant
                ), timeout=self.timeout)
        except Exception as e:
            print "GRPC请求出错:{}".format(e)
            sys.stdout.flush()
            return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            # if self._grpc_times > 2:
            #     print "GRPC请求出错:{}".format(e)
            #     return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            # else:
            #     print("GRPC请求错误10秒后重试: %s" % e)
            #     gevent.sleep(5)
            #     self.run_task()
        else:
            return {"success": response.ok, "info": response.reason}

    def run_cmd(self, jsondata):
        """
        发送取消，暂停，取消暂停等指令
        :return:
        """
        self._cmdreq_times += 1
        self.json_data = jsondata
        if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
            print "json格式不正确-样例：{'id': tsg_0001, 'msg': 'stop'}" % self.json_data
            return {"success": False, "info": "json格式错误"}
        wfm_port = self.config.wfm_port
        if "wfm_port" in self.json_data.keys():
            wfm_port = self.json_data['wfm_port']
        try:
            with grpc.insecure_channel('localhost:%s' % wfm_port) as channel:
                stub = webapi_pb2_grpc.WebApiStub(channel)
                response = stub.Control(webapi_pb2.Command(
                    id=self.json_data["id"],
                    msg=self.json_data["msg"],
                ), timeout=self.timeout)
                print response
        except Exception as e:
            # print "GRPC请求出错:{}".format(e)
            # return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            if self._cmdreq_times > 2:
                print "GRPC请求出错:{}".format(e)
                sys.stdout.flush()
                return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            else:
                print("GRPC请求错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(5)
                self.run_cmd(jsondata)
        else:
            return {"success": response.ok, "info": response.reason}

    def add_filedelete_task(self, jsondata):
        """
        添加删除对象存储上文件的任务。
        :param jsondata:
        :return:
        """
        self._filedeletereq_times += 1
        self.json_data = jsondata
        if not isinstance(self.json_data, dict) or "project_id" not in self.json_data.keys() or\
                "task_id" not in self.json_data.keys() or "unique_id" not in self.json_data.keys():
            print "json格式不正确-样例：{'project_id': 'ungnprbq2h8igqigalcufde9i8', " \
                  "'task_id': 'tsg_002', 'unique_id': 'a6173961e9d6143877decf2b6a4f37ca'}" % self.json_data
            return {"success": False, "info": "json格式错误"}
        wfm_port = self.config.wfm_port
        if "wfm_port" in self.json_data.keys():
            wfm_port = self.json_data['wfm_port']
        try:
            with grpc.insecure_channel('localhost:%s' % wfm_port) as channel:
                stub = webapi_pb2_grpc.WebApiStub(channel)
                response = stub.S3delete(webapi_pb2.Proinfo(
                    projectid=self.json_data["project_id"],
                    taskid=self.json_data["task_id"],
                    uniqueid=self.json_data["unique_id"]
                ), timeout=120)
                # print response
        except Exception as e:
            if self._filedeletereq_times > 2:
                print "GRPC请求出错:{}".format(e)
                sys.stdout.flush()
                return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            else:
                print("GRPC请求错误30秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(30)
                self.add_filedelete_task(jsondata)
        else:
            return {"success": response.ok, "info": response.reason}


def wait():
    pass


# class AddTask(object):
#     """
#     该函数功能是提交workflow的任务到workflow队列，用于发起工作流+发送工作流的暂停，取消，取消暂停的消息。
#     提交工作流的时候json_data与当前流程中的json.data是一致的
#     提交工作流的运行状态的时候，json_data:{"id": workflow_id, "msg": "stop"}
#     """
#     def __init__(self, jsondata, timeout=None):
#         self.json_data = jsondata
#         self.instant = False
#         self.client = "client03"
#         self.config = Config()
#         self.timeout = 600
#         # self._grpc_times = 0
#         self._cmdreq_times = 0
#         if timeout:
#             self.timeout = timeout
#
#     def run_task(self):
#         # self._grpc_times += 1
#         if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
#             print "add workflow %s format error!" % self.json_data
#             return {"success": False, "info": "json格式错误"}
#         self.json_data["WPM"] = True
#         wfm_port = self.config.wfm_port
#         if "wfm_port" in self.json_data.keys():
#             wfm_port = self.json_data['wfm_port']
#         print "wfm_port:{}".format(wfm_port)
#         try:
#             with grpc.insecure_channel('localhost:%s' % wfm_port) as channel:
#                 stub = webapi_pb2_grpc.WebApiStub(channel)
#                 response = stub.Submit(webapi_pb2.Task(
#                     id=self.json_data["id"],
#                     client=self.client,
#                     json=json.dumps(self.json_data),
#                     instant=self.instant
#                 ), timeout=self.timeout)
#                 # print "----------------------------"
#                 # print response
#         except Exception as e:
#             print "GRPC请求出错:{}".format(e)
#             return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
#             # if self._grpc_times > 2:
#             #     print "GRPC请求出错:{}".format(e)
#             #     return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
#             # else:
#             #     print("GRPC请求错误10秒后重试: %s" % e)
#             #     gevent.sleep(5)
#             #     self.run_task()
#         else:
#             return {"success": response.ok, "info": response.reason}
#
#     def run_cmd(self):
#         """
#         发送取消，暂停，取消暂停等指令
#         :return:
#         """
#         print "000000000000000000000000000000"
#         self._cmdreq_times += 1
#         if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
#             print "json格式不正确-样例：{'id': tsg_0001, 'msg': 'stop'}" % self.json_data
#             return {"success": False, "info": "json格式错误"}
#         try:
#             with grpc.insecure_channel('localhost:%s' % self.config.wfm_port) as channel:
#                 stub = webapi_pb2_grpc.WebApiStub(channel)
#                 response = stub.Control(webapi_pb2.Command(
#                     id=self.json_data["id"],
#                     msg=self.json_data["msg"],
#                 ), timeout=self.timeout)
#                 print response
#         except Exception as e:
#             # print "GRPC请求出错:{}".format(e)
#             # return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
#             if self._cmdreq_times > 2:
#                 print "GRPC请求出错:{}".format(e)
#                 return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
#             else:
#                 print("GRPC请求错误10秒后重试: %s" % e)
#                 gevent.sleep(5)
#                 self.run_cmd()
#         else:
#             return {"success": response.ok, "info": response.reason}


if __name__ == "__main__":
    # a = AddTask(jsondata={"id": "tsg_10017", "msg": "stop"})
    # a = AddTask(jsondata={"id": "tsg_10019", "msg": "pause"})
    worker = worker_client()
    response = worker.run_cmd(jsondata={"id": "tsg_10019", "msg": "continue"})
    print response
