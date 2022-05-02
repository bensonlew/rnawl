# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""rpc远程客户端和服务端"""

import zerorpc
from .config import Config
import datetime
import gevent
from multiprocessing import Queue
import traceback
import sys
from .core.function import MaxLengthError
import json


class Report(object):
    """
    RPCServer
    """
    def __init__(self, workflow):
        self.workflow = workflow

    def report(self, msg, is_rpc=True):
        """
        用于传递消息的方法。

        :param msg: 通过RPC服务接收的远程消息
        :param is_rpc: 是否通过RPC调用
        """
        if is_rpc is True:
            try:
                msg["data"] = json.loads(msg["data"])
            except Exception:
                pass
        self.workflow.logger.debug("接收到RPC数据: %s " % msg)
        self.workflow.last_get_message = datetime.datetime.now()
        agent = self.workflow.find_tool_by_id(msg['id'])
        if not agent:
            self.workflow.logger.error("Server在workflow中找不到对应的tool: {} !".format(msg['id']))
        else:
            # agent.send_exit_action("Server在workflow中找不到对应的tool: {} !".format(msg['id']), msg["version"])
            if (not isinstance(msg, dict)) or ('id' not in msg.keys()) or ("version" not in msg.keys()):
                self.workflow.logger.error("Server接收到不符合规范的消息: 不是字典类型或没有key值'id'或者'version'!")
            elif msg["version"] != agent.version:
                self.workflow.logger.error("接收到已经重新投递任务的历史版本信号，丢弃: %s ！" % msg)
                agent.send_exit_action("此版本号已经更新过期", msg["version"])
            else:
                self.workflow.last_update = datetime.datetime.now()
                if agent.job:
                    agent.actor.receive(msg)
                else:
                    agent.send_exit_action(reason="Job尚未开始却接收到信息", version=msg['version'])
            if not self.workflow.sheet.instant:
                return agent.get_callback_action(msg['version'])


class RPC(object):
    def __init__(self, workflow):
        self._workflow = workflow
        self._rpc_server = zerorpc.Server()
        self._report = Report(self._workflow)
        config = Config()
        if workflow.sheet.endpoint:
            self.endpoint = workflow.sheet.endpoint
        else:
            self.endpoint = "tcp://{}:{}".format(config.LISTEN_IP, config.LISTEN_PORT)
        # self._rpc_server.bind(self.endpoint)
        self._closed = False
        self._last_start = None
        self._retry_times = 0
        self.last_none_block = None

    def run(self):
        """
        开始运行RPC监听,此时会阻塞线程
        """
        if self._workflow.sheet.WPM:
            gevent.spawn(self._get_rpc_info)
        else:
            gevent.spawn(self._start_listen)
        while not self._closed:
            gevent.sleep(1)

    def _get_rpc_info(self):
        message_queue = self._workflow.rpc_message_queue
        self._workflow.logger.info("开始与WPM RPC Server通信!")
        while True:
            if self._closed:
                break
            if self.last_none_block:
                block_time = (datetime.datetime.now() - self.last_none_block).seconds
                if block_time >= 600:
                    self._workflow.logger.warning("进程被阻塞超过%s秒，可能已经导致RPC数据丢失!请减少阻塞时间！"
                                                  "多个阻塞循环可在循环中使用gevent.sleep(0)跳出阻塞!" % block_time)
                elif block_time > 30:
                    self._workflow.logger.warning("进程被阻塞超过%s秒，长时间阻塞会导致RPC数据丢失,请减少阻塞时间! "
                                                  "多个阻塞循环可在循环中使用gevent.sleep(0)跳出阻塞!" % block_time)
            try:
                self.last_none_block = datetime.datetime.now()
                message_list = []
                with message_queue.get_lock():
                    for i, v in enumerate(message_queue):
                        if v.id != "" and v.id.split(".")[0] == self._workflow.id:
                            message_list.append(v.to_dict())
                            # self._workflow.logger.debug("检测到%s!" % v)
                            message_queue[i] = ("", "", "", "", "", 0, 0)
                for msg in message_list:
                    callback_msg = self._report.report(msg)
                    if callback_msg and callback_msg["action"] != "none":
                        self._set_callback_queue(msg)
            except Exception, e:
                exstr = traceback.format_exc()
                print e
                print exstr
                sys.stdout.flush()
                self._workflow.logger.error("获取设置RPC信息出错!")
            gevent.sleep(3)

    def _set_callback_queue(self, data):
        if len(data["id"]) > 500:
            raise MaxLengthError("数据%s id 长度超过最长500个限制!" % data)
        if len(data["action"]) > 50:
            raise MaxLengthError("数据%s action 长度超过最长50个限制!" % data)
        data["data"] = json.dumps(data["data"])
        if len(data["data"]) > 2048:
            self._workflow.logger.error("数据: %s ,data长度超过最长2048个限制, 截取数据片段!" % data)
            data["data"] = data["data"].decode("utf-8")[0:2000].encode("utf-8")
        callback_queue = self._workflow.rpc_callback_queue
        with callback_queue.get_lock():
            for i, v in enumerate(callback_queue):
                if v.id == "":
                    callback_queue[i] = (data["id"], data["action"], data["data"], data["version"])
                    break

    def _start_listen(self):
        try:
            self._rpc_server = zerorpc.Server(self._report)
            self._rpc_server.bind(self.endpoint)
            self._rpc_server.run()
            self._last_start = datetime.datetime.now()
        except Exception, e:
            exstr = traceback.format_exc()
            print e
            print exstr
            sys.stdout.flush()
            self._retry_times += 1
            if self._retry_times < 3:
                gevent.sleep(3)
                self._workflow.logger.debug("启动监听错误,尝试重试！")
                self._start_listen()
            else:
                self._workflow.exit(data="启动监听错误超过3次，终止运行!")
                return

    def close(self):
        if self._workflow.sheet.WPM:
            self._workflow.logger.debug("关闭RPC通信")
            self._closed = True
        else:
            self._workflow.logger.debug("关闭RPC监听")
            self._closed = True
            self._rpc_server.close()
            self._workflow.logger.debug("关闭RPC监听成功")

    def restart(self):
        if self._last_start and (datetime.datetime.now() - self._last_start).seconds < 300:
            self._workflow.logger.debug("RPC服务器启动不超过5分钟，忽略重启!")
            return
        if self._workflow.sheet.WPM:
            self._workflow.rpc_restart_signal.value = 1
            self._last_start = datetime.datetime.now()
        else:
            self._retry_times = 0
            self._workflow.logger.warning("重启监听服务%s" % self.endpoint)
            self._rpc_server.close()
            gevent.sleep(3)
            gevent.spawn(self._start_listen)


class LocalServer(object):
    def __init__(self, workflow):
        self._report = Report(workflow)
        self._workflow = workflow
        self._close = False
        # (reader, writer) = gipc.pipe()
        # print "reader %s, writer %s" % (reader, writer)

        self._process_queue = None
        self.endpoint = ""

    @property
    def process_queue(self):
        if not self._process_queue:
            self._process_queue = Queue()
        return self._process_queue

    def run(self):
        self._workflow.logger.info("开始启动与本地WPM RPC服务通信!")
        while True:
            gevent.sleep(0.3)
            try:
                msg = self.process_queue.get_nowait()
            except Exception:
                pass
            else:
                if msg:
                    self._report.report(msg, False)
            if self._close:
                break

    def close(self):
        self._close = True


#
# class RPCClient(zerorpc.Client):
#     """
#     sent message of each tool running states
#     发送tool的状态
#     """
#     def __init__(self, endpoint, tool):
#         """
#         link to a server.
#         连接到一个服务端
#         """
#         super(RPCClient, self).__init__()
#         self.connect(endpoint)
#         self.tool = tool
#
#     def run(self, msg):
#         """
#         sent message to server.
#         发送消息给服务端, 并获得返回的消息，将获得的消息告诉Tool
#         """
#         try:
#             feedback = self.rpc_letter(msg)
#         except Exception, e:
#             self._logger.error("Client 运行出现错误: {}".format(e))
#         else:
#             if feedback is None:
#                 pass
#                 self._logger.info("Client 没有反馈消息给tool")
#             elif feedback[0] == "error":
#                 self.tool.kill()
#                 self._logger.info("Client 接到终止tool运行消息，终止{}".format(self.tool))
#             else:
#                 getattr(self.tool, feedback[1])
#                 self._logger.info("Client 传递消息给{}执行：{}".format(self.tool, feedback[1]))
