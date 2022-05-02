# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import zerorpc
from ..core.function import hostname
from ..core.singleton import singleton
import time
import gevent
import traceback
import sys
import copy
import os
from multiprocessing import Process
import setproctitle
from .logger import Logger
from ..config import Config
import json
from gevent.lock import BoundedSemaphore


class RpcServerProcess(Process):
    def __init__(self, endpoint, message_queue, callback_queue, restart_signal):
        super(RpcServerProcess, self).__init__()
        self.manager = Manager(endpoint, message_queue, callback_queue, restart_signal)
        # self.logger = Logger("")

    def run(self):
        super(RpcServerProcess, self).run()
        setproctitle.setproctitle("WPM[RPC Server]")
        self.manager.run()


@singleton
class Manager(object):
    """
    RPCServer
    """
    def __init__(self, endpoint, message_queue, callback_queue, restart_signal):
        self.endpoint = endpoint
        self.message_queue = message_queue
        self.callback_queue = callback_queue
        self.restart_signal = restart_signal
        self._messages = []
        self._callbacks = {}
        self.last_get_message = None
        self.logger = Logger("RPC Server")
        self.rpc_server = RPC(self)
        self._log_date = None
        self.config = Config()
        self.sem = BoundedSemaphore(1)

    def report(self, msg):
        self.last_get_message = int(time.time())
        self.logger.info("接收到RPC数据: %s " % msg)
        if (not isinstance(msg, dict)) or ('id' not in msg.keys()) or ("version" not in msg.keys()):
            self.logger.error("Server接收到不符合规范的消息: 不是字典类型或没有key值'id'或者'version'!")
            return
        if len(msg["state"]) > 50:
            self.logger.error("接收到信息: %s, state 长度超过最长50个限制，丢弃！" % msg)
            return
        if len(msg["id"]) > 500:
            self.logger.error("接收到信息: %s, id 长度超过最长500个限制，丢弃!" % msg)
            return
        msg["data"] = json.dumps(msg["data"])
        if len(msg["data"]) > 2048:
            self.logger.warning("接收到信息: %s ,data长度超过最长2048个字节限!" % msg)
            return
        msg["version"] = int(msg["version"])
        with self.sem:
            self._messages.append(msg)
        key = "%s|%s" % (msg["id"], msg["version"])
        if key in self._callbacks.keys():
            return self._callbacks[key]
        else:
            return {'action': 'none'}

    def _set_message_queue(self):
        while True:
            if len(self._messages) > 0:
                try:
                    with self.sem:
                        data_list = copy.copy(self._messages)
                        self._messages = []
                    with self.message_queue.get_lock():
                        for i, v in enumerate(self.message_queue):
                            if v.id == "":
                                if len(data_list) > 0:
                                    msg = data_list.pop(0)
                                    self.message_queue[i] = (msg["id"], msg["state"], msg["jobid"], msg["host"],
                                                             msg["data"], msg["version"], int(time.time()))
                                    # self.logger.debug("保存数据1: %s " % self.message_queue[i])
                                else:
                                    break
                            elif v.timestamp > 0 and (int(time.time()) - v.timestamp) > 600:
                                if len(data_list) > 0:
                                    msg = data_list.pop(0)
                                    self.logger.error("信息%s超过10分钟没有被获取，丢弃!" % v)
                                    self.message_queue[i] = (msg["id"], msg["state"], msg["jobid"], msg["host"],
                                                             msg["data"], msg["version"], int(time.time()))
                                    # self.logger.debug("保存数据2: %s " % self.message_queue[i])
                                else:
                                    break
                            else:
                                if len(data_list) > 0:
                                    found = None
                                    for m in data_list:
                                        if m["id"] == v.id and v.state == "keepalive":
                                            found = m
                                            self.message_queue[i] = (m["id"], m["state"], msg["jobid"], msg["host"],
                                                                     m["data"], m["version"], int(time.time()))
                                            # self.logger.debug("保存数据3: %s " % self.message_queue[i])
                                            break
                                    if found:
                                        data_list.remove(found)
                                else:
                                    break
                    if len(data_list) > 0:
                        self.logger.warning("队列已满，稍后再试...")
                        with self.sem:
                            self._messages.extend(data_list)
                        gevent.sleep(3)
                except Exception, e:
                    exstr = traceback.format_exc()
                    print exstr, e
                    sys.stdout.flush()
                    self.logger.error("设置RPC队列信息错误...")
            gevent.sleep(1)

    def _get_callback_queue(self):
        while True:
            try:
                msg_list = []
                with self.callback_queue.get_lock():
                    for i, v in enumerate(self.callback_queue):
                        if v.id != "":
                            msg_list.append(v.to_dict())
                            self.callback_queue[i] = ("", "", "", 0)
                timestamp = int(time.time())
                for msg in msg_list:
                    msg["timestamp"] = timestamp
                    key = "%s|%s" % (msg.id, msg.version)
                    self._callbacks[key] = msg
                for key in self._callbacks.keys():
                    if int(time.time()) - self._callbacks[key]["timestamp"] > 600:
                        self._callbacks.pop(key)
                        self.logger.error("CallBack信息%s超过10分钟没有被获取，丢弃!" % self._callbacks[key])
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
                self.logger.error("获取Callback队列信息错误...")
            gevent.sleep(1)

    def _check_rpc(self):
        while True:
            try:
                if self.restart_signal.value == 1:
                    self.restart_signal.value = 0
                    if self.last_get_message and (int(time.time()) - self.last_get_message) < 100:
                        self.logger.info("RPC服务%s秒以前还接收到RPC信息，还在正常工作，忽略重启信号!" %
                                         (int(time.time()) - self.last_get_message))
                    else:
                        self.rpc_server.restart()
                if os.getppid() == 1:
                    self.rpc_server.close()
                    self.logger.info("检测到主服务已经退出，退出RPC服务!")
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
                self.logger.error("检测RPC服务出错!...")
            gevent.sleep(1)

    def _check_date(self):
        while True:
            try:
                if self._log_date != time.strftime('%Y%m%d', time.localtime(time.time())):
                    self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
                    log = os.path.join(self.config.wpm_log_file, "%s.%s.rpc.log" % (self._log_date, hostname))
                    so = file(log, 'a+')
                    se = file(log, 'a+', 0)
                    os.dup2(so.fileno(), sys.stdout.fileno())
                    os.dup2(se.fileno(), sys.stderr.fileno())
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
            gevent.sleep(60)

    def run(self):
        self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
        log = os.path.join(self.config.wpm_log_file, "%s.%s.rpc.log" % (self._log_date, hostname))
        so = file(log, 'a+')
        se = file(log, 'a+', 0)
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())
        self.logger.info("启动RPC服务!")
        gevent.spawn(self._set_message_queue)
        gevent.spawn(self._get_callback_queue)
        gevent.spawn(self._check_rpc)
        gevent.spawn(self._check_date)
        self.rpc_server.run()


class RPC(object):
    def __init__(self, manager):
        self.manager = manager
        self._rpc_server = None
        self._closed = False
        self._last_start = None
        self._retry_times = 0

    def run(self):
        """
        开始运行RPC监听,此时会阻塞线程
        """
        gevent.spawn(self._start_listen)
        while not self._closed:
            gevent.sleep(1)

    def _start_listen(self):
        try:
            self._rpc_server = zerorpc.Server(self.manager)
            self._rpc_server.bind(self.manager.endpoint)
            self._rpc_server.run()
            self._last_start = int(time.time())
        except Exception, e:
            exstr = traceback.format_exc()
            print e
            print exstr
            sys.stdout.flush()
            self._retry_times += 1
            if self._retry_times < 3:
                gevent.sleep(3)
                self.manager.logger.debug("启动监听错误,尝试重试1")
                self._start_listen()
            else:
                self.manager.logger.error("启动监听错误超过3次，终止运行!")
                return

    def close(self):
        self.manager.logger.debug("关闭RPC监听")
        self._closed = True
        self._rpc_server.close()
        self.manager.logger.debug("关闭RPC监听成功")

    def restart(self):
        if self._last_start and (int(time.time()) - self._last_start) < 300:
            self.manager.logger.debug("RPC服务器启动不超过5分钟，忽略重启!")
            return
        self._retry_times = 0
        self.manager.logger.warning("重启监听服务%s" % self.manager.endpoint)
        self._rpc_server.close()
        gevent.sleep(3)
        gevent.spawn(self._start_listen)
