# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""actor消息处理机制"""

import gevent
from gevent import Greenlet
import datetime
import threading
import inspect
import zerorpc
import platform
import traceback
import sys
import time
import os
import pickle
import copy
import json


class State(object):
    def __init__(self, name, data=None):
        self.name = name
        self.data = data


class LocalActor(gevent.Greenlet):
    """
    继承至Greenlet,每个Agent对象都会生成一个LocalActor对象，循环处理接受的消息并负责记录最新状态更新时间,根据接收到的不同信息调用不同的处理函数,如果处理函数不存在则调用默认处理函数予以提示
    """
    def __init__(self, agent):
        # self.inbox = Queue()
        self._config = agent.config
        self._agent = agent
        # self._name = agent.name + ".Actor"
        self._update = None    # 最近接受消息时间
        self._start_time = None  # 开始运行时间
        self._kao_fire_times = 0
        self._has_firekaoout = False
        self._wto_fire_times = 0
        self._has_firewtoout = False
        self.auto_break = False  # 自动退出
        self._run_startout_times = 0
        self._recive_finish_state = False
        Greenlet.__init__(self)

    def check_time(self):
        now = datetime.datetime.now()
        workflow = self._agent.get_workflow()
        if workflow.sheet.WPM and workflow.rpc_server.last_none_block:
            if (now - workflow.rpc_server.last_none_block).seconds > 10:
                return
        if self._agent.is_end or self.auto_break:
            return
        if self._update is not None:
            if (now - self._update).seconds - self._kao_fire_times * self._config.MAX_KEEP_ALIVE_TIME\
                    > self._config.MAX_KEEP_ALIVE_TIME:
                if self._kao_fire_times >= self._config.MAX_FIRE_KAO_TIMES:
                    if self._has_firekaoout is False:
                        self._agent.fire("firekaoout", self._kao_fire_times)
                        self._has_firekaoout = True
                    return
                self._agent.logger.warning("KeepAlive通信超时%s秒,第%s次触发!" %
                                           ((self._kao_fire_times + 1) * self._config.MAX_KEEP_ALIVE_TIME,
                                            self._kao_fire_times + 1))
                self._agent.fire('keepaliveout')
                self._kao_fire_times += 1

        if self._update is None and (not self._agent.is_wait):
                if self._agent.job.submit_time is not None:
                    if (now - self._agent.job.submit_time).seconds \
                            - self._run_startout_times * self._config.RUN_START_TIMEOUT \
                            > self._config.RUN_START_TIMEOUT:
                        self._run_startout_times += 1

                        self._agent.fire('runstartout', self._run_startout_times)

                    if (now - self._agent.job.submit_time).seconds - \
                            self._wto_fire_times * self._config.MAX_WAIT_TIME > self._config.MAX_WAIT_TIME:
                        if self._wto_fire_times >= self._config.MAX_FIRE_WTO_TIMES:
                            if self._has_firewtoout is False:
                                self._agent.fire("firewtoout", self._wto_fire_times)
                                self._has_firewtoout = True
                            return
                        self._agent.logger.warning("等待运行超时%s秒,第%s次触发!" %
                                                   ((self._wto_fire_times + 1) * self._config.MAX_WAIT_TIME,
                                                    (self._wto_fire_times + 1)))
                        self._agent.fire('waittimeout')
                        self._wto_fire_times += 1

    def receive(self, message):
        """
        接收消息后解析消息数据，动态调用对应的处理方法

        :param message: message为远程rpc传递的数据,python dict类型数据 必须包含 key "event"
        """
        if self._agent.is_start is False:
            self._agent.logger.debug("尚未开始运行，丢弃接收到的消息: %s " % message)
            return
        if self._agent.is_end:
            self._agent.logger.debug("已经停止运行，丢弃接收到的消息: %s " % message)
            return
        workflow = self._agent.get_workflow()
        try:
            if not workflow.sheet.instant and self._update is None:
                if "host" in message.keys():
                    self._agent.fire('runstart', message["host"])
                else:
                    self._agent.fire('runstart', "unknown")
            self._update = datetime.datetime.now()
            if (not isinstance(message, dict)) or ('state' not in message.keys()):
                self._agent.logger.warning("接收到不符合规范的消息，丢弃!")
                return
            if "jobid" in message.keys():
                self._agent.job.id = message["jobid"]
            if message['state'] == "finish":
                if self._recive_finish_state:
                    self._agent.logger.debug("已经接收到finish State后再次接收到finish,丢弃!")
                    return
                self._recive_finish_state = True
            if message['state'] != "keepalive":
                self._agent.fire("recivestate", message)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            workflow.set_error(e)

    def _run(self):
        self._start_time = datetime.datetime.now()
        self.running = True
        while self.running:
            # message = self.inbox.get()
            # self.receive(message)
            if self._agent.is_end or self.auto_break:
                break
            self.check_time()
            gevent.sleep(3)

    def close(self):
        self.auto_break = True


class RemoteActor(threading.Thread):
    """
    每个Tool远程运行时都会产生一个RemoteActor对象,负责远端运行时的消息处理机制。
    负责将Tool运行过程中添加的State状态发送到其对应的Agent对象，并调用对应的函数进行处理。如果Agent返回了Action命令，则调用对于的Action方法
    """
    def __init__(self, tool, main_thread):
        super(RemoteActor, self).__init__()
        self._tool = tool
        self.config = tool.config
        self.mutex = tool.mutex
        self._lost_connection_count = 0
        self.main_thread = main_thread

    def run(self):
        """
        定时检测Tool的State状态，并发送状态到远端

        :return: None
        """
        super(RemoteActor, self).run()
        self._tool.logger.debug("Actor启动，开始工作!")

        def is_finished():
            if len(self._tool.states) > 0:
                return False
            elif self._tool.is_end:
                self._tool.logger.debug("Tool结束，Actor退出运行!")
                return True
            else:
                return False
        find_parent_exit = None
        while not is_finished():
            try:
                if os.getppid() == 1:
                    if find_parent_exit is None:
                        self._tool.logger.debug("检测到父进程PID为1,等待检测是否接收到Linux退出信号!")
                        find_parent_exit = int(time.time())
                        if self._tool.receive_exit_signal is True:
                            self._tool.logger.debug("接收到Linux退出信号!")
                        else:
                            self._tool.logger.debug("没有接收到Linux退出信号!")
                            self._tool.exit_handler(0, 0)
                    elif (int(time.time()) - find_parent_exit) > 10 and len(self._tool.states) == 0:
                        self._tool.logger.debug("检测到父进程PID为1,退出Actor运行!")
                        self._tool.exit(1)
                        break
                if self._tool.exit_signal and len(self._tool.states) == 0:
                    self._tool.logger.debug("接收到退出信号，终止Actor信号发送!")
                    self._tool.exit(1)
                    break
                if not self.main_thread.is_alive() and len(self._tool.states) == 0 \
                        and self._tool.exit_signal is not True:
                    self.send_state(State('error', "检测到远程主线程异常结束"))
                    self._tool.logger.debug("检测到主线程已退出，终止运行!")
                    self._tool.exit(1)
                    break
                if len(self._tool.states) > 0:
                    with self.mutex:
                        state = self._tool.states.pop(0)
                    action = self.send_state(state)
                    if isinstance(action, dict) and 'action' in action.keys():
                        if action['action'] != "none":
                            if hasattr(self._tool, action['action'] + '_action'):
                                func = getattr(self._tool, action['action'] + '_action')
                                argspec = inspect.getargspec(func)
                                args = argspec.args
                                if len(args) == 1:
                                    func()
                                elif len(args) == 2:
                                    func(action['data'])
                                else:
                                    raise Exception("action处理函数参数不能超过2个(包括self)!")
                            else:
                                self._tool.logger.warn("没有为返回action %s设置处理函数!" % action['action'])
                    if state.name in {"finish", "error", "memory_limit", "cancelled"}:
                        self._tool.logger.debug("state name: %s " % state.name)
                        self._tool.exit(0)
                        break
                else:
                    if not self._tool.instant:
                        action = self.send_state(State('keepalive'))
                        if isinstance(action, dict) and 'action' in action.keys():
                            if action['action'] != "none":
                                if hasattr(self._tool, action['action'] + '_action'):
                                    func = getattr(self._tool, action['action'] + '_action')
                                    argspec = inspect.getargspec(func)
                                    args = argspec.args
                                    if len(args) == 1:
                                        func()
                                    elif len(args) == 2:
                                        func(action['data'])
                                    else:
                                        raise Exception("action处理函数参数不能超过2个(包括self)!")
                                else:
                                    self._tool.logger.warn("没有为返回action %s设置处理函数!" % action['action'])
            except Exception, e:
                exstr = traceback.format_exc()
                print >> sys.stderr, exstr
                print >>sys.stderr, e
                sys.stderr.flush()
            time.sleep(int(self.config.KEEP_ALIVE_TIME))

    def send_state(self, state):
        """
        发送State信息到远程Actor接收

        :param state: State对象
        :return:
        """
        jobid = None
        if "SLURM_JOB_ID" in os.environ.keys():
            jobid = os.environ["SLURM_JOB_ID"]
        if "PBS_JOBID" in os.environ.keys():
            jobid = os.environ["PBS_JOBID"]
        try:
            data = json.dumps(state.data)
        except Exception:
            data = json.dumps(str(state.data))
        msg = {"id": self._tool.id,
               "state": state.name,
               "data": data,
               "jobid": jobid,
               "host": platform.uname()[1],
               "version": self._tool.version
               }

        client = None
        try:
            self._tool.logger.debug("发送消息%s" % msg)
            client = zerorpc.Client()
            client.connect(self.config.endpoint)
            result = client.report(msg)
            client.close()
        except Exception, e:
            exstr = traceback.format_exc()
            print >> sys.stderr, exstr
            sys.stderr.flush()
            if client:
                client.close()
            self._lost_connection_count += 1
            if self._lost_connection_count >= 10:
                self._tool.logger.error("网络连接出现错误，尝试10次仍然无法连接，即将退出运行:%s" % e)
                self._tool.exit_signal = True
                self._tool.exit(1)
            else:
                self._tool.logger.error("网络连接出现错误，将重新尝试连接:%s" % e)
                time.sleep(15)
                self.send_state(state)
        else:
            self._tool.logger.debug("返回消息%s" % result)
            if not isinstance(result, dict):
                self._tool.logger.error("接收到异常信息，退出运行!")
                self._tool.exit(1)
            self._lost_connection_count = 0
            return result


# class ProcessActor(RemoteActor):
#     def __init__(self, tool, main_thread):
#         super(ProcessActor, self).__init__(tool, main_thread)
#
#     def run(self):
#         self.config.KEEP_ALIVE_TIME = 1
#         super(ProcessActor, self).run()
#
#     def send_state(self, state):
#
#         msg = {"id": self._tool.id,
#                "state": state.name,
#                "data": state.data,
#                "version": self._tool.version
#                }
#         try:
#             self._tool.logger.debug("put msg %s" % msg)
#             self._tool.process_queue.put(msg)
#         except Exception, e:
#             self._tool.logger.debug("error: %s", e)
#
#         # print "Put MSG:%s" % msg
#         key = "%s" % self._tool.version
#         if key in self._tool.shared_callback_action.keys():
#             action = self._tool.shared_callback_action[key]
#             del self._tool.shared_callback_action[key]
#         else:
#             action = {'action': 'none'}
#         return action


class RemoteData(object):
    def __init__(self, parent):
        self._parent = parent

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def save(self):
        """
        保存数据到文件

        :return:
        """
        data = copy.copy(self.__dict__)
        data.pop("_parent")
        path = os.path.join(self._parent.work_dir, "remote_data.pk")
        with open(path, "w") as f:
            pickle.dump(data, f)

    def sync(self, tag="default"):
        """
        同步数据到远程agent,只支持tool
        :return:
        """
        self.save()
        self._parent.add_state("sync", data=tag)

    def load(self):
        """
        加载文件中保存的数据

        :return:
        """
        path = os.path.join(self._parent.work_dir, "remote_data.pk")
        if os.path.exists(path):
            with open(path, "r") as f:
                data = pickle.load(f)
                self.__dict__.update(data)
