# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from ..core.singleton import singleton
from multiprocessing import Event, Process, Array, Queue, Value
from .workflow import WorkflowWorker
from ..wsheet import Sheet
from .logger import Logger
from ..config import Config
from .db import WorkflowModel, CheckModel, ClientKeyModel
import os
import traceback
from pwd import getpwnam
from threading import Thread
import setproctitle
import time
import sys
from ..core.function import hostname, MessageData, ActionData, LogData, \
    add_run_queue, add_action_queue, StateData, CallbackData
import hashlib
import datetime
import re
import threading
import copy


@singleton
class WorkflowManager(object):
    def __init__(self):
        self.queue = Queue()
        self.action_queue = Array(ActionData, 200)
        self.workflows = {}
        self.event = {}
        self.multi_event = {}
        self.return_msg = {}
        self.logger = Logger()
        self.server = None
        self._key = {}
        self._port_list = {}
        self.config = Config()
        self.run_info_queue = Array(MessageData, 2000)
        self.log_info_queue = Array(LogData, 3000)
        self.end_list = []
        self.lock = threading.Lock()
        self.workflow_lock = threading.Lock()
        self.endpoint = "tcp://{}:{}".format(self.config.LISTEN_IP, self.config.get_listen_port())
        self.rpc_message_queue = Array(StateData, 3000)
        self.rpc_callback_queue = Array(CallbackData, 2000)
        self.rpc_restart_signal = Value("i", 0)
        self.e_mutex = threading.Lock()
        self.me_mutex = threading.Lock()
        self.r_mutex = threading.Lock()

    def add_task(self, json_data):
        try:
            if not isinstance(json_data, dict) or "id" not in json_data.keys():
                self.logger.error("add workflow %s format error!" % json_data)
                return {"success": False, "info": "json格式错误"}
            json_data["WPM"] = True
            # json_data["endpoint"] = "tcp://{}:{}".format(self.__get_ip(), self.__get_port(json_data["id"]))
            json_data["endpoint"] = self.endpoint
            wsheet = Sheet(data=json_data)
            model = WorkflowModel(wsheet)
            if wsheet.rerun is not True:
                if wsheet.id in self.workflows.keys():
                    self.logger.error("Workflow %s has already run! " % wsheet.id)
                    return {"success": False, "info": "Workflow %s正在运行，不能重复添加!" % wsheet.id}
                if model.find():
                    self.logger.error("workflow %s is already exists!" % wsheet.id)
                    return {"success": False, "info": "Workflow %s已经存在，不能重复运行!" % wsheet.id}
            with self.workflow_lock:
                self.workflows[wsheet.id] = model
            if wsheet.rerun is not True:
                m = re.match(r"^(.+)\.report\.", wsheet.name)
                if m:
                    model.save(report=m.group(1))
                else:
                    model.save()
            model.close()
            self.queue.put(json_data)
            self.logger.info("接收到Workflow[%s] 请求,放入队列... " % wsheet.id)
            return {"success": True, "info": "任务提交成功.", "workflow_id": wsheet.id}
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 add_task: %s" % e)

    def __get_port(self, wid):
            if wid in self._port_list.keys():
                return self._port_list[wid]
            port = self.config.get_listen_port()
            if port in self._port_list.values():
                return self.__get_port(wid)
            else:
                self._port_list[wid] = port
                return port

    def __get_ip(self):
        return self.config.LISTEN_IP

    def get_event(self, *wid):
        for eid in wid:
            if eid not in self.workflows.keys():
                raise Exception("ID%s不存在，请先添加任务!" % eid)
        if len(wid) == 1:
            my_id = wid[0]
        else:
            ids = list(wid)
            ids.sort()
            my_id = hashlib.sha1("".join([str(i) for i in ids])).hexdigest()

        if my_id in self.event.keys():
            return self.event[my_id]
        else:
            event = Event()
            with self.e_mutex:
                self.event[my_id] = event
            if len(wid) > 1:
                for i in wid:
                    with self.me_mutex:
                        self.multi_event[my_id] = {}
                        self.multi_event[my_id][i] = False
            return event

    def start(self, wid, data):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].update_pid(data["pid"], data["work_dir"])
                self.workflows[wid].close()
                self.logger.info("Workflow %s 开始运行! " % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 start: %s" % e)

    def rerun(self, wid, pid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].update_rerun(pid)
                self.workflows[wid].close()
                self.logger.info("开始重新运行Workflow %s , PID: %s! " % (wid, pid))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 start: %s" % e)

    def set_end(self, wid, msg=None):
        try:
            if wid in self.workflows.keys():
                with self.workflow_lock:
                    model = self.workflows.pop(wid)
                with self.lock:
                    if wid in self.end_list:
                            self.end_list.remove(wid)
            else:
                wsheet = Sheet(data={"id": wid})
                model = WorkflowModel(wsheet)
            if model:
                model.end()
                model.close()
            if wid in self.event.keys():
                with self.r_mutex:
                    self.return_msg[wid] = {"success": True, "info": msg}
                with self.e_mutex:
                    event = self.event.pop(wid)
                if event:
                    event.set()
            # if wid in self._port_list.keys():
            #     self._port_list.pop(wid)
            for m_id in self.multi_event.keys():
                all_over = True
                for mw_id in self.multi_event[m_id]:
                    if mw_id == wid:
                        with self.me_mutex:
                            self.multi_event[m_id][mw_id] = {"success": True, "info": "运行结束!"}
                    if self.multi_event[m_id][mw_id] is False:
                        all_over = False
                if all_over:
                    if m_id in self.event:
                        with self.e_mutex:
                            event = self.event.pop(m_id)
                        if event:
                            event.set()
            self.logger.info("Workflow %s 运行结束! " % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 set_end: %s" % e)

    def process_end(self, *wid):
        try:
            for w_id in wid:
                # if isinstance(w_id, tuple) or isinstance(w_id, list):
                m = re.match(r'([\w\-]+):(.*)', w_id)
                if m:
                    error_msg = "%s运行错误:%s " % (w_id, m.group(2))
                    w_id = m.group(1)
                else:
                    error_msg = "%s进程意外结束.." % w_id
                if wid in self.end_list:
                    continue
                if w_id in self.workflows.keys():
                    with self.workflow_lock:
                        model = self.workflows.pop(w_id)

                    if model:
                        model.error(error_msg)
                        model.close()

                    # if w_id in self._port_list.keys():
                    #     self._port_list.pop(w_id)

                    if w_id in self.event.keys():
                        with self.r_mutex:
                            self.return_msg[w_id] = {"success": False, "info": error_msg}
                        with self.e_mutex:
                            event = self.event.pop(w_id)
                        if event:
                            event.set()
                    for m_id in self.multi_event.keys():
                        all_over = True
                        for mw_id in self.multi_event[m_id]:
                            if mw_id == w_id:
                                with self.me_mutex:
                                    self.multi_event[m_id][mw_id] = {"success": False, "info": error_msg}
                            if self.multi_event[m_id][mw_id] is False:
                                all_over = False
                        if all_over:
                            if m_id in self.event:
                                with self.e_mutex:
                                    event = self.event.pop(m_id)
                                if event:
                                    event.set()
                    self.logger.error(error_msg)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 process_end: %s" % e)

    def get_msg(self, *wid):
        try:
            if len(wid) == 1:
                if wid[0] in self.return_msg.keys():
                    with self.r_mutex:
                        msg = self.return_msg.pop(wid[0])
                    return msg
                else:
                    return None
            else:
                ids = list(wid)
                ids.sort()
                my_id = hashlib.sha1("".join([str(i) for i in ids])).hexdigest()

                if my_id in self.multi_event.keys():
                    msg = {}
                    for mw_id in self.multi_event[my_id]:
                        if self.multi_event[my_id][mw_id] is not False:
                            msg[mw_id] = self.multi_event[my_id][mw_id]
                    with self.me_mutex:
                        self.multi_event.pop(my_id)
                    return msg
                else:
                    return None
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 get_msg: %s" % e)

    def keep_alive(self, wid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].update()
                self.logger.debug("接收到Workflow %s keepavlie信息! " % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 keep_alive: %s" % e)

    def set_error(self, wid, error_msg):
        try:
            if wid in self.workflows.keys():
                with self.workflow_lock:
                    model = self.workflows.pop(wid)
                if model:
                    model.error(error_msg)
                    model.close()

                # if wid in self._port_list.keys():
                #     self._port_list.pop(wid)

                if wid in self.event.keys():
                    with self.r_mutex:
                        self.return_msg[wid] = {"success": False, "info": error_msg}
                    with self.e_mutex:
                        event = self.event.pop(wid)
                    if event:
                        event.set()
                for m_id in self.multi_event.keys():
                    all_over = True
                    for mw_id in self.multi_event[m_id]:
                        if mw_id == wid:
                            with self.me_mutex:
                                self.multi_event[m_id][mw_id] = {"success": False, "info": error_msg}
                        if self.multi_event[m_id][mw_id] is False:
                            all_over = False
                    if all_over:
                        if m_id in self.event:
                            with self.e_mutex:
                                event = self.event.pop(m_id)
                            if event:
                                event.set()
                self.logger.error("Workflow %s 运行出错: %s " % (wid, error_msg))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 set_error: %s" % e)

    def set_pause(self, wid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].pause()
                self.workflows[wid].close()
                self.logger.info("Workflow %s 暂停运行..." % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 set_pause: %s" % e)

    def set_pause_exit(self, wid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].exit_pause()
                self.workflows[wid].close()
                self.logger.info("Workflow %s 退出暂停，继续运行..." % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 set_pause_exit: %s" % e)

    def pause_timeout(self, wid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].pause_timeout()
                self.workflows[wid].close()
                self.logger.info("Workflow %s 暂停超时，结束运行..." % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 pause_timeout: %s" % e)

    def set_stop(self, wid):
        try:
            if wid in self.workflows.keys():
                self.workflows[wid].stop()
                self.workflows[wid].close()
                self.logger.info("Workflow %s 接收中止运行指令..." % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("服务异常 set_stop: %s" % e)

    def get_key(self, client):
        if client not in self._key.keys():
            model = ClientKeyModel()
            self._key[client] = model.find_key(client)
            model.close()
        return self._key[client]


def get_event(*wid):
    try:
        if len(wid) == 1:
            event = WorkflowManager().get_event(wid[0])
        else:
            event = WorkflowManager().get_event(*wid)

    except Exception, e:
        exstr = traceback.format_exc()
        print e
        print exstr
        sys.stdout.flush()
        return None
    else:
        return event


class ManagerProcess(Process):
    def __init__(self, queue, action_queue, run_info_queue, log_info_queue, rpc_message_queue,
                 rpc_callback_queue, rpc_restart_signal):
        super(ManagerProcess, self).__init__()
        self.queue = queue
        self.queue.cancel_join_thread()
        self.config = Config()
        self.process = {}
        self.action_queue = action_queue
        self.run_info_queue = run_info_queue
        self.log_info_queue = log_info_queue
        self._batch_queue = []
        self.process_num = Value("i", 0)
        self.batch_queue_num = Value("i", 0)
        self.queue_num = Value("i", 0)
        self._exit = False
        self.logger = Logger(log_type="PROCESS_MANAGER")
        self._log_date = None
        self.rpc_message_queue = rpc_message_queue
        self.rpc_callback_queue = rpc_callback_queue
        self.rpc_restart_signal = rpc_restart_signal
        self.p_mutex = threading.Lock()
        self.mysql_lock = threading.Lock()

    def run(self):
        super(ManagerProcess, self).run()
        os.nice(5)
        os.setgid(getpwnam(self.config.wpm_user)[3])
        os.setuid(getpwnam(self.config.wpm_user)[2])
        setproctitle.setproctitle("WPM[Process Manager]")
        self.start_thread()
        max_process = self.config.get_wpm_limit(hostname)

        while True:
            try:
                added = False
                try:
                    self.process_num.value = len(self.process.keys())
                    self.batch_queue_num.value = len(self._batch_queue)
                    self.queue_num.value = self.queue.qsize()
                    json_data = self.queue.get_nowait()
                except Exception:
                    if os.getppid() == 1:
                        if self.process_num.value == 0 and self.batch_queue_num.value == 0:
                            self._exit = True
                            self.logger.info("检测到主线程结束，且所有任务完成，准备退出...")
                            break
                    time.sleep(1)
                else:
                    try:
                        wsheet = Sheet(data=json_data)
                        if wsheet.batch and wsheet.instant is not True:
                            self._batch_queue.append(wsheet)
                            added = True
                        else:
                            worker = WorkflowWorker(wsheet, self.action_queue, self.run_info_queue, self.log_info_queue,
                                                    self.rpc_message_queue, self.rpc_callback_queue,
                                                    self.rpc_restart_signal)
                            worker.start()
                            self.logger.info("创建进程PID%s, 运行Workflow %s ..." % (worker.pid, json_data["id"]))
                            with self.p_mutex:
                                self.process[json_data["id"]] = worker
                            if wsheet.rerun:
                                self.process_rerun(json_data["id"], worker.pid)
                            else:
                                self.process_start(json_data["id"], worker.pid, worker.work_dir)
                    except Exception, e:
                        exstr = traceback.format_exc()
                        print exstr
                        print e
                        sys.stdout.flush()
                        sys.stderr.flush()
                        self.process_end("%s:%s" % (json_data["id"], e))

                if len(self.process.keys()) <= max_process:
                    if len(self._batch_queue) > 0:
                        sheet = self._batch_queue.pop(0)
                        try:
                            worker = WorkflowWorker(sheet, self.action_queue, self.run_info_queue, self.log_info_queue,
                                                    self.rpc_message_queue, self.rpc_callback_queue,
                                                    self.rpc_restart_signal)
                            worker.start()
                            self.logger.info("创建进程PID%s, 运行一键化任务 %s ..." % (worker.pid, sheet.id))
                            with self.p_mutex:
                                self.process[sheet.id] = worker
                            if sheet.rerun:
                                self.process_rerun(sheet.id, worker.pid)
                            else:
                                self.process_start(sheet.id, worker.pid, worker.work_dir)
                        except Exception, e:
                            exstr = traceback.format_exc()
                            print exstr
                            print e
                            sys.stdout.flush()
                            sys.stderr.flush()
                            self.process_end("%s:%s" % (sheet.id, e))
                else:
                    if added:
                        self.logger.info("正在运行的任务数为%s,当前排队一键化任务数为%s" %
                                         (len(self.process.keys()), len(self._batch_queue)))
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()
                sys.stderr.flush()

    def process_end(self, *wid):
        # self.run_info_queue.put(("process_end", wid))
        # self.log_info_queue.put(("set_end", wid))
        # self._add_run_queue("process_end", ";".join(wid))
        # self._add_log_queue("set_end", ";".join(wid))
        add_run_queue(self.run_info_queue, "process_end", wid)
        # add_log_queue(self.log_info_queue, "set_end", wid)
        # class WorkerManager(BaseManager):
        #     pass
        # try:
        #     WorkerManager.register("worker")
        #     wpm_manager = WorkerManager(address=self.config.wpm_listen, authkey=self.config.wpm_authkey)
        #     wpm_manager.connect()
        #     worker = wpm_manager.worker()
        #     worker.process_end(*wid)
        #     del worker
        # except Exception, e:
        #     exstr = traceback.format_exc()
        #     print exstr
        #     print e
        #     sys.stdout.flush()

        # class LogManager(BaseManager):
        #     pass
        # try:
        #     LogManager.register("apilog")
        #     m = LogManager(address=self.config.wpm_logger_listen, authkey=self.config.wpm_logger_authkey)
        #     m.connect()
        #     log = m.apilog()
        #     log.set_end(*wid)
        #     del log
        # except Exception, e:
        #     exstr = traceback.format_exc()
        #     print exstr
        #     print e
        #     sys.stdout.flush()

    def process_start(self, wid, pid, work_dir):
        # class WorkerManager(BaseManager):
        #     pass
        #
        # WorkerManager.register("worker")
        # wpm_manager = WorkerManager(address=self.config.wpm_listen, authkey=self.config.wpm_authkey)
        # wpm_manager.connect()
        # worker = wpm_manager.worker()
        # worker.start(wid, pid)
        # del worker
        # self.run_info_queue.put(("start", (wid, pid)))
        # self._add_run_queue("start", wid, pid)
        add_run_queue(self.run_info_queue, "start", wid, {"pid": pid, "work_dir": work_dir})

    def process_rerun(self, wid, pid):
        add_run_queue(self.run_info_queue, "rerun", wid, pid)

    # def _add_action_queue(self, action, wid):
    #     # length = len(self.action_queue)
    #     time_stamp = int(time.time())
    #     is_full = True
    #     with self.action_queue.get_lock():
    #         for i, v in enumerate(self.action_queue):
    #             if v.id == wid and v.action:
    #                 return
    #         for i, v in enumerate(self.action_queue):
    #             if self.action_queue[i].action is not None:
    #                 if time_stamp - self.action_queue[i].timestamp > 200:
    #                     is_full = False
    #                     self.action_queue[i] = (action, wid, time_stamp)
    #                     break
    #             else:
    #                 is_full = False
    #                 self.action_queue[i] = (action, wid, time_stamp)
    #                 break
    #     if is_full:
    #         time.sleep(1)
    #         self._add_action_queue(action, wid)
    #         return

    def _check_process(self):
        wid_list = {}
        while True:
            if self._exit:
                break
            time.sleep(1)
            try:
                for wid, p in self.process.items():
                    if not p.is_alive():
                        p.join()
                        wid_list[wid] = datetime.datetime.now()
                        with self.p_mutex:
                            self.process.pop(wid)
                    else:
                        if p.has_start_run.value == 0:
                            if (datetime.datetime.now() - p.start_time).seconds > 300:
                                # 运行300秒未开始运行被卡住
                                if p.is_alive():
                                    # 强制终止
                                    self.logger.info("进程%s,worker %s，运行300秒尚未开始，可能被卡住，"
                                                     "尝试强制结束并重新运行" % (p.pid, wid))
                                    os.system("kill -9 %s" % p.pid)
                                    if not p.is_alive():
                                        p.join()
                                        with self.p_mutex:
                                            self.process.pop(wid)
                                        # 重新加入队列
                                        self._batch_queue.append(p.wsheet)

                send_list = []
                for wid, find_time in wid_list.items():
                    if (datetime.datetime.now() - wid_list[wid]).seconds > 100:
                        wid_list.pop(wid)
                        send_list.append(wid)
                if len(send_list) > 0:
                    self.process_end(*send_list)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()

    def _check_running(self):
        while True:
            if self._exit:
                break
            time.sleep(60)
            try:
                for wid, p in self.process.items():
                    if not p.is_alive():
                        p.join()
                        with self.p_mutex:
                            self.process.pop(wid)
                with self.p_mutex:
                    running = self.process.keys()
                model = CheckModel()
                if running:
                    model.update_running(running)
                model.update_unknown()
                model.close()
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()
                continue

    def _check_stop(self):
        model = CheckModel()
        results = model.find_stop()
        if results:
            for row in results:
                wid = row["workflow_id"]
                if wid in self.process.keys():
                    # self.process[wid].action_queue.put("tostop")
                    # self._add_action_queue("tostop", wid)
                    add_action_queue(self.action_queue, "tostop", wid)
        model.close()

    def _check_pause(self):
        model = CheckModel()
        results = model.find_pause()
        if results:
            for row in results:
                wid = row["workflow_id"]
                if wid in self.process.keys():
                    # self.process[wid].action_queue.put("pause")
                    # self._add_action_queue("pause", wid)
                    add_action_queue(self.action_queue, "pause", wid)
        model.close()

    def _check_exit_pause(self):
        model = CheckModel()
        results = model.find_exit_pause()
        if results:
            for row in results:
                wid = row["workflow_id"]
                if wid in self.process.keys():
                    # self.process[wid].action_queue.put("exit_pause")
                    add_action_queue(self.action_queue, "exit_pause", wid)
        model.close()

    def _check(self):
        while True:
            if self._exit:
                break
            time.sleep(5)
            try:
                with self.mysql_lock:
                    self._check_stop()
                    self._check_pause()
                    self._check_exit_pause()
            except Exception, e:
                print e
                sys.stdout.flush()
                continue

    def start_thread(self):
        thread = Thread(target=self._check_process, args=(), name='thread-process_check')
        thread.setDaemon(False)
        thread.start()
        thread1 = Thread(target=self._check, args=(), name='thread-stop_pause_check')
        thread1.setDaemon(False)
        thread1.start()
        thread2 = Thread(target=self._check_date, args=(), name='thread-check_date')
        thread2.setDaemon(False)
        thread2.start()
        thread3 = Thread(target=self._check_running, args=(), name='thread-check_running')
        thread3.setDaemon(False)
        thread3.start()

    def _check_date(self):
        while True:
            if self._exit:
                break
            if self._log_date != time.strftime('%Y%m%d', time.localtime(time.time())):
                self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
                log = os.path.join(self.config.wpm_log_file, "%s.%s.log" % (self._log_date, hostname))
                so = file(log, 'a+')
                se = file(log, 'a+', 0)
                os.dup2(so.fileno(), sys.stdout.fileno())
                os.dup2(se.fileno(), sys.stderr.fileno())
            time.sleep(5)

# @singleton
# class ApiLogManager(object):
#     def __init__(self):
#         # self._workers = {}
#         self._running_workers = {}
#         self.logger = Logger(log_type="API_LOG_MANAGER[%s]" % os.getpid())
#         self.config = Config()
#
#     # @property
#     # def running_number(self):
#     #     i = 0
#     #     for key, worker in self._running_workers.items():
#     #         i += worker.running_number
#     #     return i
#
#     def clear_old(self):
#         for key, worker in self._running_workers.items():
#             time_now = int(time.time())
#             if time_now - worker.last_update_time > 3600:
#                 self._running_workers.pop(key)
#
#     def add_log(self, data):
#         try:
#             if data["api"] in self.config.update_exclude_api:
#                 self.logger.info("Worker %s  API: %s API在排除更新列表中,忽略..." % (data["task_id"], data["api"]))
#                 return
#             # if data["task_id"] in self._workers.keys():
#             #     self._workers[data["task_id"]].add_log(data)
#             #     self.logger.info("Worker %s  更新进度.. " % data["task_id"])
#             elif data["task_id"] in self._running_workers.keys():
#                 self._running_workers[data["task_id"]].add_log(data)
#                 self.logger.info("Worker %s  更新进度.. " % data["task_id"])
#             else:
#                 worker = LogWorker(data["task_id"])
#                 worker.add_log(data)
#                 self._running_workers[data["task_id"]] = worker
#                 # self._workers[data["task_id"]] = worker
#                 self.logger.info("新Worker %s  API: %s 开始运行... " % (data["task_id"], data["api"]))
#         except Exception, e:
#             exstr = traceback.format_exc()
#             print exstr
#             sys.stdout.flush()
#             self.logger.error("服务异常 add_log: %s" % e)

    # def get_worker(self):
    #     try:
    #         for wid, worker in self._workers.items():
    #             self._workers.pop(wid)
    #             if not worker.is_end:
    #                 self._running_workers[wid] = worker
    #             return worker
    #         return None
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         sys.stdout.flush()
    #         self.logger.error("服务异常 get_worker: %s" % e)

    # def set_end(self, *wid):
    #     try:
    #         for w_id in wid:
    #             # if w_id in self._workers.keys():
    #             #     self._workers[w_id].set_end()
    #             if w_id in self._running_workers.keys():
    #                 # self._running_workers[w_id].set_end()
    #                 self._running_workers.pop(w_id)
    #             self.logger.info("Worker %s 运行完成." % w_id)
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         sys.stdout.flush()
    #         self.logger.error("服务异常 set_end: %s" % e)
    #
    # def report(self, data):
    #     """
    #     保存运行数据
    #     :return:
    #     """
    #     try:
    #         wid = data["id"]
    #         model = ReportModel(wid)
    #         if len(data["modules"]) > 0:
    #             model.save_modules(data["modules"])
    #         if len(data["tools"]) > 0:
    #             model.save_tools(model.auto_id, data["tools"], commit=True, under_workflow=1)
    #         self.logger.info("Workflow %s 保存运行报告完成" % wid)
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         sys.stdout.flush()
    #         self.logger.error("服务异常 report: %s" % e)


# class ApiLogProcess(Process):
#     def __init__(self, job_queue, main_last_get, count=0):
#         super(ApiLogProcess, self).__init__()
#         self.config = Config()
#         self._log_date = None
#         self.job_queue = job_queue
#         self.last_update = Value("i", 0)
#         self.last_update.value = int(time.time())
#         self.exit_signal = Value("i", 0)
#         self.count = count
#         self.logger = None
#         self.last_get_info = int(time.time())
#         self._exit = False
#         self.is_busy = Value("i", 0)
#         self.main_last_get = main_last_get
#         self.log_manager = None
#
#     def _check_date(self):
#         while True:
#             if self._exit:
#                 break
#             self.log_manager.clear_old()
#             if self._log_date != time.strftime('%Y%m%d', time.localtime(time.time())):
#                 self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
#                 log = os.path.join(self.config.UPDATE_LOG, "%s.%s.log" % (self._log_date, hostname))
#                 if not os.path.exists(self.config.UPDATE_LOG):
#                     os.mkdir(self.config.UPDATE_LOG)
#                 so = file(log, 'a+')
#                 se = file(log, 'a+', 0)
#                 os.dup2(so.fileno(), sys.stdout.fileno())
#                 os.dup2(se.fileno(), sys.stderr.fileno())
#             gevent.sleep(60)
#
#
#     def run(self):
#         super(ApiLogProcess, self).run()
#         setproctitle.setproctitle("WPM[API Log Worker Process %s]" % self.count)
#         self.logger = Logger(log_type="Log Worker Process %s[PID:%s]" % (self.count, os.getpid()))
#         self.log_manager = ApiLogManager()
#         self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
#         log = os.path.join(self.config.UPDATE_LOG, "%s.%s.log" % (self._log_date, hostname))
#         if not os.path.exists(self.config.UPDATE_LOG):
#             os.mkdir(self.config.UPDATE_LOG)
#         so = file(log, 'a+')
#         se = file(log, 'a+', 0)
#         os.dup2(so.fileno(), sys.stdout.fileno())
#         os.dup2(se.fileno(), sys.stderr.fileno())
#         gevent.spawn(self._check_date)
#
#         while True:
#                 self.last_update.value = int(time.time())
#                 self.is_busy.value = 0
#                 with self.main_last_get.get_lock():
#                     self.main_last_get.value = int(time.time())
#                 try:
#                     data = self.job_queue.get_nowait()
#                 except:
#                     if int(time.time()) - self.last_get_info > 600:
#                         self.logger.info("超过10分钟未能获得日志信息，退出工作进程%s..." % self.name)
#                         break
#                     if self.exit_signal.value == 1:
#                         self.logger.info("检测到服务退出，日志处理完毕,准备退出..")
#                         gevent.sleep(60)
#                         self._exit = True
#                         break
#                     gevent.sleep(1)
#                 else:
#                     try:
#                         self.last_get_info = int(time.time())
#                         if data[0] == "set_end":
#                             id_list = data[1]
#                             if isinstance(id_list, list) or isinstance(id_list, tuple):
#                                 self.log_manager.set_end(*id_list)
#                             else:
#                                 self.log_manager.set_end(id_list)
#                         elif data[0] == "add_log":
#                             self.log_manager.add_log(data[1])
#                         elif data[0] == "report":
#                             self.log_manager.report(data[1])
#                             self.log_manager.set_end(data[1]["id"])
#                     except Exception, e:
#                         exstr = traceback.format_exc()
#                         print exstr
#                         sys.stdout.flush()
#                         self.logger.error("日志更新错误数据: %s, 错误信息: %s " % (json.dumps(data), e))
#         # while True:
#         #
#         #     try:
#         #         job_list = copy.deepcopy(self.job_list)
#         #         self.job_list = []
#         #         for job in job_list:
#         #             if job[0] == "set_end":
#         #                 id_list = job[1]
#         #                 if isinstance(id_list, list) or isinstance(id_list, tuple):
#         #                     log_manager.set_end(*id_list)
#         #                 else:
#         #                     log_manager.set_end(id_list)
#         #             elif job[0] == "add_log":
#         #                 log_manager.add_log(job[1])
#         #             elif job[0] == "report":
#         #                 log_manager.report(job[1])
#         #     except Exception, e:
#         #         exstr = traceback.format_exc()
#         #         print exstr
#         #         print e
#         #         sys.stdout.flush()
#         #     gevent.sleep(1)
