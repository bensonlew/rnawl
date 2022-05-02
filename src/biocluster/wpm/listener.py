# -*- coding: utf-8 -*-
# __author__ = 'guoquan'


from .manager import get_event, ManagerProcess, WorkflowManager
from ..config import Config
import setproctitle
import time
import os
import sys
from threading import Thread
from ..patch.manager import NewBaseManager
import traceback
from ..core.singleton import singleton
import datetime
from ..core.function import hostname, get_run_queue
import inspect
import threading
import copy
import signal
from pwd import getpwnam
from .log import ApiLogManagerProcess
from .rpc import RpcServerProcess


def write_log(info):
    print("%s    %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), info))
    sys.stdout.flush()


def start():
    server = MainServer()
    server.start()


@singleton
class MainServer(object):
    def __init__(self):
        self.workflow_manager = WorkflowManager()
        self.manager_server = ManagerProcess(self.workflow_manager.queue,
                                             self.workflow_manager.action_queue,
                                             self.workflow_manager.run_info_queue,
                                             self.workflow_manager.log_info_queue,
                                             self.workflow_manager.rpc_message_queue,
                                             self.workflow_manager.rpc_callback_queue,
                                             self.workflow_manager.rpc_restart_signal
                                             )
        self.api_log_server = ApiLogManagerProcess(self.workflow_manager.log_info_queue)
        self.rpc_server = RpcServerProcess(self.workflow_manager.endpoint,
                                           self.workflow_manager.rpc_message_queue,
                                           self.workflow_manager.rpc_callback_queue,
                                           self.workflow_manager.rpc_restart_signal)
        self.config = Config()
        self._log_date = None
        self.job_list = []
        self.lock = threading.Lock()
        self._task_threads = []
        self._task_thread_count = 0
        self._exit = False
        self._get_list_exit = False
        self.server = None
        self._last_get = int(time.time())
        self.observer = None
        self.to_reload = True

    def _check_date(self):
        while True:
            if self._exit:
                break
            if len(self.job_list) > 30:
                write_log("状态信息排队过长，生成新的处理线程...")
                self._add_task_thread()

            if len(self._task_threads) > 0 and len(self.job_list) > 0 and (int(time.time()) - self._last_get > 30):
                write_log("状态信息超过30秒未能更新...")
                self._add_task_thread()

            if self._log_date != time.strftime('%Y%m%d', time.localtime(time.time())):
                self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
                log = os.path.join(self.config.wpm_log_file, "%s.%s.log" % (self._log_date, hostname))
                so = file(log, 'a+')
                se = file(log, 'a+', 0)
                os.dup2(so.fileno(), sys.stdout.fileno())
                os.dup2(se.fileno(), sys.stderr.fileno())
                os.chown(log, getpwnam(self.config.wpm_user)[2], getpwnam(self.config.wpm_user)[2])
            time.sleep(5)

    def _add_task_thread(self):
        thread1 = Thread(target=self._workflow_process, args=(), name='thread-state_process_%s'
                                                                      % self._task_thread_count)
        thread1.setDaemon(False)
        thread1.start()
        self._task_thread_count += 1
        self._task_threads.append(thread1)
        write_log("创建状态处理线程%s...,当前共%s个状态处理线程" %
                  (thread1.getName(), len(self._task_threads)))

    def _get_list(self):
        while True:
            # if len(self.job_list) > 30:
            #     write_log("状态信息排队过长，生成新的处理线程...")
            #     self._add_task_thread()

            data = get_run_queue(self.workflow_manager.run_info_queue)
            with self.lock:
                if data:
                    self.job_list.extend(data)
            if len(data) > 0:
                for d in data:
                    if d[0] == "set_end":
                        with self.workflow_manager.lock:
                            self.workflow_manager.end_list.append(d[1])
                if len(self._task_threads) == 0:
                    self._add_task_thread()
            else:
                if self._exit:
                    if self.manager_server.queue_num.value == 0 and self.manager_server.batch_queue_num.value == 0 \
                            and self.manager_server.process_num.value == 0:
                        self._get_list_exit = True
                        has_end = True
                        for th in self._task_threads:
                            if th.is_alive():
                                has_end = False
                            else:
                                th.join()
                        if has_end:
                            write_log("状态信息处理完成，准备退出...")
                            os.system("kill -9 %s" % os.getpid())
                        else:
                            continue
                        break
                ths = copy.copy(self._task_threads)
                for th in ths:
                    if not th.is_alive():
                        th.join()
                        self._task_threads.remove(th)
            time.sleep(1)

    def _workflow_process(self):
        # length = len(self.workflow_manager.run_info_queue)
        last_get = datetime.datetime.now()
        while True:
            try:
                with self.lock:
                    job_list = copy.deepcopy(self.job_list)
                    self.job_list = []
                self._last_get = int(time.time())
                if len(job_list) > 0:
                    last_get = datetime.datetime.now()
                    for job in job_list:
                        if job[0] == "process_end":
                            id_list = job[1]
                            func = getattr(self.workflow_manager, job[0])
                            if isinstance(id_list, list) or isinstance(id_list, tuple):
                                func(*id_list)
                            else:
                                func(id_list)
                        else:
                            if hasattr(self.workflow_manager, job[0]):
                                func = getattr(self.workflow_manager, job[0])
                                argspec = inspect.getargspec(func)
                                args = argspec.args
                                if inspect.ismethod(func):
                                    args.pop(0)
                                length = len(args)
                                if length == 0:
                                    func()
                                elif length == 1:
                                    func(job[1])
                                else:
                                    func(job[1], job[2])
                else:
                    if self._get_list_exit:
                        write_log("状态处理线程 %s 退出" % threading.current_thread().getName())
                        break
                    if (datetime.datetime.now() - last_get).seconds > 600:
                        write_log("状态处理线程 %s 超过10分钟没有任务，退出线程" % threading.current_thread().getName())
                        break
                    time.sleep(1)
                if self.manager_server.pid and (not self.manager_server.is_alive()):
                    write_log("检测到进程管理器退出，准备结束运行")
                    os.system("kill -9 %s" % os.getpid())
                    break
                if self.api_log_server.pid and (not self.api_log_server.is_alive()):
                    write_log("检测到API Log管理器退出，准备结束运行")
                    os.system("kill -9 %s" % os.getpid())
                    break
                if self.rpc_server.pid and (not self.rpc_server.is_alive()):
                    write_log("检测到API Log管理器退出，准备结束运行")
                    os.system("kill -9 %s" % os.getpid())
                    break
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()

    def start_thread(self):
        thread = Thread(target=self._get_list, args=(), name='thread-get_queue_process')
        thread.start()
        thread1 = Thread(target=self._check_date, args=(), name='thread-check_date')
        thread1.start()

    def write_pids(self):
        main_pid_file = self.config.wpm_pid_dir + "/wpm.pid"
        main_pid = str(os.getpid())
        if not os.path.exists(self.config.wpm_pid_dir):
            os.mkdir(self.config.wpm_pid_dir)
        with open(main_pid_file, 'w+') as f:
            f.write('%s\n' % main_pid)

    def stop_listen(self, signum, frame):
        write_log("接收到终止服务器指令%s，准备关闭监听..." % signum)
        self._exit = True
        self.server.to_stop = True
        self.observer.stop()
        #  self.server.listener.close()

    def start(self):
        # start check thread
        setproctitle.setproctitle("WPM[Main Server]")
        os.nice(-5)
        self.start_thread()
        time.sleep(2)
        # start process manager
        write_log("启动进程管理器...")
        self.manager_server.start()

        # start api server
        write_log("启动API LOG监听...")
        self.api_log_server.start()

        # start main server
        write_log("启动WPM主服务监听...")
        self.write_pids()
        write_log("启动RPC服务监听...")
        self.rpc_server.start()

        class ListenerManager(NewBaseManager):
            pass
        ListenerManager.register('worker', WorkflowManager)
        ListenerManager.register('get_event', get_event)
        try:
            m = ListenerManager(address=self.config.wpm_listen, authkey=self.config.wpm_authkey)
            self.server = m.get_server()
            signal.signal(signal.SIGTERM, self.stop_listen)
            signal.signal(signal.SIGINT, self.stop_listen)
            self.server.server_until_stop()
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            sys.stderr.flush()
            self.manager_server.terminate()
            self.api_log_server.terminate()
            sys.exit(1)
        write_log("WPM主监听服务器退出...")
        sys.exit(0)
