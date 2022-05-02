# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

# import importlib
from ..core.function import CJsonEncoder
import traceback
from ..config import Config
import urllib2
import json
import urllib
from .logger import Logger
from .db import ApiLogModel
import random
import hashlib
import sys
import time
from threading import Thread
from ..core.function import get_clsname_form_path, get_log_queue, hostname
import importlib
from multiprocessing import Process
from .db import ReportModel
import copy
import os
import setproctitle
import socket
import re


class ApiLogManagerProcess(Process):
    def __init__(self, log_info_queue):
        super(ApiLogManagerProcess, self).__init__()
        self.log_info_queue = log_info_queue
        self.logger = Logger(log_type="Log Process")
        self.manager = LogManager(log_info_queue)
        self.config = Config()

    def run(self):
        super(ApiLogManagerProcess, self).run()
        setproctitle.setproctitle("WPM[API Log Process]")
        self.manager.run()


class LogManager(object):

    def __init__(self, info_queue):
        self.info_queue = info_queue
        self.log_object_queue = {}
        self.log_api_queue = {}
        self.upstate_queue = {}
        self.report_queue = {}
        self.log_api_threads = {}
        self.update_threads = {}
        self.report_threads = {}
        self.main_get_data_threads = {}
        self.logger = Logger(log_type="Log Manager")
        # self.log_object_lock = Lock()
        # self.log_api_lock = Lock()
        # self.upstate_lock = Lock()
        # self.report_lock = Lock()
        self.to_exit = False
        self.config = Config()
        self._log_date = None
        self.last_get = {}
        # self._logger = None
        self._ident_id = None
        self.running_progress = {}

    def get_queue(self):
        while True:
            try:
                data = get_log_queue(self.info_queue)
                self.last_get["main"] = int(time.time())
                # self.logger.info("接收到数据%s条.." % len(data))
                if len(data) > 0:
                    for d in data:
                        self.add_log(d)
                else:
                    if self.to_exit:
                        break
                    time.sleep(1)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                sys.stdout.flush()
                self.logger.error("运行异常: %s" % e)

    def add_log(self, data):
        data_type = data[0]
        log_data = data[1]

        if data_type == "report":
            # self.logger.info("接收到Report数据:%s" % json.dumps(log_data))
            try:
                self.report_queue[log_data["id"]] = log_data
                if len(self.report_threads) == 0:
                    self.logger.info("当前Report处理线程数为0")
                    self.add_thread("report")
                # self.logger.info("Report数据添加完成")
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
                self.logger.info("添加Report数据出错:%s" % e)
        else:
            try:
                self.logger.info("接收到Log数据:%s" % json.dumps(log_data))
                class_name = get_clsname_form_path(log_data["api"], tp="")
                api_name = "mbio.api.web.%s" % log_data["api"].lower()
                if api_name in sys.modules.keys():
                    main = sys.modules[api_name]
                else:
                    main = importlib.import_module(api_name)
                if hasattr(main, class_name):
                    api = getattr(main, class_name)
                    log = api(log_data)
                    r_key = "%s_%s" % (log.task_id, log.process_id)
                    if r_key in self.running_progress.keys() \
                            and log.progress < self.running_progress[r_key]:
                        pass
                    else:
                        self.running_progress[r_key] = log.progress
                        self.log_object_queue[id(log)] = log
                        self.upstate_queue[id(log)] = getattr(log, "_update_status")
                        self.log_api_queue[id(log)] = getattr(log, "_update_web_api")
                        if len(self.log_api_threads) == 0:
                            self.logger.info("当前Web API处理线程数为0")
                            self.add_thread("web")
                        if len(self.update_threads) == 0:
                            self.logger.info("Update State处理线程数为0")
                            self.add_thread("state")
                else:
                    self.logger.error("没有找到API模块:%s" % class_name)
                # self.logger.info("Log数据添加完成")
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
                self.logger.error("task_id: %s  api: %s  导入API模块出错:%s" % (log_data['task_id'], log_data["api"], e))

    def add_thread(self, thread_type):
        new_thread = None
        info = ""
        count = 0
        if thread_type == "report":
            new_thread = ReportThread(self, name="report_thread")
            info = "添加新运行报告处理线程"
            new_thread.start()
            self.report_threads[new_thread.ident] = new_thread
            count = len(self.report_threads)
        elif thread_type == "web":
            new_thread = WebThread(self, name="Web_api_thread")
            info = "添加新Web API处理线程"
            new_thread.start()
            self.log_api_threads[new_thread.ident] = new_thread
            count = len(self.log_api_threads)
        elif thread_type == "state":
            new_thread = StateThread(self, name="State_update_thread")
            info = "添加新状态更新处理线程"
            new_thread.start()
            self.update_threads[new_thread.ident] = new_thread
            count = len(self.update_threads)
        elif thread_type == "main_get_data":
            new_thread = Thread(target=self.get_queue, args=(), name="Main_get_data_thread")
            info = "添加新Mian get data线程"
            new_thread.start()
            self.main_get_data_threads[new_thread.ident] = new_thread
            count = len(self.main_get_data_threads) + 1
        self.logger.info("%s: %s, 当前工作线程数%s" % (info, new_thread.ident, count))

    def run(self):
        self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
        log = os.path.join(self.config.UPDATE_LOG, "%s.%s.log" % (self._log_date, hostname))
        if not os.path.exists(self.config.UPDATE_LOG):
            os.mkdir(self.config.UPDATE_LOG)
        so = file(log, 'a+')
        se = file(log, 'a+', 0)
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())
        end_thread = EndThread(self, name="End_thread")
        end_thread.start()
        check_thread = CheckTreahd(self, name="Check_thread")
        check_thread.start()
        self.get_queue()


class CheckTreahd(Thread):
    def __init__(self, manager, **kwargs):
        super(CheckTreahd, self).__init__(**kwargs)
        self.manager = manager
        self.logger = Logger(log_type="Check Worker")
        self._log_date = None
        self.config = Config()

    def run(self):
        super(CheckTreahd, self).run()
        thread_dict = {
            "web": self.manager.log_api_threads,
            "report": self.manager.report_threads,
            "state": self.manager.update_threads,
            "main_get_data": self.manager.main_get_data_threads
        }
        while True:
            try:
                if "main" in self.manager.last_get.keys() and (int(time.time()) - self.manager.last_get["main"]) > 10:
                    self.logger.info("超过10秒没有获取数据，开启新的主数据处理线程..")
                    self.manager.add_thread("main_get_data")
                if "Report Worker" in self.manager.last_get.keys() \
                        and (int(time.time()) - self.manager.last_get["Report Worker"]) > 10 \
                        and len(self.manager.report_threads) > 0:
                    self.logger.info("超过10秒没有获取Report数据，开启新的Report处理线程..")
                    self.manager.add_thread("report")
                if "Web Api Worker" in self.manager.last_get.keys() \
                        and (int(time.time()) - self.manager.last_get["Web Api Worker"]) > 10 \
                        and len(self.manager.log_api_threads) > 0:
                    self.logger.info("超过10秒没有获取Web API数据，开启新的Web API处理线程..")
                    self.manager.add_thread("web")
                if "State Worker" in self.manager.last_get.keys() \
                        and (int(time.time()) - self.manager.last_get["State Worker"]) > 10 \
                        and len(self.manager.update_threads) > 0:
                    self.logger.info("超过10秒没有获取Update State数据，开启新的Update State处理线程..")
                    self.manager.add_thread("state")
                if len(self.manager.log_api_queue) > 30:
                    self.logger.info("Web API 队列长队过长，开启新处理线程..")
                    self.manager.add_thread("web")
                if len(self.manager.report_queue) > 30:
                    self.logger.info("Report 队列长队过长，开启新处理线程..")
                    self.manager.add_thread("report")
                if len(self.manager.upstate_queue) > 30:
                    self.logger.info("Update state 队列长队过长，开启新处理线程..")
                    self.manager.add_thread("state")
                for key, q_thread in thread_dict.items():
                    for tid, p in q_thread.items():
                        if not p.is_alive():
                            self.logger.info("%s处理线程%s[%s]结束!" % (key, p.name, p.ident))
                            q_thread.pop(tid)
                            p.join()
                if os.getppid() == 1:
                    self.logger.info("检测到服务退出，日志处理完毕,准备退出")
                    self.manager.to_exit = True
                    break
                if self._log_date != time.strftime('%Y%m%d', time.localtime(time.time())):
                    self._log_date = time.strftime('%Y%m%d', time.localtime(time.time()))
                    log = os.path.join(self.config.UPDATE_LOG, "%s.%s.log" % (self._log_date, hostname))
                    if not os.path.exists(self.config.UPDATE_LOG):
                        os.mkdir(self.config.UPDATE_LOG)
                    so = file(log, 'a+')
                    se = file(log, 'a+', 0)
                    os.dup2(so.fileno(), sys.stdout.fileno())
                    os.dup2(se.fileno(), sys.stderr.fileno())
                time.sleep(10)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()
                self.logger.error("运行出错:%s" % e)


class EndThread(Thread):
    def __init__(self, manager, **kwargs):
        super(EndThread, self).__init__(**kwargs)
        self.manager = manager
        self.logger = Logger(log_type="End Worker")

        self.config = Config()

    def run(self):
        super(EndThread, self).run()
        while True:
            try:
                log_queue = copy.copy(self.manager.log_object_queue.values())
                end_list = []
                for log in log_queue:
                    if log.is_end:
                        end_list.append(log)
                        self.manager.log_object_queue.pop(id(log))
                        r_key = "%s_%s" % (log.task_id, log.process_id)
                        if log.status == "failed" or log.status == "finish" or \
                                (r_key in self.manager.running_progress.keys() and
                                         self.manager.running_progress[r_key] == 100):
                            self.manager.running_progress.pop(r_key)
                    elif log.has_update_webapi is True and log.webapi_end is False and log.web_api_success is False:
                        if log.last_webapi_update and \
                                        (time.time() - log.last_webapi_update) > log.config.UPDATE_RETRY_INTERVAL:
                            self.manager.log_api_queue[id(log)] = getattr(log, "_update_web_api")

                if len(end_list) > 0:
                    for l in end_list:
                        l.model.save()
                        l.model.close()
                        l.logger.info("所有任务运行完成!")
                else:
                    time.sleep(1)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()


class WorkerThread(Thread):
    def __init__(self, manager, **kwargs):
        super(WorkerThread, self).__init__(**kwargs)
        self.logger_name = None
        self.manager = manager
        self._last_get_message = None
        self._logger = None

    @property
    def queue(self):
        return self.manager.log_api_queue

    @property
    def logger(self):
        if not self._logger:
            self._logger = Logger(log_type="%s[%s]" % (self.logger_name, self.ident))
        return self._logger

    def _get_data(self):
        data = []
        # data_queue = copy.copy(self.queue.values())
        # self.queue = {}
        for key, value in self.queue.items():
            data.append(value)
            self.queue.pop(key)
        self.manager.last_get[self.logger_name] = int(time.time())
        if len(data) > 0:
            self._last_get_message = int(time.time())
        return data

    def _check_exit(self):
        if self.manager.to_exit and len(self.queue) == 0:
            self.logger.debug("检测到主线程退出,准备退出运行!")
            return True
        if self._last_get_message and (int(time.time()) - self._last_get_message) > 600 and len(self.queue) == 0:
            self.logger.debug("超过10分钟没有获取到新任务，准备退出运行!")
            return True
        return False

    def run(self):
        super(WorkerThread, self).run()

        while True:
            try:
                data = self._get_data()
                # self.logger.info("%s接收到数据%s条.." % (self.__class__.__name__, len(data)))
                if len(data) > 0:
                    for d in data:
                            d()
                else:
                    if self._check_exit():
                        break
                    time.sleep(1)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()


class ReportThread(WorkerThread):

    def __init__(self, *args, **kwargs):
        super(ReportThread, self).__init__(*args, **kwargs)
        self.logger_name = "Report Worker"

    @property
    def queue(self):
        return self.manager.report_queue

    def report(self, data):
        """
        保存运行数据
        :return:
        """
        wid = data["id"]
        self.logger.info("开始保存%s报告..." % wid)
        model = ReportModel(wid)
        try:
            if "cpu_used" in data.keys() and "memory_used" in data.keys():
                if "rerun" in data.keys() and data["rerun"] is True:
                    model.save_workflow(data["cpu_used"], data["memory_used"], data["error_data"], True)
                else:
                    model.save_workflow(data["cpu_used"], data["memory_used"], data["error_data"])

            if len(data["modules"]) > 0:
                model.save_modules(data["modules"])
            if len(data["tools"]) > 0:
                model.save_tools(wid, data["tools"], commit=True, under_workflow=1)
            self.logger.info("Workflow %s 保存运行报告完成" % wid)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("保存运行报告异常: %s" % e)
        finally:
            model.close()

    def run(self):
        super(WorkerThread, self).run()
        while True:
            try:
                data = self._get_data()
                # self.logger.info("%s接收到数据%s条.." % (self.__class__.__name__, len(data)))
                if len(data) > 0:
                    for d in data:
                        self.report(d)
                else:
                    if self._check_exit():
                        break
                    time.sleep(1)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr, e
                sys.stdout.flush()


class WebThread(WorkerThread):

    def __init__(self, *args, **kwargs):
        super(WebThread, self).__init__(*args, **kwargs)
        self.logger_name = "Web Api Worker"

    @property
    def queue(self):
        return self.manager.log_api_queue


class StateThread(WorkerThread):

    def __init__(self, *args, **kwargs):
        super(StateThread, self).__init__(*args, **kwargs)
        self.logger_name = "State Worker"

    @property
    def queue(self):
        return self.manager.upstate_queue


class Log(object):

    def __init__(self, data):
        self._data = data["data"]
        self.task_id = data["task_id"]
        self.api = data["api"]
        self.process_id = data["process_id"] if "process_id" in self.data.keys() else 0
        self._last_update = None
        self._response = None
        self._response_code = 0
        self.failed_times = 0
        self._url = ""
        # self._post_data = ""
        self.config = Config()
        self.logger = Logger(log_type="Log[%s: %s] " % (self.task_id, json.dumps(self.data["sync_task_log"]["log"])))
        self.update_info = self._loads_update_info(data["update_info"]) if "update_info" in data.keys() else None
        self._model = None
        self._client = "client01"
        self._key = "1ZYw71APsQ"
        self._url = "http://www.sanger.com/api/add_task_log"
        self.has_update_status = False
        self.update_status_success = False
        self.has_update_webapi = False
        self.web_api_success = False
        self.webapi_end = False
        self.update_status_end = False
        self.last_webapi_update = None
        self._project_type = self.__check_project_type()
        self._no_save = data["no_save"] if "no_save" in data.keys() else False

    @property
    def progress(self):
        if "task" in self.data["sync_task_log"].keys() and "progress" in self.data["sync_task_log"]["task"].keys():
            return int(self.data["sync_task_log"]["task"]["progress"])
        else:
            return 0

    @property
    def status(self):
        if "task" in self.data["sync_task_log"].keys() and "status" in self.data["sync_task_log"]["task"].keys():
            return self.data["sync_task_log"]["task"]["status"]

    def __check_project_type(self):
        module_name = self.__module__
        mlist = module_name.split(".")
        mlist.pop()
        m = re.match(r'mbio.api.web.([\w_]+)', ".".join(mlist))
        if m:
            return m.group(1)
        else:
            return None

    @property
    def _mongo_client(self):
        return self.config.get_mongo_client(self._project_type)

    @property
    def db(self):
        return self._mongo_client[self.config.get_mongo_dbname(self._project_type)]

    @property
    def model(self):
        if not self._model:
            self._model = ApiLogModel(self)
        return self._model

    def _loads_update_info(self, data):
        if data == "" or data == "null" or data is None:
            return None
        try:
            info = json.loads(data)
        except Exception, e:
            self.logger.info("update_info信息%s格式不正确:%s" % (data, e))
            return None
        else:
            return info

    @property
    def is_end(self):
        if self.webapi_end and self.update_status_end:
            return True
        else:
            return False

    @property
    def data(self):
        return self._data

    @property
    def response(self):
        return self._response

    @property
    def response_code(self):
        return self._response_code

    @property
    def post_data(self):
        my_content = self.data["sync_task_log"]
        my_data = dict()
        if "files" in my_content.keys() and len(my_content["files"]) > 5000:
            files_len = len(my_content["files"])
            base_info = dict()
            if "task" in my_content.keys():
                base_info["task"] = my_content["task"]
            if "log" in my_content.keys():
                base_info["log"] = my_content["log"]
            if "base_path" in my_content.keys():
                base_info["base_path"] = my_content["base_path"]
            dirs = {}
            if "dirs" in my_content.keys():
                # split_data = {
                #     "dirs": my_content["dirs"]
                # }
                # split_data.update(base_info)
                # my_data["sync_task_log"] = json.dumps(split_data, cls=CJsonEncoder)
                # yield urllib.urlencode(my_data)
                dirs = my_content["dirs"]
            start_index = 0
            if files_len % 1000 > 0:
                count = files_len/1000 + 1
            else:
                count = files_len/1000
            for i in xrange(count):
                end_index = start_index + 1000
                if end_index > files_len:
                    split_data = {
                        "files": my_content["files"][start_index:]
                    }
                else:
                    split_data = {
                        "files": my_content["files"][start_index:end_index]
                    }
                if dirs:
                    split_dirs = []
                    for f in split_data["files"]:
                        for d in dirs:
                            rel_path = os.path.relpath(f["path"], d["path"])
                            if not rel_path.startswith("../"):
                                split_dirs.append(d)
                        # if "path" in f.keys():
                        #     dir_name = os.path.dirname(f["path"])
                        #     if dir_name in dirs.keys():
                        #         split_dirs[dir_name] = dirs[dir_name]
                    split_data["dirs"] = split_dirs
                start_index = end_index
                split_data.update(base_info)
                my_data["sync_task_log"] = json.dumps(split_data, cls=CJsonEncoder)
                yield urllib.urlencode(my_data)
        else:
            my_data["sync_task_log"] = json.dumps(my_content, cls=CJsonEncoder)
            yield urllib.urlencode(my_data)

    def save(self):
        if self._no_save:
            return
        self.model.save()
        self.model.close()

    def send(self):
        self.has_update_webapi = True
        for p_data in self.post_data:
            http_handler = urllib2.HTTPHandler(debuglevel=1)
            https_handler = urllib2.HTTPSHandler(debuglevel=1)
            opener = urllib2.build_opener(http_handler, https_handler)
            urllib2.install_opener(opener)
            post_data = "%s&%s" % (self.get_sig(), p_data)
            request = urllib2.Request(self._url, post_data)
            response = urllib2.urlopen(request, timeout=60)
            self._response_code = response.getcode()
            self._response = response.read()
            response.close()
            self.logger.info("Return page:%s" % self._response)
            sys.stdout.flush()
        return self._response

    def _update_web_api(self):
        if not hasattr(self, "update"):
            self.webapi_end = True
            return
        try:
            self.last_webapi_update = int(time.time())
            getattr(self, "update")()
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.failed_times += 1
            if self.failed_times < self.config.UPDATE_MAX_RETRY:
                self.logger.error(" 更新Web api状态出错:%s,重试..." % e)
            else:
                self.webapi_end = True
                self.logger.error(" 更新Web api状态出错:%s,超过最大错误次数，终止尝试...." % e)
        else:
            if not self.has_update_webapi:
                self.webapi_end = True
                return
            if self.response:
                try:
                    response_json = json.loads(self.response)
                except Exception, e:
                    exstr = traceback.format_exc()
                    print exstr
                    sys.stdout.flush()
                    self.failed_times += 1
                    if self.failed_times < self.config.UPDATE_MAX_RETRY:
                        self.logger.error("提交失败: 提交发送错误 %s ，%s后重试..." %
                                          (e, self.config.UPDATE_RETRY_INTERVAL))
                    else:
                        self.webapi_end = True
                        self.logger.error("提交失败:提交发送错误 %s ，超过最大错误次数，终止尝试..." % e)
                else:
                    self.webapi_end = True
                    if response_json["success"] == "true" \
                            or response_json["success"] is True or response_json["success"] == 1:
                        self.web_api_success = True
                        self.logger.info("提交Web api成功")
                    else:
                        self.failed_times += 1
                        self.logger.error("提交被拒绝，终止提交:%s" % response_json["message"])
            else:
                self.webapi_end = True

    def get_sig(self):
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [self._key, timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "client": self._client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig
        }
        return urllib.urlencode(signature)

    def _update_status(self):
        if not hasattr(self, "update_status"):
            self.update_status_end = True
            return
        self.has_update_status = True
        try:
            getattr(self, "update_status")()
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            self.logger.error("更新数据库状态出错:%s" % e)
        else:
            self.logger.info("更新数据库状态成功!")
            self.update_status_success = True
        self.update_status_end = True
