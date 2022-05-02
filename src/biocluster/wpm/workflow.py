# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from multiprocessing import Process, Value
from ..core.function import load_class_by_path, hostname, add_log_queue, add_run_queue, filter_error_info,\
    get_clsname_form_path
import os
import json
import traceback
import datetime
from ..config import Config
import time
import sys
import setproctitle
from .logger import Logger
from ..core.function import CJsonEncoder
# from .client import worker_client, log_client
import re
from ..core.exceptions import CodeError


class WorkflowWorker(Process):
    """
    流程新建进程
    """
    def __init__(self, wsheet, action_queue, run_info_queue, log_info_queue, rpc_message_queue,
                 rpc_callback_queue, rpc_restart_signal):
        """

        :param wsheet:  运行流程所需的Sheet对象
        :param log_queue:
        :param action_queue:
        """
        super(WorkflowWorker, self).__init__()
        self.wsheet = wsheet
        self.wsheet.action_queue = action_queue
        self.wsheet.run_info_queue = run_info_queue
        self.wsheet.log_info_queue = log_info_queue
        self.config = Config()
        self.logger = Logger(log_type="Worker[%s]" % wsheet.id)
        self.wsheet.rpc_message_queue = rpc_message_queue
        self.wsheet.rpc_callback_queue = rpc_callback_queue
        self.wsheet.rpc_restart_signal = rpc_restart_signal
        self._work_dir = None
        self.has_start_run = Value("i", 0)
        self.wsheet.has_start_run = self.has_start_run
        self.start_time = datetime.datetime.now()

    @property
    def work_dir(self):
        if not self._work_dir:
            timestr = str(time.strftime('%Y%m%d', time.localtime(time.time())))
            self._work_dir = self.config.WORK_DIR + "/" + timestr + "/" + self._get_min_name() + "_" + self.wsheet.id
        return self._work_dir

    def _get_min_name(self):
        if self.wsheet.type != "workflow":
            return "Single"
        class_name = get_clsname_form_path(self.wsheet.name)
        base = ['Basic', 'Module', 'Tool', 'Agent', 'Workflow']
        for b in base:
            if re.search((b + "$"), class_name):
                # return re.sub((b + "$"), '', class_name).lower()
                return re.sub((b + "$"), '', class_name)
        return class_name

    def run(self):
        super(WorkflowWorker, self).run()
        setproctitle.setproctitle("WPM[worker %s]" % self.wsheet.id)
        # 输出重定向
        timestr = time.strftime('%Y%m%d', time.localtime(time.time()))
        log_dir = os.path.join(self.config.wpm_log_file, timestr)
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        log = os.path.join(log_dir, "%s_%s.log" % (self.wsheet.id, hostname))
        so = file(log, 'a+')
        se = file(log, 'a+', 0)
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())
        workflow = None
        try:
            if self.wsheet.type == "workflow":
                path = self.wsheet.name
            else:
                path = "single"

            workflow_module = load_class_by_path(path, "Workflow")
            workflow = workflow_module(self.wsheet)
            self.logger.debug("模块初始化成功，开始运行...")
            with self.has_start_run.get_lock():
                self.has_start_run.value = 1
            file_path = os.path.join(workflow.work_dir, "data.json")
            with open(file_path, "w") as f:
                json.dump(self.wsheet.data, f, indent=4, cls=CJsonEncoder)
            workflow.run()
        except CodeError, e:
            exstr = traceback.format_exc()
            print exstr
            if workflow:
                e.bind_object = workflow
                workflow.set_error(e)
            else:
                self.update_error(e)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            self.update_error(e)
        sys.exit()

    def update_error(self, e):
        if isinstance(e, CodeError):
            error = str(e)
            error_info = e.json()
        else:
            error = "运行异常:%s: %s" % (e.__class__.__name__, filter_error_info(e))
            error_info = {"error_type": "running", "code": "R001"}
        json_data = self.wsheet.data
        if json_data and "UPDATE_STATUS_API" in json_data.keys():
            json_obj = {"task": {
                "task_id": json_data["id"],
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "status": "failed",
                "run_time": 0,
                "progress": 0,
                "cpu_used": 0,
                "memory_used": 0
            },
                "log": [
                    {
                        "status": "failed",
                        "run_time": 0,
                        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "desc": "运行异常: %s" % e,
                        "error": error_info
                    }
                ]
            }
            post_data = {
                "sync_task_log": json_obj
            }
            data = {
                "task_id": json_data["id"],
                "api": json_data["UPDATE_STATUS_API"],
                "data": post_data
            }
            if "update_info" in self.wsheet.options().keys():
                data["update_info"] = self.wsheet.option('update_info')
            # log_client().add_log(data)
            # self.log_info_queue.put(("add_log", (data, )))
            add_log_queue(self.wsheet.log_info_queue, "add_log", data)
        # worker_client().set_error(json_data["id"], error)
        # self.run_info_queue.put(("set_error", (json_data["id"], error)))
        add_run_queue(self.wsheet.run_info_queue, "set_error", json_data["id"], error)
        self.logger.error("运行异常: %s " % error)


