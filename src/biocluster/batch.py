# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from .core.event import EventObject
from .logger import Wlog
from .config import Config
import time
import random
import gevent
import json
from .core.function import get_web_data
from .iofile import FileBase


class Batch(EventObject):
    def __init__(self, workflow, path, **kwargs):
        super(Batch, self).__init__()
        self._workflow = workflow
        self._config = Config()
        self._options = {}
        self._path = path
        self._upload_dir = kwargs["upload_dir"] if "upload_dir" in kwargs.keys() else None
        self._batch_type = kwargs["batch_type"] if "batch_type" in kwargs.keys() else "workflow"
        self._up_api_name = kwargs["up_api_name"] if "up_api_name" in kwargs.keys() else None
        self._ignore_error = kwargs["ignore_error"] if "ignore_error" in kwargs.keys() else None
        self._import_report_data = kwargs["import_report_data"] if "import_report_data" in kwargs.keys() else None
        self._import_report_after_end = kwargs["import_report_after_end"] if "import_report_after_end" \
                                                                             in kwargs.keys() else None
        self._logger = None
        self._work_dir = ""
        self._cpu_used = 0
        self._mem_used = 0
        self._id = None
        self._has_run = False
        self.__init_events()
        self._submit_times = 0

    @property
    def id(self):
        if not self._id:
            self._id = "%s_%s_%s" % (self._workflow.id, time.strftime('%m%d%H%M%S', time.localtime(time.time())),
                                     random.randint(1000, 10000))
        return self._id

    @property
    def has_run(self):
        return self._has_run

    @property
    def work_dir(self):
        """
        返回当前对象的工作目录,start事件触发后可用
        :return:
        """
        return self._work_dir

    @work_dir.setter
    def work_dir(self, value):
        self._work_dir = value

    @property
    def logger(self):
        if self._logger:
            return self._logger
        else:
            self._logger = Wlog(self._workflow).get_logger("Batch:" + self._path + "(" + self.id + ")")
            return self._logger

    def set_options(self, options):
        if not isinstance(options, dict):
            raise Exception("参数格式错误!")
        for name, value in options.items():
            if isinstance(value, FileBase):
                self._options[name] = value.path
            else:
                self._options[name] = value

    def __init_events(self):
        self.add_event("start")
        self.add_event("end")
        self.on("end", self._event_end)
        self.add_event("error")
        self.on("error", self._event_error)

    def _event_end(self):
        self._workflow.fire('childend', self)

    def _event_error(self, data):
        self._workflow.fire('childerror', self)
        if not self._ignore_error:
            self._workflow.exit(data=data)
        else:
            self._workflow.fire('childend', self)

    def _create_json(self):
        data = {
            "id": self.id,
            "name": self._path,
            "type": self._batch_type,
            "batch": True,
            "batch_id": self._workflow.id,
            "UPDATE_STATUS_API": self._up_api_name,
            "IMPORT_REPORT_DATA": self._import_report_data,
            "output": self._upload_dir,
            "IMPORT_REPORT_AFTER_END": self._import_report_after_end,
            "client": self._config.WEB_INTERFACE_CLIENT,
            "options": self._options
        }
        return json.dumps(data)

    def run(self):
        self.start_listener()
        try:
            json_data = self._create_json()
        except:
            self._workflow.exit(data="生成投递配置失败！")
            return
        self._post_pipeline(json_data)
        gevent.sleep(3)

    def _post_pipeline(self, json_data):
        self._submit_times += 1
        return_data = get_web_data({"json": json_data}, self.logger, method="post", api="pipeline")
        try:
            d = json.loads(return_data)
        except Exception, e:
            if self._submit_times > 3:
                self._workflow.exit(data="重试超过3次任然不成功，退出运行！")
            else:
                self.logger.error("接口返回内容不正确: %s, 30秒后重试!" % e)
                gevent.sleep(30)
                return self._post_pipeline(json_data)
        else:
            if d["success"]:
                self.logger.info("提交到远程API成功!")
                return True
            else:
                if self._submit_times > 3:
                    self._workflow.exit(data="重试超过3次任然不成功，退出运行！")
                else:
                    self.logger.error("提交失败: %s, 30秒后重试!" % d["info"])
                    gevent.sleep(30)
                    return self._post_pipeline(json_data)

    def count_used(self):
        return self._cpu_used, self._mem_used

    def end(self, data):
        if self.is_end:
            raise Exception("一个对象不能多次被结束，不能多次运行end方法!")
        self.logger.info("检测到运行结束,运行节点: %s,  PID: %s , 工作目录: %s, CPU时: %s, 内存时: %s"
                         % (data["server"], data["pid"], data["work_dir"], data["cpu_used"], data["memory_used"]))
        self._cpu_used = data["cpu_used"]
        self._mem_used = data["memory_used"]
        self.work_dir = data["work_dir"]
        self.set_end()
        self.fire("end")

    def error(self, data):
        if self.is_end:
            raise Exception("一个对象不能多次被结束，不能多次运行error方法!")
        self.logger.warning("检测到运行错误,运行节点: %s,  PID: %s , 工作目录: %s, CPU时: %s, 内存时: %s"
                            % (data["server"], data["pid"], data["work_dir"], data["cpu_used"], data["memory_used"]))
        self._cpu_used = data["cpu_used"]
        self._mem_used = data["memory_used"]
        self.work_dir = data["work_dir"]
        self.set_end()
        self.fire("error", data["error_data"])

    def set_run(self, data):
        if self.has_run:
            raise Exception("一个对象不能多次被开始，不能多次运行set_error方法!")
        self._has_run = True
        self.work_dir = data["work_dir"]
        self.logger.info("检测到开始运行,运行节点: %s,  PID: %s , 工作目录: %s"
                         % (data["server"], data["pid"], data["work_dir"]))
        self.fire("start")
