# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from .core.event import EventObject
from .logger import Wlog
from .config import Config
# import time
# import random
import gevent
import json
# from .core.function import get_web_data
from .iofile import FileBase
import grpc
from .proto import workflow_guide_pb2_grpc, workflow_guide_pb2
import traceback
import sys
from .core.actor import LocalActor


class Batch(EventObject):
    def __init__(self, workflow, path, **kwargs):
        super(Batch, self).__init__()
        self._workflow = workflow
        # self._config = Config()
        self.config = Config()
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
        self._id = None
        self._has_run = False
        self.__init_events()
        self._submit_times = 0
        self.actor = LocalActor(self)
        self.version = 0

    @property
    def id(self):
        # if not self._id:
        #     self._id = "%s_%s_%s" % (self._workflow.id, time.strftime('%m%d%H%M%S', time.localtime(time.time())),
        #                              random.randint(1000, 10000))
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
        self.add_event("recivestate", loop=True)
        self.on("recivestate", self._call_state_callback)  # 收到state状态时调用对应的处理函数

    def _event_end(self):
        self._workflow.fire('childend', self)

    def _event_error(self, data):
        self._workflow.fire('childerror', self)
        if not self._ignore_error:
            self._workflow.exit(data=data)
        else:
            self._workflow.fire('childend', self)

    # def _create_json(self):
    #     data = {
    #         "id": self.id,
    #         "name": self._path,
    #         "type": self._batch_type,
    #         "batch": True,
    #         "batch_id": self._workflow.id,
    #         "UPDATE_STATUS_API": self._up_api_name,
    #         "IMPORT_REPORT_DATA": self._import_report_data,
    #         "output": self._upload_dir,
    #         "IMPORT_REPORT_AFTER_END": self._import_report_after_end,
    #         "client": self._config.WEB_INTERFACE_CLIENT,
    #         "options": self._options
    #     }
    #     return json.dumps(data)

    def run(self):
        self.start_listener()
        # try:
        #     json_data = self._create_json()
        # except:
        #     self._workflow.exit(data="生成投递配置失败！")
        #     return
        # self._post_pipeline(json_data)
        # gevent.sleep(3)
        dbversion = self.config.DBVersion if self.config.DBVersion else None
        workflow_data = workflow_guide_pb2.Batch(
            parent=self._workflow.id,
            name=self._path,
            type=self._batch_type,
            UPDATE_STATUS_API=self._up_api_name,
            IMPORT_REPORT_DATA=self._import_report_data,
            IMPORT_REPORT_AFTER_END=self._import_report_after_end,
            output=self._upload_dir,
            # cluster=self._workflow.sheet.CLUSTER,
            options=json.dumps(self._options),
            dbversion=dbversion,
        )
        self.actor.start()
        self._add_batch(workflow_data)

    def _add_batch(self, data):
        try:
            self._submit_times += 1
            with grpc.insecure_channel('localhost:%s' % self._workflow.sheet.wfm_port) as channel:
                stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
                response = stub.AddBatch(data)
                self._id = response.batch_id
                if response.skip:
                    self.logger.info("提交Batch成功,检测到Batch %s 已经完成,跳过运行" % response.batch_id)
                    self._work_dir = response.workdir
                    d = dict()
                    d["work_dir"] = response.workdir
                    self.end(d)
                elif response.ok:
                    self.logger.info("提交Batch成功,Batch ID : %s" % response.batch_id)
                else:
                    self.logger.error("提交Batch失败,原因 : %s" % response.reason)
                    self._workflow.exit(data="提交Batch失败")
            self._submit_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._submit_times > 3:
                self.logger.error("提交Batch发生错误超过3次,退出运行: %s" % e)
                self._submit_times = 0
                raise e
            else:
                self.logger.error("提交Batch失败运行发生错误, 30秒后重试: %s" % e)
                gevent.sleep(30)
                self._add_batch(data)

    # def _post_pipeline(self, json_data):
    #     self._submit_times += 1
    #     if self._workflow.sheet.CLUSTER:
    #         api = "pipeline?cluster=%s" % self._workflow.sheet.CLUSTER
    #     else:
    #         api = "pipeline"
    #     return_data = get_web_data(params={"json": json_data}, logger=self.logger, method="post", api=api)
    #     try:
    #         d = json.loads(return_data)
    #     except Exception, e:
    #         if self._submit_times > 3:
    #             self._workflow.exit(data="重试超过3次任然不成功，退出运行！")
    #         else:
    #             self.logger.error("接口返回内容不正确: %s, 30秒后重试!" % e)
    #             gevent.sleep(30)
    #             return self._post_pipeline(json_data)
    #     else:
    #         if d["success"]:
    #             self.logger.info("提交到远程API成功!")
    #             return True
    #         else:
    #             if self._submit_times > 3:
    #                 self._workflow.exit(data="重试超过3次任然不成功，退出运行！")
    #             else:
    #                 self.logger.error("提交失败: %s, 30秒后重试!" % d["info"])
    #                 gevent.sleep(30)
    #                 return self._post_pipeline(json_data)

    # def count_used(self):
    #     return self._cpu_used, self._mem_used

    def end(self, data):
        # if self.is_end:
        #     raise Exception("一个对象不能多次被结束，不能多次运行end方法!")
        # self.logger.info("检测到运行结束,运行节点: %s,  PID: %s , 工作目录: %s, CPU时: %s, 内存时: %s"
        #                  % (data["server"], data["pid"], data["work_dir"], data["cpu_used"], data["memory_used"]))
        # self._cpu_used = data["cpu_used"]
        # self._mem_used = data["memory_used"]
        self.logger.info("batch end")
        # self.work_dir = data["work_dir"]
        self.set_end()
        self.fire("end")
        self.actor.close()

    def error(self, data):
        # if self.is_end:
        #     raise Exception("一个对象不能多次被结束，不能多次运行error方法!")
        # self.logger.warning("检测到运行错误,运行节点: %s,  PID: %s , 工作目录: %s, CPU时: %s, 内存时: %s"
        #                     % (data["server"], data["pid"], data["work_dir"], data["cpu_used"], data["memory_used"]))
        # self._cpu_used = data["cpu_used"]
        # self._mem_used = data["memory_used"]
        # self.work_dir = data["work_dir"]
        self.logger.info("batch error")
        self.set_end()
        # self.actor.close()
        # self.fire("error", data["info"])
        if not self._ignore_error:
            self.fire("error", data["info"])
        else:
            self.fire("end")
        self.actor.close()

    def set_run(self):
        # if self.has_run:
        #     raise Exception("一个对象不能多次被开始，不能多次运行set_error方法!")
        self._has_run = True
        # self.logger.info("检测到开始运行,运行节点: %s,  PID: %s , 工作目录: %s"
        #                  % (data["server"], data["pid"], data["work_dir"]))
        self.fire("start")

    def _call_state_callback(self, message):
        """
        将自定义的callback函数转换为事件处理函数

        :param message:   接收到的消息
        :return:
        """
        # if "host" in message.keys():
        #     self._run_host = message['host']
        info = ""
        try:
            d = json.loads(message.data)
            if isinstance(d, dict) and "info" in d.keys():
                info = d["info"]
        except:
            self.logger.info("message.data:%s", message.data)
        data = {
            "process_id": message.jobid,
            "host": message.host,
            "info": info
        }
        # self.work_dir = d.work_dir
        self.logger.info("message:%s", message)
        self.logger.info("收到Batch 状态更新, ProcessID: %s, Host: %s " % (data["host"], data["process_id"]))
        if message.state == "runstart":
            self.set_run()
        elif message.state == "error":
            self.error(data)
        elif message.state == "end" or message.state == "finish":
            self.end(data)

    def _get_workflow(self, obj):
        # if obj.parent:
        #     return self._get_workflow(obj.parent)
        # else:
        #     return obj
        return obj
        # return self._workflow

    def get_workflow(self):
        """
        获取当前workflow对象

        :return:  :py:class:`biocluster.workflow.Workflow` 对象
        """
        return self._get_workflow(self)
