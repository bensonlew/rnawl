# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""Tool远程代理"""

from .basic import Basic
from .core.actor import LocalActor, RemoteData
import pickle
import os
from .config import Config
# from .scheduling.job import JobManager
from .core.function import get_classpath_by_object, load_class_by_path
# import datetime
import inspect
import gevent
import re
import types
from .iofile import FileBase
import traceback
import sys
import grpc
from .proto import workflow_guide_pb2_grpc, workflow_guide_pb2
import json

PY3 = sys.version_info[0] == 3


class PickleConfig(object):
    """
    保存配置信息，用于发送到远程 :py:class:`biocluster.tool.Tool` 对象
    """
    def __init__(self):
        self._name = ""
        self._full_name = ""
        self._id = ""
        self._work_dir = ""
        # self.endpoint = ""
        self._options = ""
        self._output_path = ""
        self.version = ""
        self.SOFTWARE_DIR = ""
        self.PACKAGE_DIR = ""
        # self.KEEP_ALIVE_TIME = ""
        # self.MAX_KEEP_ALIVE_TIME = ""
        self.PROJECT_TYPE = None
        # self.RGW_ENABLE = None
        # self._parent_id = ""
        # self._parent_path = ""
        # self._parent_work_dir = ""
        # self._parent_start_time = ""
        # self._grandparent_id = ""
        self.JOB_PLATFORM = None
        self.SCRIPT_DIR = None
        self.wfm_port = 7321
        self.ntm_port = 7322
        self.instant = True
        self.current_workflow_id = ""
        self.current_tool_id = ""
        self.DBVersion = None
        self._memory_increase_step = 10

    def clone(self, agent):
        """
        从agent的属性中克隆同名属性值到自身

        :param agent:   :py:class:`Agent` 对象
        """
        for name in vars(self).keys():
            if hasattr(agent, name):
                setattr(self, name, getattr(agent, name))
            else:
                if hasattr(agent.config, name):
                    setattr(self, name, getattr(agent.config, name))

    def save(self, file_handler):
        """
        保存自身到文件对象

        :param file_handler:   文件句柄
        """
        pickle.dump(self, file_handler)


class Agent(Basic):
    """
    提供 :py:class:`biocluster.tool.Tool` 远程代理，通过Agent-Tool对，使远程Tool信息得到相应的处理
    """
    def __init__(self, parent):
        super(Agent, self).__init__(parent)
        self.config = Config()
        # self.config.init_current_tool_id(self._id)
        self.actor = LocalActor(self)
        self._queue = self.config.JOB_QUEUE
        self.current_tool_id = self._id
        self._host = ""
        self._cpu = 0
        self._memory = "1G"
        # self._default_callback_action = {'action': 'none'}
        # self._callback_action = {}
        self._status = "W"                # W 等待 Q 排队 R 运行 E 完成
        # self.add_event('keepaliveout', loop=True)   # 保持连接超时
        # self.on('keepaliveout', self._event_keepaliveout)
        # self.add_event('waittimeout', loop=True)    # 等待超时
        # self.on('waittimeout', self._event_waittimeout)
        # self.add_event('firekaoout', loop=True)  # keepaliveout超过最大次数
        # self.on('firekaoout', self._event_firekaoout)
        # self.add_event('firewtoout')  # waittimeout超过最大次数
        # self.on('firewtoout', self._event_firewtoout)
        # self.add_event('runstart')       # 远端开始运行
        # self.on("runstart", self._event_runstart)
        # self.add_event('runstartout', loop=True)       # 远端开始运行等待超时
        # self.on("runstartout", self._event_runstart_out)
        self.on("error", self._agent_event_error)
        self.add_event("recivestate", loop=True)
        self.on("recivestate", self._call_state_callback)  # 收到state状态时调用对应的处理函数
        self.add_event("sync", loop=True)
        # self.endpoint = self.get_workflow().rpc_server.endpoint
        # self.job = None
        self._run_mode = "Auto"      # 运行模式 Auto表示由 main.conf中的 platform参数决定
        # self._job_manager = JobManager()
        # self._run_time = None
        # self._start_run_time = datetime.datetime.now()
        # self._end_run_time = None
        self._rerun_time = 0
        # self.is_wait = False
        self._remote_data = RemoteData(self)
        self.version = 0  # 重投递运行次数
        # self.shared_callback_action = {}  # 进程共享变量，
        self._run_host = "unknown"
        # self._report_data = None
        # self.cpu_used = 0
        # self.memory_used = 0
        # self._is_skip = False
        self._max_rerun_time = 3
        self._memory_increase_step = 10
        # self._parent_id = self.parent.id
        # self._parent_path = get_classpath_by_object(self.parent)
        # self._parent_work_dir = self.parent.work_dir
        # self._parent_start_time = getattr(self.parent, "_start_time", datetime.datetime.now())
        # if self.parent.parent:
        #     self._grandparent_id = self.parent.parent.id
        # else:
        #     self._grandparent_id = ""
        # self._not_found_times = 0
        # self._load_report_error = 0
        self._grpc_times = 0
        # self.instant = self.get_workflow().sheet.instant
        self.instant = True
        self._submit_times = 0

    def _call_state_callback(self, message):
        """
        将自定义的callback函数转换为事件处理函数

        :param message:   接收到的消息
        :return:
        """
        self.logger.info("执行recivestate的callback函数")
        self.logger.info("message.state：{}".format(message.state))
        # if "host" in message.keys():
        #     self._run_host = message['host']
        if message.state == "runstart":
            self.logger.info("远程任务开始运行，任务ID:%s, 任务类型: %s,远程主机:%s" %
                             (message.jobid, message.jobtype, message.host))
            return
        if hasattr(self, message.state + '_callback'):
            func = getattr(self, message.state + '_callback')
            args = []
            if PY3:
                spec = inspect.signature(func)
                for i in spec.parameters.keys():
                    args.append(i)
            else:
                spec = inspect.getargspec(func)
                args = spec.args
                if inspect.ismethod(func):
                    args.pop(0)
            if len(args) == 0:
                func()
            elif len(args) == 1:
                try:
                    data = json.loads(message.data)
                    self.logger.info("data:%s" % data)
                    self.logger.info("type data:%s" % type(data))
                    func(data)
                except Exception:
                    self.logger.info("message.data:%s" % message.data)
                    self.logger.info("type message.data:%s" % type(message.data))
                    func(message.data)

            else:
                raise Exception("状态回调函数参数不能超过2个(包括self)!")
        else:
            self._default_callback(message)

    def _default_callback(self, message):
        """
        消息处理函数不存在时对默认的处理方法

        :param message:   接收到的消息
        :return:
        """
        self.logger.warning(self.name + "没有定义消息对应的处理函数" + message.state + "!")
        

    @property
    def queue(self):
        """
        获取队列名
        """
        return self._queue

    @queue.setter
    def queue(self, queue):

        self._queue = queue

    @property
    def cpu(self):
        return int(self._cpu)

    @cpu.setter
    def cpu(self, cpu):
        self._cpu = cpu

    @property
    def memory(self):
        """
        返回Job申请的内存大小，单位为字节数
        :return:
        """
        if PY3:
            check = isinstance(self._memory, str)
        else:
            check = isinstance(self._memory, types.StringTypes)

        if not check:
            raise Exception("%s的内存设置只能为字符串，不能为其他类型!")
        pattern = re.compile(r'([\d.]+)([gmk])', re.I)
        match = re.match(pattern, self._memory)
        if match:
            unit = match.group(2)
            if unit.upper() == "G":
                return float(match.group(1)) * 1024 * 1024 * 1024
            elif unit.upper() == "M":
                return float(match.group(1)) * 1024 * 1024
            elif unit.upper() == "K":
                return float(match.group(1)) * 1024
        else:
            pattern = re.compile(r'([\d.]+)', re.I)
            match = re.match(pattern, self._memory)
            if match:
                return float(match.group(1))
            else:
                if self._memory == "":
                    self._memory = "1G"
                else:
                    raise Exception("%s的内存设置必须为数字加单位G/M/K!" % self.fullname)
        return 0

    @memory.setter
    def memory(self, memory):
        self._memory = memory

    @property
    def mode(self):
        """
        返回运行模式，默认"Auto"  表示由 main.conf中的 platform参数决定，可以在子类中重写_run_mode属性来定义运行模式
        """
        return self._run_mode

    @property
    def remote(self):
        """
        获取远程数据对象

        :return:
        """
        return self._remote_data

    @property
    def status(self):
        return self._status

    # @property
    # def report_data(self):
    #     """
    #     获取运行报告, 只在运行结束后有效
    #     :return:
    #     """
    #     return self._report_data

    # def add_remote_data(self, name, value):
    #     """
    #     添加需要传递到远程Tool的数据
    #
    #     :param name:
    #     :param value:
    #     :return:
    #     """
    #     if name in self._remote_data.keys():
    #         raise Exception("远程数据名称%s已经存在，请勿重复添加" % name)
    #     if not isinstance(name, types.StringType):
    #         raise Exception("远程数据名称必须为字符串")
    #     elif not name.islower():
    #         raise Exception("命令名称必须都为小写字母！")
    #     if not (isinstance(value, types.StringTypes) or isinstance(value, types.BooleanType) or
    #             isinstance(value, types.IntType) or isinstance(value, types.LongType) or
    #             isinstance(value, types.FloatType) or isinstance(value, types.TupleType) or
    #             isinstance(value, types.ListType) or isinstance(value, types.DictType)):
    #         raise Exception("远程数据值必须为Python内置数据类型: 字符串，数字，布尔，list,tuple,dict！")
    #     self._remote_data[name] = value

    def set_queue(self, queue):
        """
        设置队列名

        :param queue:
        :return: None
        """
        self._queue = queue

    def get_resource(self):
        """
        获取所需资源，必须通过在子类中重写 :py:func:`set_resource` ，并定义 self._cpu 和 self._memory 属性

        :return: cpu,memory
        """
        if self._cpu < 1:
            raise Exception("必须重写方法set_resource,并指定所需资源")
        return self._cpu, self._memory

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory

        :return:
        """
        self._cpu = 0
        self._memory = "1GB"

    def save_config(self):
        """
        保存远程 :py:class:`biocluster.tool.Tool` 运行所需参数

        :return: 文件路径
        """
        path = os.path.join(self.work_dir, self.name + ".pk")
        for option in self._options.values():
            option.bind_obj = None
            if isinstance(option.value, FileBase):
                option.value.option = None
        # with open(path, "w") as f:
        with open(path, "wb") as f:
            config = PickleConfig()
            config.clone(self)
            config.save(f)
        for option in self._options.values():
            option.bind_obj = self
        return path

    def save_class_path(self):
        """
        保存加载运行远程 :py:class:`biocluster.tool.Tool` 所需的类包路径

        :return: path  文件路径
        """
        path = os.path.join(self.work_dir, self.name + "_class.pk")
        class_name = str(type(self))
        tool_path = get_classpath_by_object(self)
        class_list = {"tool": get_classpath_by_object(self)}   # 类文件路径
        file_class_paths = []  #
        for option in self._options.values():
            if option.type in {'outfile', 'infile'}:
                if option.format:
                    file_class_paths.append(option.format)
                else:
                    for f in option.format_list:
                        file_class_paths.append(f)
        class_list['files'] = file_class_paths
        # with open(path, "w") as f:
        with open(path, "wb") as f:
            pickle.dump(class_list, f)
        return path

    def load_output(self):
        """
        从远端 :py:class:`biocluster.tool.Tool` 保存的pk文件中读取Option输出值，并赋值给自身Option参数

        :return:
        """
        output_path = os.path.join(self.work_dir, self.name + "_output.pk")
        if os.path.exists(output_path):
            try:
                # with open(output_path, "r") as f:
                with open(output_path, "rb") as f:
                    output = pickle.load(f)
                for option in self._options.values():
                    if option.type in {'outfile', 'infile'}:
                        if option.format:
                            load_class_by_path(option.format, "File")
                        else:
                            for f in option.format_list:
                                load_class_by_path(f, "File")
                for name, value in output.items():
                    self.option(name, value)
            except Exception as e:
                exstr = traceback.format_exc()
                print(exstr)
                print(e)
                sys.stdout.flush()

    # def _load_report(self, success=1, info=""):
    #     path = os.path.join(self.work_dir, self.name + "_run_report.pk")
    #     if not os.path.exists(path):
    #         return
    #     with open(path, "r") as f:
    #         try:
    #             output = pickle.load(f)
    #         except EOFError:
    #             self._load_report_error += 1
    #             if self._load_report_error > 3:
    #                 return
    #             gevent.sleep(2)
    #             return self._load_report(success, info)
    #
    #     self._report_data = {
    #         "run_id": self.id,
    #         "path": get_classpath_by_object(self),
    #         "work_dir": self.work_dir,
    #         "run_host": self._run_host,
    #         "job_type": self.mode if self.mode != "Auto" else self.config.JOB_PLATFORM,
    #         "job_id": self.job.id,
    #         "request_cpu": self._cpu,
    #         "request_memory": self._memory,
    #         "start_time": int(time.mktime(self._run_time.timetuple())),
    #         "wait_spend_time": (self.job.submit_time - self._run_time).seconds,
    #         "queue_spend_time": (self._start_run_time - self.job.submit_time).seconds,
    #         "run_spend_time": (self._end_run_time - self._start_run_time).seconds,
    #         "end_time": int(time.mktime(self._end_run_time.timetuple())),
    #         "run_times": self._rerun_time + 1,
    #         "success": success,
    #         "info": info,
    #         "process_id": output["process_id"],
    #         "sub_process_num": output["sub_process_num"],
    #         "max_cpu_use": output["max_cpu_use"],
    #         "max_rss": output["max_rss"],
    #         "average_cpu_use": output["average_cpu_use"],
    #         "average_rss": output["average_rss"],
    #         "max_vms": output["max_vms"],
    #         "average_vms": output["average_vms"],
    #         "commands": output["cmd"],
    #         "params": self.__get_params_data()
    #     }
    #     self.parent.add_tool_report(self._report_data)
    #     if output["max_cpu_use"]:
    #         max_use = float(output["max_cpu_use"])
    #         if max_use < 1:
    #             max_use = 100
    #         max_use = float(max_use) / 100
    #         if self.cpu > max_use:
    #             cpu_max = self.cpu
    #         else:
    #             cpu_max = max_use
    #     else:
    #         cpu_max = self.cpu
    #     run_time = (self._end_run_time - self._start_run_time).seconds
    #     if run_time < 1:
    #         run_time = 1
    #     self.cpu_used = float(cpu_max) * (float(run_time)/3600)
    #     if output["max_vms"]:
    #         max_vms = long(output["max_vms"])
    #     else:
    #         max_vms = self.memory
    #     if self.memory > max_vms:
    #         max_vms = self.memory
    #     if max_vms > 0:
    #         max_vms = float(max_vms) / (1024 * 1024 * 1024)
    #     else:
    #         max_vms = 0.1
    #     self.memory_used = float(max_vms) * (float(run_time)/3600)
    #     path = os.path.join(self.work_dir, self.name + ".resource_used.pk")
    #     with open(path, "w") as f:
    #         pickle.dump({"cpu": self.cpu_used, "memory": self.memory_used}, f)

    # def canceled_report(self):
    #     path = os.path.join(self.work_dir, self.name + "_run_report.pk")
    #     if os.path.exists(path):
    #         with open(path, "r") as f:
    #             try:
    #                 output = pickle.load(f)
    #             except EOFError:
    #                 output = None
    #     else:
    #         output = None
    #     self._end_run_time = datetime.datetime.now()
    #     self._report_data = {
    #         "run_id": self.id,
    #         "path": get_classpath_by_object(self),
    #         "work_dir": self.work_dir,
    #         "run_host": self._run_host,
    #         "job_type": self.mode if self.mode != "Auto" else self.config.JOB_PLATFORM,
    #         "job_id": self.job.id,
    #         "request_cpu": self._cpu,
    #         "request_memory": self._memory,
    #         "start_time": int(time.mktime(self._run_time.timetuple())),
    #         "wait_spend_time": (self.job.submit_time - self._run_time).seconds,
    #         "queue_spend_time": (self._start_run_time - self.job.submit_time).seconds,
    #         "run_spend_time": (self._end_run_time - self._start_run_time).seconds,
    #         "end_time": int(time.mktime(self._end_run_time.timetuple())),
    #         "run_times": self._rerun_time + 1,
    #         "success": 0,
    #         "info": "流程错误被取消!",
    #         "process_id": output["process_id"] if output else 0,
    #         "sub_process_num": output["sub_process_num"] if output else 0,
    #         "max_cpu_use": output["max_cpu_use"] if output else 0,
    #         "max_rss": output["max_rss"] if output else 0,
    #         "average_cpu_use": output["average_cpu_use"] if output else 0,
    #         "average_rss": output["average_rss"] if output else 0,
    #         "max_vms": output["max_vms"] if output else 0,
    #         "average_vms": output["average_vms"] if output else 0,
    #         "commands": output["cmd"] if output else 0,
    #         "params": self.__get_params_data()
    #     }
    #     self.parent.add_tool_report(self._report_data)

    def run(self):
        """
        开始运行，并投递任务

        :return:
        """
        super(Agent, self).run()
        self.set_resource()
        # self._run_time = datetime.datetime.now()
        # workflow = self.get_workflow()
        # if workflow.sheet.rerun and \
        #         (workflow.is_skip or workflow.sheet.skip_all_success or
        #          (workflow.sheet.skip_tools and self.id in workflow.sheet.skip_tools)):
        #     self.__set_rerun_dir()
        #     is_sucess = False
        #     path = os.path.join(self.work_dir, self.name + ".status.pk")
        #     if os.path.exists(path):
        #         is_sucess = True
        #     if is_sucess:
        #         self._is_skip = True
        #         self.fire('runstart', "pass")
        #         gevent.spawn_later(3, self.finish_callback)
        #     else:
        #         self.__run()
        # else:
        self.__run()

    def __run(self):
        self.step.start()
        self.remote.save()
        # if self.get_workflow().sheet.instant:
        #     self._run_mode = "process"
        self.save_class_path()
        self.save_config()

        # self.job = self._job_manager.add_job(self)
        # self._job_manager.run_job(self.job)
        self._status = "Q"
        # if not self.get_workflow().sheet.instant:
        self.actor.start()

        self._submit_job()
        gevent.sleep(0)

    # def __set_rerun_dir(self):
    #     dir_path = self._work_dir
    #     path = os.path.join(self.work_dir, self.name + ".status.pk")
    #     if os.path.exists(path):
    #         return
    #     for i in range(1, self._max_rerun_time+1):
    #         new_path = "%s__%s" % (dir_path, i)
    #         if os.path.exists(new_path):
    #             self._work_dir = new_path
    #             self._output_path = new_path + "/output"
    #             path = os.path.join(new_path, self.name + ".status.pk")
    #             if os.path.exists(path):
    #                 return

    def rerun(self):
        """
        删除远程任务并重新运行

        :return:
        """
        self._rerun_time += 1
        if self._rerun_time > self._max_rerun_time:
            self.fire("error", "重运行超过%s次仍未成功!" % self._max_rerun_time)
        else:
            self.stop_listener()
            self.restart_listener()
            self.version += 1
            # match = re.match(r"^(.*)__\d$", self._work_dir)
            # if match:
            #     self._work_dir = "%s__%s" % (match.group(1), self.version)
            # else:
            #     self._work_dir = "%s__%s" % (self._work_dir, self.version)
            self._output_path = self._work_dir + "/output"
            self.create_work_dir()
            # if not self.get_workflow().sheet.instant:
            self.actor.kill()
            self.save_class_path()
            self.save_config()
            # self._run_time = datetime.datetime.now()
            self._status = "Q"
            self.actor.auto_break = True
            # self.actor.kill()
            self.logger.info("开始重新投递任务!")
            self.send_exit_action(reason="重运行", version=self.version)
            self.actor = LocalActor(self)
            self.actor.start()
            self._submit_job()

    def _submit_job(self):
        cpu, mem = self.get_resource()
        if mem == "":
            mem = "1G"
        if not (("g" in mem) or ("G" in mem)):
            mem = "1G"
        if int(mem.rstrip("G")) > 1000:
            mem = "1000G"
        if int(mem.rstrip("G")) < 1:
            mem = "1G"
        if int(cpu) > 20:
            cpu = 20
        py_version = "py2"
        if PY3:
            py_version = "py3"
        class_name = str(type(self))
        class_path = get_classpath_by_object(self)
        class_path_split = class_path.split("mbio.tools.")
        if len(class_path_split) > 1:
            path = class_path_split[1]
        else:
            path = class_path
        workflow = self.get_workflow()
        data = workflow_guide_pb2.Job(
            instant=self.instant,
            workflow_id=workflow.id,
            process_id=int(os.getpid()),
            tool_id=self.id,
            # job_name=self.fullname,
            job_name=path,
            work_dir=self.work_dir,
            cpu=int(cpu),
            memory=mem,
            queue=self.queue,
            version=self.version,
            toolname=self.name,
            py_version=py_version,
        )
        try:
            self._submit_times += 1
            with grpc.insecure_channel('localhost:%s' % workflow.sheet.wfm_port) as channel:
                stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
                response = stub.Submit(data)
                if response.skip:
                    self.logger.info("提交Job成功,接收到信息跳过Job运行...")
                    # self.finish_callback()
                    self.actor.skip = True
                elif response.ok:
                    self.logger.info("提交Job成功")
                else:
                    self.logger.error("提交Job失败,原因: %s" % response.reason)
                    workflow.exit(data="提交Job失败,原因: %s" % response.reason)
            self._submit_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._submit_times > 5:
                self.logger.error("提交Job发生错误超过3次,退出运行: %s" % e)
                raise e
            else:
                self.logger.error("提交Job失败发生错误30秒后重试: %s" % e)
                gevent.sleep(30)
                self._submit_job()

    # def _set_callback_action(self, action, data=None, version=None):
    #     """
    #     设置需要发送给远程的Action指令,只会被远程获取一次
    #
    #     :param action: string
    #     :param data: python内置数据类型,需要被传递给远程的数据
    #     :param version: agent版本编号
    #     :return:
    #     """
    #     if version is None:
    #         version = self.version
    #     if self.get_workflow().sheet.instant:
    #         self._callback_action = self.shared_callback_action
    #     self._callback_action["%s" % version] = {"id": self.id, 'action': action, 'data': data, "version": version}

    # def get_callback_action(self, version):
    #     """
    #     获取需要返回给远程的Action信息,此信息只会被获取一次
    #
    #     :param version: 远程tool版本号
    #     :return: None
    #     """
    #     key = "%s" % version
    #     if key in self._callback_action.keys():
    #         action = self._callback_action[key]
    #         del self._callback_action[key]
    #     else:
    #         action = self._default_callback_action
    #     return action

    def send_exit_action(self, reason="", version=None):
        """
        发送exit指令到远程tool,此指令将导致远程Tool自动退出

        :param reason: 发送exit指令的原因
        :param version: agent版本编号
        :return:
        """
        if version is None:
            version = self.version
        self.logger.info("发送GRPC请求结束远程Tool运行 版本%s ...." % version)

        self._grpc_times += 1
        workflow = self.get_workflow()
        try:
            tool = workflow_guide_pb2.Tool(tool_id=self.id, workflow_id=workflow.id, process_id=int(os.getpid()),
                                           reason=reason, version=version)
            with grpc.insecure_channel('localhost:%s' % workflow.sheet.wfm_port) as channel:
                stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
                response = stub.Delete(tool)
                if response.ok:
                    self.logger.info("GRPC删除Tool运行成功, 版本%s" % version)
                else:
                    self.logger.info("GRPC删除Tool运行失败, 版本%s ,原因: %s" % (version, response.reason))
                self._grpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._grpc_times > 5:
                self.logger.error("GRPC删除Tool运行发生错误超过3次,退出运行: %s" % e)
                self._grpc_times = 0
                raise e
            else:
                self.logger.error("GRPC删除Tool运行发生错误, 10秒后重试: %s" % e)
                gevent.sleep(20)
                self.send_exit_action(reason, version)

    # def send_rerun_action(self, reason, version=None):
    #     """
    #     发送rerun指令到远程tool,此指令讲导致远程Tool重新运行。
    #     一般情况下，需要先修改agent自身的option参数，然后再发送此指令，修改后的参数值将在远程Tool中生效。
    #     警告：未经修改option而直接发送此指令，将导致远程Tool终止并重新运行，并不会产生其他的效果。
    #
    #     :param reason: 发送rerun指令的原因
    #     :param version: agent版本编号
    #     :return:
    #     """
    #     if version is None:
    #         version = self.version
    #     self.logger.info("发送rerun指令到远程Tool,版本%s ...." % version)
    #     self.save_class_path()
    #     self.save_config()
    #     # self._set_callback_action("rerun", data=reason, version=version)

    def finish_callback(self):
        """
        收到远程发送回的 :py:class:`biocluster.core.actor.State` end状态时的处理函数，设置当前Agent状态为结束

        :return:
        """
        self.load_output()
        self.remote.load()
        self._status = "E"

        # if self._is_skip:
        #     # self._load_resource_used_file()
        #     self.end()
        # else:
        # self.job.set_end()
        # self._end_run_time = datetime.datetime.now()
        # secends = (self._end_run_time - self._start_run_time).seconds if self._start_run_time else 0
        # self._load_report()
        self.actor.close()
        self.end()
        # self._write_finish_file()
        self.logger.info("任务运行结束")

    # def _write_finish_file(self):
    #     """
    #     将正常完成的状态写入文件，用于重运行
    #     :return:
    #     """
    #     path = os.path.join(self.work_dir, self.name + ".status.pk")
    #     if not os.path.exists(path):
    #         with open(path, "w") as f:
    #             pickle.dump(True, f)

    # def _load_resource_used_file(self):
    #     path = os.path.join(self.work_dir, self.name + ".resource_used.pk")
    #     if os.path.exists(path):
    #         with open(path, "r") as f:
    #             try:
    #                 output = pickle.load(f)
    #                 self.memory_used = output["memory"]
    #                 self.cpu_used = output["cpu"]
    #             except EOFError:
    #                 gevent.sleep(2)
    #                 self._load_resource_used_file()
    #     else:
    #         self.memory_used = 0
    #         self.cpu_used = 0

    def error_callback(self, data):
        """
        收到远程发送回的 :py:class:`biocluster.core.actor.State` error 错误状态时的处理函数, 默认触发error事件，并发送Action指令命令远程Tool退出::

            如不希望发现错误时退出，可以重写此方法

        :param data: 远程state信息中的data信息
        :return: None
        """
        # self._end_run_time = datetime.datetime.now()
        # self.job.set_end()
        # if isinstance(data, dict) and "info" in data.keys():
        #     self._load_report(success=0, info=data["info"])
        # else:
        #     self._load_report(success=0, info=data)
        self.actor.close()
        self._status = "E"
        self.fire("error", data)

    # def memory_limit_callback(self, data=""):
    #     """
    #     接收到Slurm 任务调度系统因内存超标而终止运行时发送的memory_limit状态时执行的状态处理函数
    #
    #     :return:
    #     """
    #     self._increase_memory(self._memory_increase_step)
    #     self.logger.warning("发现Tool出现exceeded virtual memory limit! ,调整内存申请大小到%s重新运行!\n%s" %
    #                         (self._memory, data))
    #     self.parent.fire("childrerun", self)

    # def killed_callback(self, data=""):
    #     self.logger.warning("接收到Killed状态:%s, 30秒后开始检测远程任务输出..." % data)
    #     gevent.sleep(30)
    #     if self.job.memory_limit():
    #         self.memory_limit_callback()
    #     elif self.job.cancelled():
    #         self.cancelled_callback()
    #     else:
    #         self.logger.warning("远程任务原因未知:%s,尝试重运行..." % data)
    #         if self.parent:
    #             self.parent.fire("childrerun", self)

    def sync_callback(self, data):
        """
        接收到同步数据请求时加载最新的远程数据

        :return:
        """
        self.logger.debug("接收到远程Tool同步数据请求,tag:%s，加载最新数据.." % data)
        self.remote.load()
        self.fire("sync", data)

    # def cancelled_callback(self, data=""):
    #     self.logger.warning("发现Tool远程任务被取消:%s,尝试重新投递新任务!" % data)
    #     if self.parent:
    #         self.parent.fire("childrerun", self)

    def _increase_memory(self, size=10):
        """
        增加任务申请内存量

        :param size: 增加大小 单位G

        :return:
        """
        pattern = re.compile(r'([\d.]+)([gmk])', re.I)
        match = re.match(pattern, self._memory)
        if match:
            unit = match.group(2)
            if unit.upper() == "G":
                self._memory = "%sG" % (int(match.group(1)) + size)
            elif unit.upper() == "M":
                self._memory = "%sM" % (int(match.group(1)) + size * 1024)
            elif unit.upper() == "K":
                self._memory = "%sK" % (int(match.group(1)) + size * 1024 * 1024)
        else:
            pattern = re.compile(r'([\d.]+)', re.I)
            match = re.match(pattern, self._memory)
            if match:
                self._memory = "%s" % (int(match.group(1)) + size*1024*1024*1024)

    # def sleeping_callback(self, data):
    #     """
    #
    #     :return:
    #     """
    #     self.logger.warning("Command %s长时间处于Sleep/Zombie状态，尝试重新运行!" % data)

    def _agent_event_error(self, data):
        """
        Agent发生错误时默认的处理方式

        :return:
        """
        # self._end_run_time = datetime.datetime.now()
        # secends = (self._end_run_time - self._start_run_time).seconds if self._start_run_time else 0
        if isinstance(data, dict) and "info" in data.keys():
            self.logger.error("发现运行错误:%s, " % data["info"])
        else:
            self.logger.error("发现运行错误:%s" % data)
        # self.job.set_end()
        self.actor.close()
        # self._job_manager.remove_all_jobs()
        self.get_workflow().exit(data=data)

    # def _event_keepaliveout(self):
    #     """
    #     当远程 :py:class:`biocluster.tool.Tool` 通信间隔时间超过规定时，执行此方法 在main.conf max_keep_alive_time中配置此时间
    #
    #     默认情况下将会删除远程Tool所属的任务，并重新投递
    #
    #     :return:
    #     """
    #     workflow = self.get_workflow()
    #     if workflow.rpc_server.check_rpc(self.id):
    #         return
    #     if self.is_end:
    #         return
    #     self.job.check_state()
    #     if self.job.is_error():
    #         self.logger.debug("远程任务运行出错，准备重新运行!")
    #         self.job.delete()
    #         if self.parent:
    #             if self.job.memory_limit():
    #                 self.memory_limit_callback()
    #             else:
    #                 if not workflow.rpc_server.check_rpc(self.id):
    #                     self.parent.fire("childrerun", self)
    #     elif self.job.is_running():
    #         self.logger.debug("远程任务正在运行，但是未能更新状态，尝试重启RPC服务!")
    #         self.get_workflow().fire("rpctimeout")
    #     elif self.job.is_completed():
    #         self.logger.debug("远程任务已结束，但是未正常更新状态，一分钟后重新检测!")
    #         gevent.sleep(60)
    #         if self.is_end:
    #             return
    #         if not self.job.is_end:
    #             if not workflow.rpc_server.check_rpc(self.id):
    #                 self.logger.debug("远程任务已结束，但是仍然未能更新状态，尝试重新运行,并重启RPC服务,并重新运行!")
    #                 self.get_workflow().rpc_server.restart()
    #         self.job.delete()
    #         if self.parent:
    #             if self.job.memory_limit():
    #                 self.memory_limit_callback()
    #             else:
    #                 if not workflow.rpc_server.check_rpc(self.id):
    #                     if self.is_end:
    #                         return
    #                     self.parent.fire("childrerun", self)
    #     else:
    #         if self.job.state == "Not_found" and self._not_found_times <= 5:
    #             self._not_found_times += 1
    #             self.logger.debug("远程任务状态未知: %s ，可能是已经运行结束，检测RPC信息!" % self.job.state)
    #             workflow.rpc_server.check_rpc(self.id)
    #         elif self.parent and (not workflow.rpc_server.check_rpc(self.id)):
    #             self.logger.debug("远程任务状态未知: %s ，尝试重新运行!" % self.job.state)
    #             self.parent.fire("childrerun", self)

    # def _event_waittimeout(self):
    #     """
    #
    #     当远程 :py:class:`biocluster.tool.Tool` 超过规定时间未能启动通信时执行此方法，在main.conf max_wait_time中配置此时间
    #
    #     默认情况下将会删除远程Tool所属的任务，并重新投递
    #
    #     :return:
    #     """
    #     workflow = self.get_workflow()
    #     if workflow.rpc_server.check_rpc(self.id):
    #         return
    #     if self.is_end:
    #         return
    #     self.job.check_state()
    #     if self.job.is_queue():
    #         self.logger.error("当前任务正在排队，继续等待...")
    #     elif self.job.is_error():
    #         self.logger.error("远程任务运行出错，准备重新运行!")
    #         self.job.delete()
    #         if self.parent:
    #             if self.job.memory_limit():
    #                 self.memory_limit_callback()
    #             else:
    #                 if not workflow.rpc_server.check_rpc(self.id):
    #                     self.parent.fire("childrerun", self)
    #     elif self.job.is_running():
    #         self.logger.error("远程任务已开始运行，等待1分钟后重新开始查询状态!")
    #         gevent.sleep(60)
    #         if self.is_end:
    #             return
    #         if not self._start_run_time:
    #             self.logger.error("远程任务已开始运行，但是未正常更新状态，尝试重新运行!")
    #             self.job.delete()
    #             if self.parent:
    #                 if self.job.memory_limit():
    #                     self.memory_limit_callback()
    #                 else:
    #                     if not workflow.rpc_server.check_rpc(self.id):
    #                         if self.is_end:
    #                             return
    #                         self.parent.fire("childrerun", self)
    #     elif self.job.is_completed():
    #         self.logger.error("远程任务已结束，但是未正常更新状态，尝试重新运行!")
    #         if self.parent and (not workflow.rpc_server.check_rpc(self.id)):
    #             self.parent.fire("childrerun", self)
    #     else:
    #         if self.job.state == "Not_found" and self._not_found_times <= 5:
    #             self._not_found_times += 1
    #             self.logger.debug("远程任务状态未知: %s ，可能是已经运行结束，检测RPC信息!" % self.job.state)
    #             workflow.rpc_server.check_rpc(self.id)
    #         elif self.parent and (not workflow.rpc_server.check_rpc(self.id)):
    #             self.logger.error("远程任务状态未知: %s ，尝试重新运行!" % self.job.state)
    #             self.parent.fire("childrerun", self)

    # def _event_runstart_out(self, times):
    #     """
    #     等待运行超时
    #
    #     :param times:
    #     :return:
    #     """
    #     if times > 5:
    #         self.logger.warning("等待运行超时第%s此触发，开始检查任务状态!" % times)
    #         self._event_waittimeout()
    #     else:
    #         self.logger.warning("等待运行超时第%s此触发，忽略且继续等待..." % times)

    # def _event_firekaoout(self, times):
    #     if self.job.is_running():
    #         self.logger.warning("KeepAlive触发超过%s次，但是远程任务仍在运行，继续等待...!" % times)
    #         self.actor._kao_fire_times = 0
    #         self.actor._has_firekaoout = False
    #     else:
    #         self.logger.warning("KeepAlive触发超过%s次，尝试重新运行!" % times)
    #         if self.parent:
    #             self.parent.fire("childrerun", self)

    # def _event_firewtoout(self, times):
    #     self.logger.warning("WaitTimeOut触发超过%s次，尝试重新运行!" % times)
    #     if self.parent:
    #         self.parent.fire("childrerun", self)

    # def count_used(self):
    #     return self.cpu_used, self.memory_used
