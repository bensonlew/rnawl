# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
"""workflow工作流类模块"""

from .core.function import load_class_by_path, CJsonEncoder
from .basic import Basic
from .config import Config
import os
import sys
from .rpc import RPC
from .logger import Wlog
# from .agent import Agent
# from .module import Module
import datetime
import gevent
import time
from biocluster.api.file.remote import RemoteFileManager
from biocluster.api.database.base import ApiManager
import re
import importlib
import types
import traceback
# from .core.watcher import Watcher
# from .scheduling.job import JobManager
# from .wpm.db import ReportModel
from .core.exceptions import RunningError, ExitError
from .batch import Batch
# import json
import copy
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import glob
# from boto.s3.bucket import Bucket
import grpc
from .proto import workflow_guide_pb2, workflow_guide_pb2_grpc, public_pb2
import json
import signal

PY3 = sys.version_info[0] == 3


class Workflow(Basic):
    """
    工作流程基类
    """

    def __init__(self, wsheet, **kwargs):
        if "debug" in kwargs.keys():
            self.debug = kwargs["debug"]
        elif wsheet.debug:
            self.debug = True
        else:
            self.debug = False

        # self._queue_to_stop_check = False

        super(Workflow, self).__init__(**kwargs)
        self.sheet = wsheet
        print(self)
        print("id {} name {}".format(self._id, self._name))

        self._return_msg = []  # 需要返回给任务调用进程的值,支持常用数据类型
        # self.last_update = datetime.datetime.now()
        if "parent" in kwargs.keys():
            self._parent = kwargs["parent"]
        else:
            self._parent = None

        self.pause = False
        self._pause_time = None
        # self.USE_DB = False
        self.__json_config()
        self._id = wsheet.id
        self.config = Config()
        self.config.current_mode = "workflow"
        os.environ['current_mode'] = "workflow"
        os.environ["NTM_PORT"] = str(wsheet.ntm_port)
        os.environ["WFM_PORT"] = str(wsheet.wfm_port)
        # 初始化配置
        # self.config.init_config(wsheet.PROJECT_TYPE, self._id, wsheet.wfm_port, wsheet.instant, wsheet.ntm_port)
        if wsheet.packagedir is not None:
            self.config.init_config_from_data(wsheet,self._id)
        else:
            self.config.init_config(wsheet.PROJECT_TYPE, self._id, wsheet.wfm_port, wsheet.instant, wsheet.ntm_port)

        if wsheet.DBVersion is not None:
            self.config.DBVersion = wsheet.DBVersion
        if wsheet.dydb == 1:
            self.config.ProjectID = wsheet.project_sn
        if not self.debug:
            self._work_dir = self.__work_dir()
            self._output_path = self._work_dir + "/output"
        # self._tools_report_data = []
        # self.add_event('rpctimeout', True)  # RPC接收信息超时
        # self.on("rpctimeout", self._event_rpctimeout)
        self.last_get_message = None
        if not self.debug:
            if not os.path.exists(self._work_dir):
                os.makedirs(self._work_dir)
            if not os.path.exists(self._output_path):
                os.makedirs(self._output_path)
            self.step_start()
            self.__check_to_file_option()
        self._logger = Wlog(self).get_logger("")
        # if self.sheet.instant is True:
        #     # self.USE_DB = False
        #     self.logger.debug("初始化本地RPC服务..!")
        #     # self.rpc_server = LocalServer(self)
        # else:
        #     self.logger.debug("初始化RPC服务器..!")
        self.rpc_server = RPC(self)
        signal.signal(41, self.rpc_server.get_tool_state)
        # self._stop_timeout_check = False
        self.is_skip = False
        self._start_time = None
        self.last_event_fire = None
        self.task_option_data = {}
        self._send_rpc_times = 0
        if wsheet.has_start_run:
            with wsheet.has_start_run.get_lock():
                wsheet.has_start_run.value = 1

    def __json_config(self):
        # if self.sheet.USE_DB is True:
        #     self.USE_DB = True
        if self.sheet.UPDATE_STATUS_API is not None:
            self.UPDATE_STATUS_API = self.sheet.UPDATE_STATUS_API
        if self.sheet.IMPORT_REPORT_DATA is True:
            self.IMPORT_REPORT_DATA = True
        if self.sheet.IMPORT_REPORT_AFTER_END is True:
            self.IMPORT_REPORT_AFTER_END = True

    # def _event_rpctimeout(self):
    #     """
    #     RPC信息接收时间超时时触发
    #
    #     :return:
    #     """
    #     if self.last_get_message:
    #         time_last = (datetime.datetime.now() - self.last_get_message).seconds
    #         if time_last < 30:
    #             self.logger.debug("RPC服务%s秒钟以前接收到信息，应该还在正常工作，忽略重启!" % time_last)
    #             return
    #     self.rpc_server.restart()

    # def add_tool_report(self, data):
    #     """
    #     添加Tool运行报告
    #     :return:
    #     """
    #     self._tools_report_data.append(data)

    def revise_infiles(self):
        pass
    def step_start(self):
        """
        流程开始api更新
        :return:
        """
        self.logger.info("开始更新步骤信息...")
        self.step.start()
        self.step.update()

    def __work_dir(self):
        """
        获取并创建工作目录
        """
        if self.sheet.work_dir:
            return self.sheet.work_dir
        if self.sheet.rerun and self.sheet.run_time:
            timestr = self.sheet.run_time.strftime('%Y%m%d')
        else:
            timestr = str(time.strftime('%Y%m%d', time.localtime(time.time())))
        if self.sheet.work_base:
            work_dir = self.sheet.work_base
        else:
            work_dir = self.config.WORK_DIR

        work_dir = work_dir + "/" + timestr + "/" + self.name + "_" + self._id
        return work_dir

    def __check_to_file_option(self):
        """
        转换内置的参数为文件参数

        :return:
        """
        if "to_file" in self._sheet.data.keys():
            data = self._sheet.data["to_file"]
            if PY3:
                check = isinstance(data, str)
            else:
                check = isinstance(data, types.StringTypes)
            if check:
                to_files = [data]
            else:
                to_files = data
            for opt in to_files:
                self.logger.debug("开始to_file参数文件转换%s..." % opt)
                m = re.match(r"([_\w.]+)\((.*)\)", opt)
                if m:
                    func_path = m.group(1)
                    f_list = func_path.split(".")
                    func_name = f_list.pop()
                    lib_path = ".".join(f_list)
                    options = m.group(2)
                    opt_list = re.split(r"\s*,\s*", options)
                    for optl in opt_list:
                        if re.match(r"^biocluster\.api\.to_file", lib_path):
                            imp = importlib.import_module("%s" % lib_path)
                        else:
                            imp = importlib.import_module("mbio.api.to_file.%s" % lib_path)
                            if hasattr(imp, "project_type"):
                                ptype = imp.project_type
                                if hasattr(imp, "db"):
                                    cli = self.config.get_mongo_client(mtype=ptype)
                                    imp.db = cli[self.config.get_mongo_dbname(ptype)]
                        func = getattr(imp, func_name)
                        self._sheet.set_option(optl, func(self._sheet.option(optl), optl, self.work_dir, self))
                else:
                    imp = importlib.import_module("biocluster.api.to_file.parameter")
                    func = getattr(imp, "json_to_file")
                    self._sheet.set_option(opt, func(self._sheet.option(opt), opt, self.work_dir, self))

                # json_data = self._sheet.option(opt)
                # file_path = os.path.join(self.work_dir, "%s_input.json" % opt)
                # with open(file_path, "w") as f:
                #     json.dump(json_data, f, indent=4)

    def add_module(self, path):
        """
        添加下属 :py:class:`biocluster.module.Module`

        :param path:  :py:class:`biocluster.module.Module` 对应的自动加载path路径，请参考教程中对应的说明

        """
        module = load_class_by_path(path, tp="Module")(self)
        self.add_child(module)
        return module

    def add_tool(self, path):
        """
        直接添加下属 :py:class:`biocluster.agent.Agent`

        :param path: String   :py:class:`biocluster.agent.Agent` 动态加载path路径
        :return:  Agent 返回Tool对应的 :py:class:`biocluster.agent.Agent` 对象
        """
        tool = load_class_by_path(path, tp="Agent")(self)
        self.add_child(tool)
        return tool

    def find_tool_by_id(self, toolid):
        """
        通过id搜索所属Tool

        :param toolid:  :py:class:`biocluster.tool.Tool` 对象的ID
        """
        # ids = toolid.split(".")
        # length = len(ids)
        # if length < 2 or length > 3:
        #     return False
        # if ids[0] != self.id:
        #     return False
        modules = copy.copy(self.children)
        # if length == 2:
        for md in modules:
            if str(md.id) == str(toolid):
                return md
            else:
                if hasattr(md, "find_tool_by_id"):
                    tool = md.find_tool_by_id(toolid)
                    if tool:
                        return tool
        return False

    def run(self):
        """
        开始运行

        :return:
        """
        super(Workflow, self).run()
        self.logger.info("worklow is start")
        # self._start_time = datetime.datetime.now()
        # watcher = Watcher()
        # watcher.add(self._check_batch, 10)
        # if self.sheet.instant is not True and self.sheet.WPM:
        #     watcher.add(self.__check, 2)
        # 停止
        signal.signal(42, self._stop)
        # 暂停
        signal.signal(43, self._pause)
        # 退出暂停
        signal.signal(44, self._exit_pause)
        self._update("start")  # 告诉WPM信号监听已经起来了
        self.rpc_server.run()

    def end(self):
        """
        停止workflow运行

        :return:
        """
        super(Workflow, self).end()
        try:
            self.stop_timeout_check()
            self._upload_result()
            self._import_report_data()
            # self._save_report_data()
            # self._update("end")
            self.step.finish()
            self.logger.info("运行step finish 成功!")
            self.step.update()
            self.logger.info("运行step update 成功!")
            self._update("end")
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            self.logger.error("结束异常流程: %s" % e)
            self._update("error", "结束异常流程: %s" % e)
        self.rpc_server.close()
        self.logger.info("workflow运行结束!")

    def set_return_msg(self, msg):
        """
        设置返回信息给WPM，使用WorkflowManager.get_msg方法可在运行完成后获取此信息，只能获取一次


        :param msg: 需要传递的信息，支持常用数据类型
        :return:
        """
        self._return_msg = msg

    def add_return_mongo_id(self, collection_name, table_id, desc='', add_in_sg_status=True):
        """
        为了兼容旧版本,新版本请使用set_return_msg

        :param collection_name:
        :param table_id:
        :param desc:
        :param add_in_sg_status:
        :return:
        """
        return_dict = dict()
        return_dict['id'] = table_id
        return_dict['collection_name'] = collection_name
        return_dict['desc'] = desc
        return_dict['add_in_sg_status'] = add_in_sg_status
        self._return_msg.append(return_dict)

    def _upload_result(self):
        """
        上传结果文件到远程路径

        :return:
        """
        if self._sheet.output:
            self.logger.info("开始上传结果文件")
            dirs = re.split(r"\s*;\s*", self._sheet.output)
            fm = MultiFileTransfer()
            count = 0
            for dir_path in dirs:
                for up in self.upload_dir:
                    # target_dir = os.path.join(self.sheet.output, os.path.dirname(up.upload_path))
                    target_dir = dir_path.rstrip("/") + '/'
                    from_dir = up.path
                    if not up.path.endswith("/"):
                        from_dir = up.path + "/"
                    fm.add_upload(from_dir, target_dir, up.path)
                    count += 1
                    # 去掉target目录加上了上传相对路径的路径差 即 上传目录 与工作目录的路径差
                    # remote_file = RemoteFileManager(target_dir)
                    # self.logger.info("开始上传%s到%s" % (up.path, target_dir))
                    # if remote_file.type != "local":
                    #     remote_file.upload(up.path)
                    #     self.logger.info("上传%s到%s完成" % (up.path, target_dir))
            if count > 0:
                fm.perform()

    def _import_report_data(self):
        if self.IMPORT_REPORT_DATA is True and self.IMPORT_REPORT_AFTER_END is True:
            api_call_list = self.api.get_call_records_list()
            if len(api_call_list) > 0:
                self.logger.info("开始导入Report数据...")
                api_manager = ApiManager(self, play_mod=True, debug=self.debug)
                api_manager.load_call_records_list(api_call_list)
                api_manager.play()
                self.logger.info("导入Report数据结束！")

    def exit(self, exitcode=1, data=None, terminated=False):
        """
        立即退出当前流程

        :param exitcode:
        :param data:
        :param terminated:
        :return:
        """
        # exstr = traceback.format_exc()
        # print exstr
        sys.stdout.flush()
        # self.end_unfinish_job()
        if isinstance(data, dict) and "error_type" in data.keys() and "info" in data.keys():
            error_str = json.dumps(data)
            if PY3:
                error_str = bytes(json.dumps(data), encoding='utf-8')
            # self._save_report_data(data)
        else:
            error_str = str(data)
            # self._save_report_data()
        if terminated:
            self.step.terminated(data)
        else:
            self.step.failed(data)
        self.step.update()
        self._update("error", error_str)  # 退出流程不再给wfm发送error状态
        self.logger.error("Failed: %s " % error_str)
        self.rpc_server.close()

        raise ExitError(error_str)

    # def __update_service(self):
    #     """
    #     每隔30s定时更新数据库last_update时间
    #
    #     :return:
    #     """
    #     # while self.is_end is False:
    #     #     gevent.sleep(60)
    #     if self.is_end is True:
    #         return "exit"
    #     self._update("online")
    #     # with self.db_sem:
    #     # try:
    #     #     with self.db.transaction():
    #     #         self.db.query("UPDATE workflow SET last_update=CURRENT_TIMESTAMP where workflow_id=$id",
    #     #                       vars={'id': self._id})
    #     # except Exception, e:
    #     #     exstr = traceback.format_exc()
    #     #     print exstr
    #     #     self.logger.debug("数据库更新异常: %s" % e)
    #     # finally:
    #     #     self.close_db_cursor()

    def _update(self, status, msg=None):
        """
        插入数据库，更新流程运行状态,只在后台服务调用时生效

        :param status: 更新信息类型
        :return:
        """
        # if self.USE_DB:
        #     try:
        #         myvar = dict(id=self._id)
        #         self.db.update("workflow", vars=myvar, where="workflow_id = $id", **data)
        #     except Exception, e:
        #         exstr = traceback.format_exc()
        #         print exstr
        #         self.logger.debug("数据库更新异常: %s" % e)
        #     finally:
        #         self.close_db_cursor()

        # if self.sheet.WPM:
        #     try:
        #         worker = worker_client()
        #         if type == "end":
        #             self.run_info_queue.put((type, (self.sheet.id, self._return_msg)))
        #         elif type == "keepalive":
        #             # worker.keep_alive(self.sheet.id)
        #             self.run_info_queue.put((type, (self.sheet.id, )))
        #         elif type == "error":
        #             # worker.set_error(self.sheet.id, msg)
        #             self.run_info_queue.put((type, (self.sheet.id,)))
        #         elif type == "pause":
        #             worker.set_pause(self.sheet.id)
        #         elif type == "stop":
        #             worker.set_stop(self.sheet.id)
        #         elif type == "pause_exit":
        #             worker.set_pause_exit(self.sheet.id)
        #         elif type == "pause_timeout":
        #             worker.pause_timeout(self.sheet.id)
        #         elif type == "report":
        #             worker.report(self.sheet.id, msg)
        #     except Exception, e:
        #         exstr = traceback.format_exc()
        #         print exstr
        #         sys.stdout.flush()
        #         self.logger.error("连接WPM服务异常: %s," % e)
        #     else:
        #         del worker
        #     if action == "set_end":
        #         add_run_queue(self.run_info_queue, "set_end", self.sheet.id, self._return_msg)
        #     elif action == "set_error":
        #         add_run_queue(self.run_info_queue, action, self.sheet.id, msg)
        #     else:
        #         add_run_queue(self.run_info_queue, action, self.sheet.id)
        if status == "end":
            msg = json.dumps(self._return_msg)
            if PY3:
                msg = bytes(json.dumps(self._return_msg), encoding='utf-8')
        statu = workflow_guide_pb2.Status(
            statu=status,
            workflow=public_pb2.Workflow(instant=self.sheet.instant, workflow_id=self.id, process_id=int(os.getpid())),
            message=msg,
        )
        try:
            pass

            # with grpc.insecure_channel('localhost:%s' % self.config.wfm_port) as channel:
            #     stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
            #     response = stub.Update(statu)
                # self.logger.info('response:%s' % response)
                # if response.ok:
                #     self.logger.info("发送status %s 成功" % status)
                # else:
                #     self.logger.debug("wfm服务拒绝接受status %s, 原因: %s" % (status, response.reason))
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            print(e)
            sys.stdout.flush()
            self.logger.error("发送status %s 失败" % status)

    # def __add_run_queue(self, action, data, msg=None):
    #     # length = len(self.run_info_queue)
    #     is_full = True
    #     with self.run_info_queue.get_lock():
    #         # if self.run_info_queue[length - 1].action:
    #         #     time.sleep(0.3)
    #         #     self.__add_run_queue(action, data, msg)
    #         #     return
    #         for i, v in enumerate(self.run_info_queue):
    #             if self.run_info_queue[i].action is None:
    #                 is_full = False
    #                 if msg is not None:
    #                     msg = json.dumps(msg, cls=CJsonEncoder)
    #                 self.run_info_queue[i] = (action, data, msg)
    #                 break
    #     if is_full:
    #         time.sleep(0.3)
    #         self.__add_run_queue(action, data, msg)
    #         return

    # def _add_log_queue(self, action, data, msg=None):
    #     # length1 = len(self.log_info_queue)
    #     is_full = True
    #     with self.log_info_queue.get_lock():
    #         # if self.log_info_queue[length1 - 1].action:
    #         #     time.sleep(0.3)
    #         #     self._add_log_queue(action, data, msg)
    #         #     return
    #         for i, v in enumerate(self.log_info_queue):
    #             if self.log_info_queue[i].action is None:
    #                 is_full = False
    #                 if msg is not None:
    #                     msg = json.dumps(msg, cls=CJsonEncoder)
    #                 self.log_info_queue[i] = (action, data, msg)
    #                 break
    #     if is_full:
    #         time.sleep(0.3)
    #         self._add_log_queue(action, data, msg)
    #         return

    # def send_log(self, data):
    #     """
    #     发送API LOG信息到WPM API LOG管理器
    #
    #     ;:param data: 需要发送的数据
    #     :return:
    #     """
    # if self.sheet.WPM:
            #
            # try:
            #     log = log_client()
            #     log.add_log(data)
            # except Exception, e:
            #     exstr = traceback.format_exc()
            #     print exstr
            #     sys.stdout.flush()
            #     self.logger.error("连接WPM log服务异常: %s" % e)
            # else:
            #     del log
            # self.log_info_queue.put(("add_log", (data, )))
            # json_str = json.dumps(data, cls=CJsonEncoder)
            # add_log_queue(self.log_info_queue, "add_log", data)

    # def __check(self):
    #     if self.is_end is True:
    #         return "exit"
    #     now = datetime.datetime.now()
    #     if self.pause:
    #         if (now - self._pause_time).seconds > self.config.MAX_PAUSE_TIME:
    #             self._update("pause_timeout")
    #             e = RunningError("暂停超过规定的时间%ss,自动退出运行!", (self.config.MAX_PAUSE_TIME, ), "002")
    #             e.bind_object = self
    #             self.exit(exitcode=1, data=e.json(), terminated=True)
    #
    #     if self.sheet.WPM:
    #         action = get_action_queue(self.action_queue, self.sheet.id)
    #         if action == "tostop":
    #             self._stop()
    #         elif action == "pause":
    #             self._pause()
    #         elif action == "exit_pause":
    #             self._exit_pause()
    #         if action is not None:
    #             self.logger.debug("接收到指令:%s" % action)
    #     manager = JobManager()
    #     if self.config.JOB_PLATFORM.lower() == "slurm" and len(manager.run_jobs) > 0:
    #         if manager.get_slurm_queue_job_number() > 0:  # 有slurm 任务在排队
    #             self._queue_to_stop_check = True
    #         else:
    #             if self._queue_to_stop_check is True:
    #                 # self.last_update = datetime.datetime.now()
    #                 self._queue_to_stop_check = False
    #     if len(manager.run_jobs) == 0 and len(manager.queue_jobs) > 0:  # 所有任务都在排队等待投递
    #         self._queue_to_stop_check = True
    #     else:
    #         if self._queue_to_stop_check is True:
    #             # self.last_update = datetime.datetime.now()
    #             self._queue_to_stop_check = False
    #     if not (self._stop_timeout_check is True or self._queue_to_stop_check is True or self.pause):
    #         # count = 0
    #         # for job in manager.run_jobs:
    #         #     if not job.is_end:
    #         #         count += 1
    #         # time_out = self.config.MAX_WAIT_TIME * 3 + 100
    #         # if count > 0 and (datetime.datetime.now() - self.last_update).seconds > time_out:
    #         #     self.set_error("超过 %s s没有任何运行更新，退出运行！", time_out, "003")
    #         if self.rpc_server.last_none_block \
    #                 and (datetime.datetime.now() - self.rpc_server.last_none_block).seconds < 180:
    #             if self.last_event_fire and \
    #                             (datetime.datetime.now() - self.last_event_fire).seconds > self.config.MAX_WAIT_TIME:
    #                 self.set_error("进程未被阻塞的情况下超过 %s s没有触发任何事件，退出运行！",
    #                                self.config.MAX_WAIT_TIME, "003")

    def stop_timeout_check(self):
        # self._stop_timeout_check = True
        self._send_rpc_times += 1
        workflow_data = public_pb2.Workflow(instant=self.sheet.instant,
                                            workflow_id=self.id,
                                            process_id=int(os.getpid()))
        try:
            with grpc.insecure_channel('localhost:%s' % self.config.wfm_port) as channel:
                stub = workflow_guide_pb2_grpc.WorkflowGuideStub(channel)
                response = stub.StopTimeoutCheck(workflow_data)
                if response.ok:
                    self.logger.info("发送stop_timeout_check成功")
                else:
                    self.logger.debug("ntm服务拒绝接受stop_timeout_check, 原因: %s" % response.reason)
                self._send_rpc_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._send_rpc_times > 5:
                self.logger.error("发送stop_timeout_check发生错误超过3次,退出运行: %s" % e)
                raise e
            else:
                self.logger.error("发送stop_timeout_check发生错误10秒后重试: %s" % e)
                gevent.sleep(10)
                self.stop_timeout_check()

    def _stop(self, signum, frame):
        self.logger.info("接收到终止运行指令。")
        self._update("stop")
        e = RunningError("接收到终止运行指令，退出运行!", None, "004")
        e.bind_object = self
        self.exit(exitcode=1, data=e.json(), terminated=True)

    def _pause(self, signum, frame):
        self.pause = True
        self._pause_time = datetime.datetime.now()
        self.step.pause()
        self.step.update()
        self._update("pause")
        self.logger.warning("检测到暂停指令，暂停所有新模块运行")

    def _exit_pause(self, signum, frame):
        self.pause = False
        self._pause_time = None
        self.step.start()
        self.step.update()
        self._update("continue")
        self.logger.info("检测到恢复运行指令，恢复所有模块运行!")

    # def _check_batch(self):
    #     if self.is_end is True:
    #         return "exit"
    #     has_batch_run = False
    #     if len(self._batch_list) > 0:
    #         for b in self._batch_list:
    #             if b.is_start and not b.is_end:
    #                 has_batch_run = True
    #     if not has_batch_run:
    #         return
    #     try:
    #         batch_list = json.loads(get_web_data({"batch_id": self.id}, self.logger, method="get"))
    #     except Exception, e:
    #         self.logger.warning("Web API返回内容格式错误: %s " % e)
    #     else:
    #         self.last_event_fire = datetime.datetime.now()
    #         for b in batch_list:
    #             batch = self.get_batch(b["workflow_id"])
    #             if batch:
    #                 if batch.is_end or not batch.is_start:
    #                     continue
    #                 if b["has_run"] == 1:
    #                     if not batch.has_run:
    #                         batch.set_run(b)
    #                     if b["is_error"] == 1:
    #                         batch.error(b)
    #                     else:
    #                         if b["is_end"] == 1:
    #                             batch.end(b)

    # def _save_report_data(self, error_data=None):
    #     try:
    #         if self.sheet.WPM:
    #             self.logger.debug("开始发送运行数据!")
    #             data = dict()
    #             data["modules"] = []
    #             for child in self.children:
    #                 if isinstance(child, Module):
    #                     data["modules"].append(child.get_report_data())
    #             data["tools"] = self._tools_report_data
    #             data["id"] = self.sheet.id
    #             (cpu, memory) = self.count_used()
    #             data["cpu_used"] = cpu
    #             data["memory_used"] = memory
    #             data["error_data"] = json.dumps(error_data)
    #             if self.sheet.rerun:
    #                 data["rerun"] = True
    #             file_path = os.path.join(self.work_dir, "run_report_data.json")
    #             with open(file_path, "w") as f:
    #                 json.dump(data, f, indent=4, cls=CJsonEncoder)
    #             try:
    #                 add_log_queue(self.log_info_queue, "report", data)
    #                 self.logger.debug("运行数据发送完成!")
    #             except MaxLengthError, e:
    #                 self.logger.debug("保存到Log共享内存出错，直接保存到Mysql: %s" % e)
    #                 self.__report_to_mysql(data)
    #             except Exception, e:
    #                 exstr = traceback.format_exc()
    #                 print exstr
    #                 sys.stdout.flush()
    #                 self.logger.error("保存运行报告异常: %s" % e)
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         sys.stdout.flush()
    #         self.logger.error("保存运行报告异常: %s" % e)

    # def __report_to_mysql(self, data):
    #     """
    #     保存运行数据
    #     :return:
    #     """
    #     wid = data["id"]
    #     self.logger.info("开始保存%s报告..." % wid)
    #     model = ReportModel(wid)
    #     try:
    #         if "cpu_used" in data.keys() and "memory_used" in data.keys():
    #             if self.sheet.rerun:
    #                 model.save_workflow(data["cpu_used"], data["memory_used"], data["error_data"], True)
    #             else:
    #                 model.save_workflow(data["cpu_used"], data["memory_used"], data["error_data"])
    #         if isinstance(data["modules"], list) and len(data["modules"]) > 0:
    #             model.save_modules(data["modules"])
    #         if isinstance(data["tools"], list) and len(data["tools"]) > 0:
    #             model.save_tools(wid, data["tools"], commit=True, under_workflow=1)
    #         self.logger.info("Workflow %s 保存运行报告完成" % wid)
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         sys.stdout.flush()
    #         self.logger.error("保存运行报告异常: %s" % e)
    #     finally:
    #         model.close()

    # def __check_tostop(self):
    #     """
    #     检查数据库的停止指令，如果收到则退出流程
    #
    #     :return:
    #     """
    #     # while self.is_end is False:
    #     # if self.is_end is True:
    #     #     return "exit"
    #     # if (datetime.datetime.now() - self.last_update).seconds > self.config.MAX_WAIT_TIME:
    #     #     self.exit(data="超过 %s s没有任何运行更新，退出运行！" % self.config.MAX_WAIT_TIME)
    #     # gevent.sleep(10)
    #     # myvar = dict(id=self._id)
    #     # try:
    #     #     results = self.db.query("SELECT * FROM tostop "
    #     #                             "WHERE workflow_id=$id and done  = 0", vars={'id': self._id})
    #     #     if isinstance(results, long) or isinstance(results, int):
    #     #         self.close_db_cursor()
    #     #         gevent.sleep(10)
    #     #         return
    #     #     if len(results) > 0:
    #     #         data = results[0]
    #     #         update_data = {
    #     #             "stoptime": datetime.datetime.now(),
    #     #             "done": 1
    #     #         }
    #     #         self.db.update("tostop", vars=myvar, where="workflow_id = $id", **update_data)
    #     #         self.exit(data="接收到终止运行指令,%s" % data.reson, terminated=True)
    #     #
    #     # except Exception, e:
    #     #     exstr = traceback.format_exc()
    #     #     print exstr
    #     #     self.logger.info("查询数据库异常: %s" % e)
    #     # finally:
    #     #     self.close_db_cursor()
    #
    # def __check_pause(self):
    #     """
    #     检查暂停指令或终止暂停指令
    #
    #     :return:
    #     """
    #     if self.is_end is True:
    #         return "exit"
    #     myvar = dict(id=self._id)
    #     try:
    #         results = self.db.query("SELECT * FROM pause WHERE workflow_id=$id and "
    #                                 "has_continue  = 0 and timeout = 0", vars={'id': self._id})
    #         if isinstance(results, long) or isinstance(results, int):
    #             self.close_db_cursor()
    #             gevent.sleep(10)
    #             return
    #         if len(results) > 0:
    #             data = results[0]
    #             if data.has_pause == 0:
    #                 self.pause = True
    #                 self._pause_time = datetime.datetime.now()
    #                 update_data = {
    #                     "pause_time": datetime.datetime.now(),
    #                     "has_pause": 1
    #                 }
    #                 self.db.update("pause", vars=myvar, where="workflow_id = $id", **update_data)
    #                 self.db.query("UPDATE workflow SET paused = 1 where workflow_id=$id", vars={'id': self._id})
    #                 self.step.pause()
    #                 self.step.update()
    #                 self.logger.info("检测到暂停指令，暂停所有新模块运行: %s" % data.reason)
    #             else:
    #                 if data.exit_pause == 0:
    #                     now = datetime.datetime.now()
    #                     if self.pause:
    #                         if (now - self._pause_time).seconds > self.config.MAX_PAUSE_TIME:
    #                             update_data = {
    #                                 "timeout_time": datetime.datetime.now(),
    #                                 "timeout": 1
    #                             }
    #                             self.db.update("pause", vars=myvar, where="workflow_id = $id", **update_data)
    #                             self.db.query("UPDATE workflow SET paused = 0 where workflow_id=$id",
    #                                           vars={'id': self._id})
    #                             self.exit(data="暂停超过规定的时间%ss,自动退出运行!" %
    #                                            self.config.MAX_PAUSE_TIME, terminated=True)
    #                 else:
    #                     if data.has_continue == 0 and data.timeout == 0:
    #                         self.pause = False
    #                         self._pause_time = None
    #                         update_data = {
    #                                 "continue_time": datetime.datetime.now(),
    #                                 "has_continue": 1
    #                         }
    #                         self.db.update("pause", vars=myvar, where="workflow_id = $id", **update_data)
    #                         self.db.query("UPDATE workflow SET paused = 0 where workflow_id=$id",
    #                                       vars={'id': self._id})
    #                         self.step.start()
    #                         self.step.update()
    #                         self.logger.info("检测到恢复运行指令，恢复所有模块运行!")
    #     except Exception, e:
    #         exstr = traceback.format_exc()
    #         print exstr
    #         self.logger.info("查询数据库异常: %s" % e)
    #     finally:
    #         self.close_db_cursor()

    # def end_unfinish_job(self):
    #     """
    #     结束所有未完成的job任务
    #
    #     :return:
    #     """
    #     manager = JobManager()
    #     for job in manager.get_unfinish_jobs():
    #         job.delete()
    #         job.agent.canceled_report()
    #
    #     if hasattr(self, "process_share_manager"):
    #         self.process_share_manager.shutdown()

    # def close_db_cursor(self):
    #     cursor = self.db._db_cursor()
    #     cursor.close()

    # def count_used(self):
    #     childs = self.children
    #     cpu_used = 0
    #     memory_used = 0
    #     for c in childs:
    #         (x, y) = c.count_used()
    #         cpu_used += x
    #         memory_used += y
    #     for b in self._batch_list:
    #         (x, y) = b.count_used()
    #         try:
    #             cpu_used += x
    #             memory_used += y
    #         except:
    #             pass
    #     return cpu_used, memory_used

    def add_batch(self, path, ignore_error=False, batch_type="workflow", upload_dir=None, up_api_name=None,
                  import_report_data=True, import_report_after_end=True):
        if batch_type not in ["workflow", "module", "tool"]:
            raise Exception("batch_type类型错误!")
        batch = Batch(self, path, ignore_error=ignore_error, batch_type=batch_type,
                      upload_dir=upload_dir, up_api_name=up_api_name,
                      import_report_data=import_report_data, import_report_after_end=import_report_after_end)
        # if self.get_batch(batch.id):
        #     return self.add_batch(path, ignore_error=ignore_error, batch_type=batch_type,
        #                           upload_dir=upload_dir, up_api_name=up_api_name,
        #                           import_report_data=import_report_data,
        #                           import_report_after_end=import_report_after_end)
        # else:
        self._batch_list.append(batch)
        return batch

    def get_batch(self, batch_id):
        for b in self._batch_list:
            if batch_id == b.id:
                return b

    # def download_from_s3(self, from_file, to_path="download/", cover=True):
    #     """
    #     从s3对象存储下载数据到本地, 为了避免堵塞进程，此功能应该放置在流程最后执行。
    #     :param from_file: 需要下载的文件路径或文件路径, 必须是类似s3region://bucket/key写法。
    #     因为对象存储中没有文件夹的概念，需要下载文件夹必须使用"/"结尾，以明确表明下载的是文件夹
    #     :param to_path: 下载文件相对于当前工作目录的存放目录。
    #     当路径为"/"结尾时，表示下载文件存放在此文件夹下，否者为下载完整路径。
    #     当from_file为文件夹时，此参数也必须以"/"结尾。目录层级与下载的s3目录层级结构相同。
    #     默认情况下放置在当前模块工作目录的download目录下。
    #     :param cover: 对已存在的文件是否覆盖
    #     :return:
    #     """
    #     if re.match(r"^/|^\.\.", to_path):
    #         raise Exception("不能使用绝对路径或切换到其他目录!")
    #     if os.path.basename(to_path) == ".":
    #         raise Exception("目标文件不能叫\".\"!")
    #     target_dir = False
    #     if re.match(r"/$", to_path):
    #         target_dir = True
    #     fm =
    #     # self.s3transfer = TransferManager()
    #     # self.s3transfer.base_path = self.work_dir
    #     # self.s3transfer.overwrite = cover
    #     # m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
    #     # if not m:
    #     #     raise Exception("下载路径%s格式不正确!" % from_file)
    #     # else:
    #     #     region = m.group(1)
    #     #     bucket_name = m.group(2)
    #     #     key_name = m.group(3)
    #     #     if re.match(r"/$", key_name):
    #     #         if not target_dir:
    #     #             raise Exception("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
    #     #         conn = self.s3transfer.config.get_rgw_conn(region, bucket_name)
    #     #         bucket = Bucket(connection=conn, name=bucket_name)
    #     #         for key in bucket.list(prefix=key_name):
    #     #             source = os.path.join(from_file, key.name)
    #     #             target = os.path.join(target_dir, os.path.relpath(key.name, key_name))
    #     #             self.s3transfer.add(source, target)
    #     #     else:
    #     #         if not target_dir:  # 处理已存在文件的情况
    #     #             target = os.path.join(self.work_dir, to_path)
    #     #             if os.path.exists(target):
    #     #                 if cover:
    #     #                     if os.path.isdir(target):
    #     #                         shutil.rmtree(target)
    #     #                     else:
    #     #                         os.remove(target)
    #     #                 else:
    #     #                     raise Exception("目标文件夹%s已经存在!" % target)
    #     #         else:
    #     #             target = os.path.join(self.work_dir, to_path, os.path.basename(key_name))
    #     #         self.s3transfer.add(from_file, target)
    #     # self.s3transfer.wait()

    def upload_to_s3(self, from_file, to_path, cover=True):
        """
        从本地上传数据到S3对象存储, 为了避免堵塞进程，此功能应该放置在流程最后执行。
        :param from_file: 需要上传的文件，相对于当前模块工作目录的相对路径，可以是文件或文件夹，也可使用通配符。
        当为文件夹时会上传文件夹的目录结构。通配符匹配到的文件夹也会上传目录结构。
        :param to_path: 上传文件的路径，必须是类似s3region://bucket/key写法。
        当路径为"/"结尾时，表示上传文件存放在此文件夹下，否者为上传完整路径。
        当使用通配符时必须为"/"结尾，且通配符匹配到的文件名不能重复，否则会相互覆盖。
        :param cover: 对已存在的文件是否覆盖
        :return:
        """
        if not re.match(r"^\w+://\S+/.+$", to_path):
            raise Exception("上传路径%s格式不正确!" % to_path)
        if os.path.basename(to_path) == ".":
            raise Exception("目标文件不能叫\".\"!")
        fm = MultiFileTransfer()
        fm.add_upload(from_file, to_path, from_file)
        fm.perform()
        # target_dir = False
        # if re.match(r"/$", to_path):
        #     target_dir = True
    #     self.s3transfer = TransferManager()
    #     self.s3transfer.base_path = self.work_dir
    #     self.s3transfer.overwrite = cover
    #     source = os.path.join(self.work_dir, from_file)
    #     for f in glob.glob(os.path.join(source)):
    #         if target_dir:
    #             target = os.path.join(to_path, os.path.basename(f))
    #         else:
    #             target = to_path
    #         if os.path.isdir(f):
    #             for root, dirs, files in os.walk(f):
    #                 rel_path = os.path.relpath(root, f)
    #                 for i_file in files:
    #                     if rel_path == ".":
    #                         key_path = os.path.join(target, i_file)
    #                     else:
    #                         key_path = os.path.join(target, rel_path, i_file)
    #                     i_file_path = os.path.join(root, i_file)
    #                     self.s3transfer.add(i_file_path, key_path)
    #         else:
    #             self.s3transfer.add(f, target)
    #     self.s3transfer.wait()

    @property
    def project_type(self):
        """
        获取和设置所属项目类型名称，用于控制database api查找相应的数据库配置
        """
        if hasattr(self.config, "PROJECT_TYPE"):
            return self.config.PROJECT_TYPE
        else:
            return None

    @project_type.setter
    def project_type(self, value):
        self.config.PROJECT_TYPE = value

    def save_error_manual_post_data(self):
        if len(self.upload_dir) > 0 and self.sheet.output:
            spend_time = (datetime.datetime.now() - self._start_time).seconds
            json_obj = {
                "task": {
                    "task_id": self.sheet.id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    # "status": "run",
                    "run_time": spend_time,
                    "progress": 100
                },
                "log": [{
                    "status": "finish",
                    "run_time": spend_time,
                    "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "desc": "运行完成",
                }]
            }
            files, dirs = self.step.upload_files()
            m = re.match(r"^([\w\-]+)://(.*)$", self.sheet.output)
            if m:
                json_obj['region'] = m.group(1)
                json_obj['base_path'] = m.group(2)
            else:
                m1 = re.match(r"^\w+:(.*)$", self.sheet.output)
                if m1:
                    json_obj['base_path'] = m1.group(1)
                else:
                    json_obj['base_path'] = self.sheet.output
            json_obj['files'] = files
            json_obj['dirs'] = dirs
            # (cpu, memory) = self.count_used()
            # json_obj['task']['cpu_used'] = cpu
            # json_obj['task']['memory_used'] = memory
            post_data = {
                'sync_task_log': json_obj
            }
            # for k, v in self.step._api_data:
            #     post_data[k] = v
            wpm_log_data = {
                "task_id": self.sheet.id,
                'api': self.step.api_type,
                "data": post_data,
                "process_id": os.getpid()
            }
            # if "update_info" in self.sheet.options().keys():
            #     wpm_log_data["update_info"] = self.sheet.option('update_info')
            file_path = os.path.join(self.work_dir, "error_manual_post_data.json")
            with open(file_path, "w") as f:
                json.dump(wpm_log_data, f, indent=4, cls=CJsonEncoder)
        else:
            return

    def add_task_option(self, name, value):
        if name in ["status", "memory_used", "created_ts", "task_id", "run_time", "progress", "cpu_used"]:
            raise Exception("名称不能为%s！" % name)
        if name in self.task_option_data.keys():
            raise Exception("名称%s已经存在，不能重复添加！" % name)
        self.task_option_data[name] = value
