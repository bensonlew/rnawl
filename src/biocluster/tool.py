# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from .core.actor import RemoteActor, RemoteData
from .core.actor import State
import types
from .command import Command
from .logger import Wlog
import sys
import pickle
import os
import threading
import time
import re
import importlib
from .core.function import stop_thread, friendly_size
from .api.database.base import ApiManager
import signal
import inspect
import psutil
import traceback
from .package import PackageManager, Package
import shutil
# from biocluster.api.file.lib.s3 import S3TransferManager
import glob
from .core.exceptions import OptionError, RunningError, CodeError
from boto.s3.bucket import Bucket
from .config import Config
from .iofile import FileBase


class Tool(object):
    """
    远程运行工具，与Agent对应
    """

    def __init__(self, config):
        """
        初始化并加载config

        :param config: PickleConfig对象
        :return:
        """
        super(Tool, self).__init__()
        self.config = config
        self._name = ""
        self._full_name = ""
        self._id = ""
        self._work_dir = ""
        self._commands = {}
        self._states = []
        self._output_path = ""
        self._run = False
        self._end = False
        self._options = {}
        self.version = 0
        self.instant = False  # 本地进程模式
        self.load_config()  
        self.logger = Wlog(self).get_logger('')
        self.main_thread = threading.current_thread()
        self.mutex = threading.Lock()
        self.exit_signal = False
        self._rerun = False
        self.process_queue = None
        self.shared_callback_action = None
        # self._has_record_commands = []
        self._remote_data = RemoteData(self)
        self._remote_data.load()
        self.is_wait = False
        self.receive_exit_signal = False
        self.api = ApiManager(self)
        if self.instant is True:
            self.actor = None
        else:
            self.actor = RemoteActor(self, self.main_thread)
            if self.config.DEBUG is not True:
                self.actor.start()
                self.logger.debug("启动Actor线程!")

        self._cpu_max_use = 0
        self._cup_avg_use = 0
        self._mem_max_rss = 0
        self._mem_avg_rss = 0
        self._mem_max_vms = 0
        self._mem_avg_vms = 0
        self._set_resource_count = 0
        self._all_processes = []
        self.package = PackageManager(self)
        self.s3transfer = None

    @property
    def remote(self):
        """
        设置远程agent数据
        :return:
        """
        return self._remote_data

    @property
    def name(self):
        """
        工具名称，与其Agent名称一致
        :return:
        """
        return self._name

    @property
    def is_end(self):
        """
        是否运行完成
        :return:
        """
        return self._end

    @property
    def id(self):
        """
        id，与其Agent id一致
        :return:
        """
        return self._id

    @property
    def states(self):
        """
        返回发送到远程Agent的状态列表,list State对象列表
        :return:
        """
        return self._states

    @property
    def work_dir(self):
        """
        工作目录

        :return:
        """
        return self._work_dir

    @property
    def output_dir(self):
        return self._output_path

    @property
    def commands(self):
        return self._commands

    def add_command(self, name, cmd, script_dir=False, default_return_code=0, ignore_error=False, shell=False):
        """
        执行命令生成Command对象，
        可以在自类中定义命令检测函数.，函数名为 commandname_check ,此类函数将在self.run调用后自动加入微线程运行,3秒钟执行一次
        commandname_check传入参数为名称为commandname的Command对象

        :param name:  命令名称
        :param cmd: 需要执行的命令，路径必须为相对于配置文件main.conf的software_dir的相对路径
        :param script_dir: 是否是脚本目录，脚本程序的路径由main.conf中的[Command] script_dir指定
        :param default_return_code 正常结束时的返回值
        :param ignore_error  非正常结束时是否忽略不报错
        :param shell  是否是shell，shell命令不会添加路径前缀，也可以使用管道符和重定向符
        :return: 返回Command对象
        """
        if name in self._commands.keys():
            raise Exception("命令名称已经存在，请勿重复添加")
        if not isinstance(name, types.StringType):
            raise Exception("命令名称必须为字符串")
        elif not name.islower():
            raise Exception("命令名称必须都为小写字母！")
        else:
            cmd = Command(name, cmd, self, default_return_code, ignore_error, shell)
            if script_dir:
                cmd.software_dir = self.config.SCRIPT_DIR.rstrip("/")
            self._commands[name] = cmd
            return cmd

    def get_option_object(self, name=None):
        """
        通过参数名获取当前对象的参数 :py:class:`biocluster.option.Option` 对象

        :param name: string 参数名，可以不输入，当不输出时返回当前对象的所有参数对象 list输出
        :return: :py:class:`biocluster.option.Option` object 或 :py:class:`biocluster.option.Option` object数组
        """
        if not name:
            return self._options
        elif name not in self._options.keys():
            raise Exception("参数%s不存在，请先添加参数" % name)
        else:
            return self._options[name]

    def option(self, name, value=None):
        """
        获取/设置对象的参数值

        :param name: 参数名称
        :param value: 当value==None时，获取参数值 当value!=None是，设置对应的参数值
        :return: 参数对应的值
        """
        if name not in self._options.keys():
            raise Exception("参数%s不存在，请先添加参数" % name)
        if value is None:
            return self._options[name].value
        else:
            self._options[name].value = value

    def set_options(self, options):
        """
        批量设置参数值

        :param options: dict key是参数名,value是参数值
        :return:
        """
        if not isinstance(options, dict):
            raise Exception("参数格式错误!")
        for name, value in options.items():
            self.option(name, value)

    def wait(self, *cmd_name_or_objects):
        """
        暂停当前线程并等待,指定的Command对象或Package对象执行完成,如果cmdname未指定，则等待所有Command运行完成

        :cmdname:  一个或多个cmd名称 或Commnad对象
        """
        cmds = []
        if len(cmd_name_or_objects) > 0:
            for c in cmd_name_or_objects:
                if isinstance(c, Command) or isinstance(c, Package):
                    cmds.append(c)
                else:
                    if c not in self._commands.keys():
                        raise Exception("Commnad名称不存在！")
                    cmds.append(self._commands[c])
        else:
            cmds.extend(self._commands.values())
            cmds.extend(self.package.packages)
        while True:
            time.sleep(1)
            if self.exit_signal:
                self.logger.info("接收到退出信号，终止程序运行!")
                break
            if len(cmds) == 0:
                break
            is_running = False
            for command in cmds:
                if command.has_run:
                    if command.is_end is False:
                        is_running = True
                    else:
                        command.join_thread()
                        if command.return_code != command.default_return_code:
                            if command.ignore_error is False:
                                self.set_error("命令 %s 返回码: %s , 运行未能正常完成!",
                                               (command.name, command.return_code), "006")
                else:
                    is_running = True
            if is_running is False:
                break

    def add_state(self, name, data=None):
        """
        添加状态,用于发送远程状态，

        :param name: string State名称
        :param data: 需要传递给Agent相关的数据,data必须为python内置的简单数据类型
        :return: self
        """
        if not isinstance(name, types.StringType):
            raise Exception("状态名称必须为字符串")
        elif not name.islower():
            raise Exception("状态名称必须都为小写字母！")
        else:
            if self.instant:
                action = self._send_local_state(State(name, data))
                if isinstance(action, dict) and 'action' in action.keys():
                    if action['action'] != "none":
                        if hasattr(self, action['action'] + '_action'):
                            func = getattr(self, action['action'] + '_action')
                            argspec = inspect.getargspec(func)
                            args = argspec.args
                            if len(args) == 1:
                                func()
                            elif len(args) == 2:
                                func(action['data'])
                            else:
                                raise Exception("action处理函数参数不能超过2个(包括self)!")
                        else:
                            self.logger.warn("没有为返回action %s设置处理函数!" % action['action'])
            else:
                with self.mutex:
                    self._states.append(State(name, data))
        return self

    def _send_local_state(self, state):
        self.save_report()
        msg = {"id": self.id,
               "state": state.name,
               "data": state.data,
               "version": self.version
               }
        try:
            self.process_queue.put(msg)
        except Exception, e:
            self.logger.debug("error: %s", e)

        # print "Put MSG:%s" % msg
        key = "%s" % self.version
        action = {'action': 'none'}
        if key in self.shared_callback_action.keys():
            action = self.shared_callback_action.pop(key)
        return action

    def run(self):
        """
        开始运行,此方法应该在子类中被扩展

        :return:
        """
        threading.Thread(target=self.check_command, args=(), name='thread-check-command').start()
        self.logger.debug("启动Check Command线程!")
        if not self.instant:
            signal.signal(signal.SIGTERM, self.exit_handler)
            signal.signal(signal.SIGHUP, self.exit_handler)
            self.logger.info("注册信号处理函数!")
        self._run = True
        self.logger.info("开始运行!")

    def exit_handler(self, signum, frame):
        self.receive_exit_signal = True
        self.save_report()
        if signum == 0:
            self.logger.debug("检测到父进程终止，准备退出!")
        else:
            self.logger.debug("接收到Linux signal %s 信号，终止运行!" % signum)
        if "SLURM_JOB_ID" in os.environ.keys():
            time.sleep(3)
            sys.stderr.flush()
            sys.stdout.flush()
            self.logger.debug("开始检测SLURM STDERR输出...")
            slurm_error_path = os.path.join(self.work_dir, "%s_%s.err" % (self.name, os.environ["SLURM_JOB_ID"]))
            error_msg = "程序被终止运行,原因未知.."
            exceeded = False
            cancelled = False
            if not os.path.exists(slurm_error_path):
                self.logger.debug("未发现STDERR文件%s!" % slurm_error_path)
                return
            with open(slurm_error_path, "r") as f:
                f.seek(0, 2)
                size = os.path.getsize(slurm_error_path)
                point = 5000 if size > 5000 else size
                f.seek(-point, 2)
                lines = f.readlines()
                for line in lines:
                    if re.match(r"^slurmstepd:", line):
                        error_msg += line
                        if re.search(r"memory limit", line):
                            exceeded = True
                        elif re.search(r"CANCELLED", line):
                            cancelled = True
            self.logger.debug("检测完毕!")
            if exceeded:
                self.logger.info("检测到内存使用超过申请数被系统杀死!")
                self.add_state("memory_limit", error_msg)
            elif cancelled:
                self.logger.info("检测到任务被取消!")
                self.add_state("cancelled", error_msg)
            else:
                self.add_state("killed", error_msg)
        else:
            self.set_error("检测到终止信号，但是原因未知！")

    def resource_record(self, command):
        """
        记录运行资源

        :param command: 需要监控的Commnd对象
        :return:
        """
        if command == "":
            filepath = os.path.join(self.work_dir, "All_resource.txt")
            resource = self.check_resource()
            if self.is_end:
                return
        else:
            filepath = os.path.join(self.work_dir, command.name+"_resource.txt")
            if command.is_running:
                resource = command.check_resource()
            else:
                return
        time_now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
        if resource:
            with open(filepath, "a") as f:
                if len(resource) > 0:
                    r = resource.pop(0)
                    f.write("%s\tmain_pid:%s\tcpu_percent:%s\tmemory_rss:%s\tmemory_vms:%s\tcmd:%s\n" %
                            (time_now, r[0], r[2], friendly_size(r[3][0]), friendly_size(r[3][1]), r[1]))
                    if len(resource) > 0:
                        for r in resource:
                            f.write("\t\t\tChild pid:%s\tcpu_percent:%s\tmemory_rss:%s\tmemory_vms:%s\tcmd:%s\n"
                                    % (r[0], r[2], friendly_size(r[3][0]), friendly_size(r[3][1]), r[1]))

    def _set_resource_use(self, cpu, memory_rss, memory_vms):

        if cpu > self._cpu_max_use:
            self._cpu_max_use = cpu
        if memory_rss > self._mem_max_rss:
            self._mem_max_rss = memory_rss
        if memory_vms > self._mem_max_vms:
            self._mem_max_vms = memory_vms
        self._cup_avg_use = (self._set_resource_count * self._cup_avg_use + cpu) / (self._set_resource_count + 1)
        self._mem_avg_rss = (self._set_resource_count * self._mem_avg_rss + memory_rss) / (self._set_resource_count + 1)
        self._mem_avg_vms = (self._set_resource_count * self._mem_avg_vms + memory_vms) / (self._set_resource_count + 1)
        self._set_resource_count += 1

    def check_resource(self):
        try:
            main_process = psutil.Process(os.getpid())
            child_processes = main_process.children(recursive=True)
            pid = main_process.pid
            cmd = " ".join(main_process.cmdline())
            cpu_percent = main_process.cpu_percent(interval=1)
            memory_info = main_process.memory_info()
            if len(child_processes) == 0:
                self._set_resource_use(cpu_percent, memory_info[0], memory_info[1])
                return [(pid, cmd, cpu_percent, memory_info)]
            else:
                resource_list = [(pid, cmd, cpu_percent, memory_info)]
                cpu = cpu_percent
                memory_rss = memory_info[0]
                memory_vms = memory_info[1]
                for child in child_processes:
                    if child not in self._all_processes:
                        self._all_processes.append(child)
                    pid = child.pid
                    cmd = " ".join(child.cmdline())
                    cpu_percent = child.cpu_percent()
                    memory_info = child.memory_info()
                    resource_list.append((pid, cmd, cpu_percent, memory_info))
                    cpu += cpu_percent
                    memory_rss += memory_info[0]
                    memory_vms += memory_info[1]
                self._set_resource_use(cpu, memory_rss, memory_vms)
                return resource_list
        except Exception, e:
            self.logger.debug("Tool %s 监控资源时发生错误: %s" % (self.name, e))
            return None

    def get_resource_use(self):
        """
        返回最近一次运行使用资源情况

        :return:
        """
        r = (self._cpu_max_use, self._cup_avg_use, self._mem_max_rss, int(self._mem_avg_rss),
             int(self._mem_max_vms), self._mem_avg_vms)
        return r

    def exit_action(self, data):
        """
        处理远程agent端发回的exit action

        :param data:  退出指令说明
        :return:
        """
        self.logger.info("接收到action退出指令:%s" % str(data))
        self.exit()

    def rerun_action(self, data):
        """
        处理远程agent端发回的rerun指令

        :param data:  退出指令说明
        :return:
        """
        self.logger.info("接收到rerun退出指令:%s" % str(data))
        self._rerun = True
        self.kill_all_commonds()
        os.chdir(self.work_dir)
        script = sys.executable
        args = " ".join(sys.argv)
        self.logger.warning("终止主线程运行...")
        stop_thread(self.main_thread)
        self.logger.warning("开始重新运行...")
        self.exit_signal = True
        os.system("%s %s" % (script, args))
        self.exit()

    def end(self):
        """
        设置Tool已经完成，设置完成后,actor将在发送finish状态后终止发送信息

        :return:
        """
        self.save_report()
        self.save_output()
        self.remote.save()
        self.add_state('finish')
        self._end = True
        self.logger.info("Tool程序运行完成")

    def exit(self, status=1):
        """
        不发送任何信号，立即退出程序运行并终止所有命令运行

        :param status: 退出运行时的exitcode
        :return:
        """
        self.save_report()
        self.kill_all_commonds()
        if not self.is_end and self.main_thread.is_alive():
            self._end = True
            self.exit_signal = True
            # stop_thread(self.main_thread)
            # os._exit(status)
            if status > 0:
                # os.kill(os.getpid(), signal.SIGTERM)
                os.system("kill -9 %s" % os.getpid())
            os._exit(status)

        # if self.main_thread.is_alive():
        #     stop_thread(self.main_thread)
        self._end = True
        self.exit_signal = True
        self.logger.info("退出运行")
        os._exit(status)

    def save_output(self):
        """
        将输出参数Option对象写入Pickle文件，便于传递给远程Agent

        :return:
        """
        path = os.path.join(self.work_dir, self.name + "_output.pk")
        output = {}
        for name, option in self._options.items():
            if option.type == 'outfile':
                output[name] = self.option(name)
                output[name].option = None
        with open(path, "w") as f:
            pickle.dump(output, f)
        return path

    def save_report(self):
        path = os.path.join(self.work_dir, self.name + "_run_report.pk")
        all_resource_use = self.get_resource_use()
        output = {
            "process_id": os.getpid(),
            "sub_process_num": len(self._all_processes),
            "max_cpu_use": all_resource_use[0],
            "max_rss": all_resource_use[2],
            "average_cpu_use": all_resource_use[1],
            "average_rss": all_resource_use[3],
            "max_vms": all_resource_use[4],
            "average_vms": all_resource_use[5],
            "cmd": []
        }
        for name, command in self.commands.items():
            resource_use = command.get_resource_use()
            cmd_list = {
                "name": name,
                "cmd": command.cmd,
                "start_time": command.start_time,
                "end_time": command.end_time,
                "run_times": command.run_times,
                "main_pid": command.pid,
                "sub_process_num": command.sub_process_num,
                "max_cpu_use": resource_use[0],
                "max_rss": resource_use[2],
                "average_cpu_use": resource_use[1],
                "average_rss": resource_use[3],
                "return_code": command.return_code,
                "max_vms": resource_use[4],
                "average_vms": resource_use[5],
            }
            output["cmd"].append(cmd_list)
        with open(path, "w") as f:
            try:
                pickle.dump(output, f)
            except MemoryError, e:
                exstr = traceback.format_exc()
                print >> sys.stderr, exstr
                print >> sys.stderr, e
                sys.stderr.flush()
                self.add_state("memory_limit", "检测到内存使用时出现错误!")
            except Exception, e:
                exstr = traceback.format_exc()
                print >> sys.stderr, exstr
                print >> sys.stderr, e
                sys.stderr.flush()
        return path

    def load_config(self):
        """
        从Config对象中加载属性

        :return:
        """
        c = Config()
        for name in vars(self.config).keys():
            if hasattr(self, name):
                setattr(self, name, getattr(self.config, name))
            else:
                setattr(c, name, getattr(self.config, name))

        for option in self._options.values():
            option.bind_obj = self
            if isinstance(option.value, FileBase):
                option.value.option = option

    def set_error(self, error_data, variables=None, code="001"):
        """
        设置错误状态,设置完成后,actor将在发送error状态后终止发送信息

        :param error_data:
        :param variables:
        :param code:
        :return:
        """
        self.save_report()
        if isinstance(error_data, CodeError):
            e = error_data
        else:
            e = RunningError(str(error_data), variables, code)
            e.bind_object = self
        self.add_state('error', e.json())
        self.exit_signal = True
        self.logger.error("运行出错:%s" % e)
        sys.exit(0)

    @staticmethod
    def set_environ(**kwargs):
        """
        设置环境变量,清除环境变量使用os.unset_environ

        :param kwargs:  一个或多个带名称的参数,参数名为变量名 value为变量值
        """
        for (key, value) in kwargs.items():
            if key not in os.environ.keys():
                os.environ[key] = value
            else:
                os.environ[key] = value + ":" + os.environ[key]

    @staticmethod
    def unset_environ(*varname):
        """
        清除环境变量

        :param varname:  环境变量名
        :return:   返回被删除的环境变量值
        """
        if varname in os.environ.keys():
            return os.environ.pop(varname)

    @staticmethod
    def load_package(path):
        """
        动态加载mbio.packges下自定义模块的类对象

        :param path: 自动以模块path路径
        :return: class对象
        """
        path = re.sub(r'[^a-zA-Z0-9\._]', '', path)
        name = path.split(".").pop()
        l = name.split("_")
        l = [el.capitalize() for el in l]
        class_name = "".join(l)
        imp = importlib.import_module("mbio.packages." + path)
        return getattr(imp, class_name)

    def kill_all_commonds(self):
        """
        杀死所有正在运行的Commond
        :return:
        """
        for name, command in self._commands.items():
            if command.is_running:
                self.logger.warning("终止命令%s运行..." % name)
                command.kill()
                os.system("kill -9 %s" % command.pid)
        main_process = psutil.Process(os.getpid())
        childs = main_process.children(recursive=True)
        for p in childs:
            p.kill()
            os.system("kill -9 %s" % p.pid)

    def check_command(self):
        while not self.is_end:
            if not os.path.exists(self.work_dir):
                self.logger.error("工作目录被删除，退出运行!")
                self.exit()
            try:
                if not self.main_thread.is_alive():
                    break
                if self.is_end or self.exit_signal:
                    break
                self.resource_record("")
                for name, cmd in self.commands.items():
                    if cmd.is_running:
                        self.resource_record(cmd)
            except Exception, e:
                exstr = traceback.format_exc()
                print >> sys.stderr, exstr
                print >> sys.stderr, e
                sys.stderr.flush()
            time.sleep(15)

    def check_options(self):
        pass

    def run_check_options(self):
        try:
            self.check_options()
        except OptionError as e:
            e.bind_object = self
            self.add_state('error', e.json())
            self.exit_signal = True
            self.logger.error("运行出错:%s" % e)

    def download_from_s3(self, from_file, to_path="download/", cover=True):
        """
        从s3对象存储下载数据到本地
        :param from_file: 需要下载的文件路径或文件路径, 必须是类似s3region://bucket/key写法。
        因为对象存储中没有文件夹的概念，需要下载文件夹必须使用"/"结尾，以明确表明下载的是文件夹
        :param to_path: 下载文件相对于当前工作目录的存放目录。
        当路径为"/"结尾时，表示下载文件存放在此文件夹下，否者为下载完整路径。
        当from_file为文件夹时，此参数也必须以"/"结尾。目录层级与下载的s3目录层级结构相同。
        默认情况下放置在当前模块工作目录的download目录下。
        :param cover: 对已存在的文件是否覆盖
        :return:
        """
        if re.match(r"^/|^\.\.", to_path):
            raise Exception("不能使用绝对路径或切换到其他目录!")
        if os.path.basename(to_path) == ".":
            raise Exception("目标文件不能叫\".\"!")
        target_dir = False
        if re.match(r"/$", to_path):
            target_dir = True
        # self.s3transfer = S3TransferManager()
        # self.s3transfer.base_path = self.work_dir
        m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
        if not m:
            raise Exception("下载路径%s格式不正确!" % from_file)
        else:
            region = m.group(1)
            bucket_name = m.group(2)
            key_name = m.group(3)
            if re.match(r"/$", key_name):
                if not target_dir:
                    raise Exception("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
                conn = self.s3transfer.config.get_rgw_conn(region, bucket_name)
                bucket = Bucket(connection=conn, name=bucket_name)
                for key in bucket.list(prefix=key_name):
                    source = os.path.join(from_file, key.name)
                    target = os.path.join(target_dir, os.path.relpath(key.name, key_name))
                    self.s3transfer.add(source, target)
            else:
                if not target_dir:  # 处理已存在文件的情况
                    target = os.path.join(self.work_dir, to_path)
                    if os.path.exists(target):
                        if cover:
                            if os.path.isdir(target):
                                shutil.rmtree(target)
                            else:
                                os.remove(target)
                        else:
                            raise Exception("目标文件夹%s已经存在!" % target)
                else:
                    target = os.path.join(self.work_dir, to_path, os.path.basename(key_name))
                self.s3transfer.add(from_file, target)
            self.s3transfer.wait(end=True)

    def upload_to_s3(self, from_file, to_path, cover=True):
        """
        从本地上传数据到S3对象存储
        :param from_file: 需要上传的文件，相对于当前模块工作目录的相对路径，可以是文件或文件夹，也可使用通配符。
        当为文件夹时会上传文件夹的目录结构。通配符匹配到的文件夹也会上传目录结构。
        :param to_path: 上传文件的路径，必须是类似s3region://bucket/key写法。
        当路径为"/"结尾时，表示上传文件存放在此文件夹下，否者为上传完整路径。
        当上传文件使用通配符或为文件夹时，此参数必须为"/"结尾，且通配符匹配到的文件名不能重复，否则会相互覆盖。
        :param cover: 对已存在的文件是否覆盖
        :return:
        """
        if not re.match(r"^\w+://\S+/.+$", to_path):
            raise Exception("上传路径%s格式不正确!" % to_path)
        if os.path.basename(to_path) == ".":
            raise Exception("目标文件不能叫\".\"!")
        target_dir = False
        if re.match(r"/$", to_path):
            target_dir = True
        # self.s3transfer = S3TransferManager(base_path=self.work_dir, overwrite=cover)
        source = os.path.join(self.work_dir, from_file)
        for f in glob.glob(os.path.join(source)):
            if target_dir:
                target = os.path.join(to_path, os.path.basename(f))
            else:
                target = to_path
            if os.path.isdir(f):
                for root, dirs, files in os.walk(f):
                    rel_path = os.path.relpath(root, f)
                    for i_file in files:
                        if rel_path == ".":
                            key_path = os.path.join(target, i_file)
                        else:
                            key_path = os.path.join(target, rel_path, i_file)
                        i_file_path = os.path.join(root, i_file)
                        self.s3transfer.add(i_file_path, key_path)
            else:
                self.s3transfer.add(f, target)
        self.s3transfer.wait(end=True)
