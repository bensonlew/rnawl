# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""
任务命令监控
"""

import shlex
import subprocess
import os
import psutil
import threading
import datetime
import inspect
import time
from .core.function import friendly_size


class Command(object):
    """
    命令模块，在 :py:class:`biocluster.agent.Tool` 中调用add_commond方法将生成此类的一个对象.
    Command对象负责运行一个命令，并监控其运行过程

    :param name:  string 名称
    :param cmd:  需要执行的命令，可执行文件路径必须为相对于配置文件main.conf的software_dir的相对路径
    :param tool: 调用此对象的 :py:class:`biocluster.agent.Tool` 实例
    """
    def __init__(self, name, cmd, tool, default_return_code=0, ignore_error=False, shell=False):
        """
        初始化生成Command对象,将会新建一个线程负责运行cmd
        """
        super(Command, self).__init__()
        # Command.count += 1
        self._pid = ""
        self.cmd = cmd
        self._last_run_cmd = None
        self.tool = tool
        self.config = tool.config
        self.work_dir = tool.work_dir
        self.software_dir = self.config.SOFTWARE_DIR.rstrip("/")
        self._name = name
        self._subprocess = None
        self.threading_lock = threading.Lock()
        # self._psutil_process = None
        # self._all_processes = []
        self._is_error = False
        # self.process_status = None
        # self.check_sleep = check_sleep
        # self.max_sleep_time = max_sleep_time
        # self.max_run_times = max_run_times
        self._run_times = 0
        # self.to_rerun = False
        self._start_run_time = None  # 开始运行时间
        self._end_run_time = None  # 结束运行时间
        self._cpu_max_use = 0
        self._cup_avg_use = 0
        self._mem_max_rss = 0
        self._mem_avg_rss = 0
        self._mem_max_vms = 0
        self._mem_avg_vms = 0
        self._set_resource_count = 0
        self._all_processes = []
        self._has_run = False
        self._is_end = False
        self._thread = None
        self.default_return_code = default_return_code
        self.ignore_error = ignore_error
        self._is_sh = shell

    @property
    def is_end(self):
        return self._is_end

    @property
    def is_error(self):
        return self._is_error

    @property
    def pid(self):
        """
        命令运行的系统pid
        :return:
        """
        return self._pid

    @property
    def run_times(self):
        """
        返回运行次数
        :return:
        """
        return self._run_times

    @property
    def start_time(self):
        """
        返回最近一次运行的开始时间戳
        :return:
        """
        return self._start_run_time

    @property
    def sub_process_num(self):
        """
        返回子进程数量
        :return:
        """
        return len(self._all_processes)

    @property
    def end_time(self):
        """
        返回最近一次运行的结束时间戳
        :return:
        """
        return self._end_run_time

    @property
    def name(self):
        """
        命令名称
        :return:
        """
        return self._name

    @property
    def has_run(self):
        """
        返回是否已经开始运行
        """
        return self._has_run
        # if self._subprocess is not None or self.is_error:
        #     return True
        # else:
        #     return False

    @property
    def is_running(self):
        """
        返回是否在运行状态
        """
        if not self._subprocess or self.is_error:
            return False
        try:
            if self._subprocess.poll() is None:
                return True
            else:
                return False
        except Exception, e:
            print e
            return False

    @property
    def return_code(self):
        """
        获取命令退出时返回的状态编码
        """
        # self.tool.logger.debug("获取returncode: %s " % self._subprocess)
        if self._subprocess:
            self._subprocess.poll()
            return self._subprocess.returncode
        else:
            return None

    def set_cmd(self, cmd):
        """
        重设需要运行的命令

        :param cmd:  命令文本，可执行文件路径必须为相对于配置文件main.conf的software_dir的相对路径
        :return:
        """
        self.cmd = cmd

    # def get_psutil_processes(self):
    #     """
    #     返回当前命令及其子进程psutil.Process对象
    #
    #     :return: list of all child psutil_process objects
    #     """
    #     if self.is_running:
    #         try:
    #             if not self._psutil_process:
    #                 self._psutil_process = psutil.Process(self._pid)
    #                 self._all_processes = [self._psutil_process]
    #             chidrens = self._psutil_process.children(recursive=True)
    #             chidrens.insert(0, self._psutil_process)
    #             all_process = copy.copy(self._all_processes)
    #             for child in all_process:   # 删除已完成的进程
    #                 found = False
    #                 for p in chidrens:
    #                     if p.pid == child.pid:
    #                         found = True
    #                 if not found:
    #                     self._all_processes.remove(child)
    #             for p in chidrens:                # 添加新进程
    #                 found = False
    #                 for child in self._all_processes:
    #                     if p.pid == child.pid:
    #                         found = True
    #                 if not found:
    #                     self._all_processes.append(p)
    #         except Exception, e:
    #             self.tool.logger.debug("获取命令%s进程时发生错误: %s" % (self.name, e))
    #     return self._all_processes

    def _run(self):
        """
        运行命令

        :return: self
        """
        try:
            self.tool.logger.debug("开始启动Command %s 子进程subprocess " % self.name)
            if not self._is_sh:
                command = self.software_dir + "/" + self.cmd
                self._start_run_time = int(time.time())
                args = shlex.split(command)
                self._subprocess = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                                    env=os.environ, universal_newlines=True, bufsize=0)
            else:
                self._subprocess = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE,
                                                    stderr=subprocess.STDOUT,
                                                    env=os.environ, universal_newlines=True, bufsize=0)
            self.tool.logger.debug("Command %s subprocess生成: %s, PID: %s " %
                                   (self.name, self._subprocess, self._subprocess.pid))

        except Exception, e:
            if self._subprocess:
                self._subprocess.communicate()
            self._is_error = True
            self.tool.set_error(e)
        else:
            try:
                self._pid = self._subprocess.pid
                self._start_run_time = int(time.time())
                self._last_run_cmd = self.cmd
                func = None
                if hasattr(self.tool, self.name + '_check'):
                    func = getattr(self.tool, self.name + '_check')
                    argspec = inspect.getargspec(func)
                    args = argspec.args
                    if len(args) != 3:
                        raise Exception("状态监测函数参数必须为3个(包括self)!")
                count = 0
                tmp_file = os.path.join(self.work_dir, self._name + ".o")
                self.tool.logger.debug("开始检测Command %s输出" % self.name)
                with open(tmp_file, "w") as f:
                    f.write("%s\t运行开始\n" % datetime.datetime.now())
                    f.flush()
                    last_output = []
                    while True:
                        line = self._subprocess.stdout.readline()
                        if line == b''and self._subprocess.poll() is not None:
                                break
                        if line:
                            count += 1
                            if count < 5000:
                                f.write(line)
                            elif count == 5000:
                                f.write("输出过大，后续省略...\n")
                                f.flush()
                                last_output.append(line)
                            else:
                                last_output.append(line)
                                if len(last_output) > 100:
                                    last_output.pop(0)
                            if func is not None:
                                line = line.strip()
                                func(self, line)   # check function(toolself, command, line)  single line
                    if len(last_output) > 0:
                        f.writelines(last_output)
                        f.flush()
                    self._end_run_time = int(time.time())
                    use_time = self._end_run_time - self._start_run_time
                    # time_now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
                    f.write("%s\t运行结束，运行时长:%ss,exitcode:%s\n" % (datetime.datetime.now(), use_time,
                            self.return_code))
            except IOError, e:
                self._is_error = True
                self.tool.set_error("运行命令%s出错: %s" % (self.name, e))
        self._subprocess.communicate()
        self._subprocess.stdout.close()
        self._subprocess.poll()
        self._is_end = True
        return self

    def run(self):
        """
        新建线程开始运行

        :return: self
        """
        if self._is_sh:
            command = self.cmd
        else:
            command = self.software_dir + "/" + self.cmd
        self.tool.logger.info("命令内容为{}".format(command))
        # if self.is_running or self.has_run:
        #     raise OSError("命令已经运行，不能重复运行!")
        if not self._is_sh:
            if "|" in self.cmd or ">" in self.cmd or "<" in self.cmd:
                raise Exception("不能使用管道符或重定向符!")
        args = shlex.split(command)
        if not os.path.isfile(args[0]):
            self.tool.set_error("运行的命令文件不存在")
            raise Exception("你所运行的命令文件不存在，请确认!")
        self._thread = threading.Thread(target=self._run)
        self._thread.start()
        self._has_run = True
        self._run_times += 1
        return self

    def rerun(self):
        """
        重运行命令

        :return:
        """
        if self.cmd == self._last_run_cmd:
            if self._run_times > 3:
                raise Exception("重复运行相同的命令不能超过3次！命令:%s" % self.cmd)
            self.tool.logger.warning('重新运行了相同的命令')  # shenghe modified 20161215
        else:
            self._run_times = 0

        if self.is_running:
            self.kill()
        self._subprocess = None
        self._pid = ""
        self._is_error = False
        self._start_run_time = None
        self._end_run_time = None
        self._has_run = False
        self._is_end = False
        self._thread = None
        self.run()
        return self

    def kill(self):
        """
        结束命令运行
        :return:
        """
        if self.is_running:
            chidrens = psutil.Process(self.pid).children(recursive=True)
            for p in chidrens:
                p.kill()
            self._subprocess.kill()
            self.join_thread()
            self._is_error = True
            self._is_end = True

    def join_thread(self):
        if threading.currentThread() != self._thread:
            self._thread.join()

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
            if self.is_running:
                main_process = psutil.Process(self._pid)
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
            self.tool.logger.debug("Command %s 监控资源时发生错误: %s" % (self.name, e))
            return None

    def get_resource_use(self):
        """
        返回最近一次运行使用资源情况

        :return:
        """
        r = (self._cpu_max_use, self._cup_avg_use, self._mem_max_rss, int(self._mem_avg_rss),
             int(self._mem_max_vms), self._mem_avg_vms)
        return r
