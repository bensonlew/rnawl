# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from ..core.singleton import singleton
from ..config import Config
import gevent
import importlib
import datetime
from ..core.watcher import Watcher
import copy
import subprocess
import re


@singleton
class JobManager(object):
    """
    管理单个流程中的所有Job对象::

        singleton单例类，只会实例化一次
    """
    def __init__(self):
        config = Config()
        self.config = config
        self.run_jobs = []
        self.default_mode = config.JOB_PLATFORM
        self.max_job_number = config.MAX_JOB_NUMBER
        self.max_cpu_used = config.MAX_CPU_USED
        self.max_memory_used = config.MAX_MEMORY_USED
        self.queue_jobs = []
        Watcher().add(self._watch_waiting_jobs, 10)
        Watcher().add(self._watch_processing_jobs, 10)
        self._slurm_queue = []

    def add_job(self, agent):
        """
        根据agent信息生成一个Job子类对象，并且加入job队列开始运行

        :return:
        """
        mode = self.default_mode.lower()
        if agent.mode.lower() != "auto":
            mode = agent.mode.lower()
        module = importlib.import_module("biocluster.scheduling.%s" % mode)
        job = getattr(module, mode.upper())(agent)
        return job

    def run_job(self, job):
        agent = job.agent
        workflow = agent.get_workflow()
        paused = False
        while workflow.pause:
            if not paused:
                self.logger.info("流程处于暂停状态，排队等待恢复运行!")
            paused = True
            gevent.sleep(3)

        if len(self.get_unfinish_jobs()) >= self.max_job_number:
            self.queue_jobs.append(job)
            agent.logger.warning("任务队列达到最大上限%s个，排队等待运行!" % self.max_job_number)
            agent.is_wait = True
        elif self.get_running_used_cpu() > self.max_cpu_used:
            self.queue_jobs.append(job)
            agent.logger.warning("任务队列使用的CPU核心数超过最大限制%s个，排队等待运行!" % self.max_cpu_used)
            agent.is_wait = True
        elif self.get_running_used_memory() > self.max_memory_used:
            self.queue_jobs.append(job)
            agent.logger.warning("任务队列使用的内存数超过最大限制%sG，排队等待运行!" % self.max_memory_used)
            agent.is_wait = True
        elif len(self.run_jobs) > 10:
            mode = self.default_mode.lower()
            if agent.mode.lower() != "auto":
                mode = agent.mode.lower()
            if mode == "slurm":
                if len(self.get_slurm_queue_list()) > 100:
                    self.queue_jobs.append(job)
                    agent.logger.warning("slurm所有任务排队数量超过100个，暂停投递!")
                    agent.is_wait = True
                elif self.get_slurm_queue_job_number() > 10:
                    self.queue_jobs.append(job)
                    agent.logger.warning("当前流程slurm排队数量超过10个，暂停投递!")
                    agent.is_wait = True
                else:
                    self._run_job(job)
            else:
                self._run_job(job)
        else:
            self._run_job(job)

    def _run_job(self, job):
        agent = job.agent
        mode = self.default_mode.lower()
        if agent.mode.lower() != "auto":
            mode = agent.mode.lower()
        agent.is_wait = False
        agent.logger.info("开始投递远程任务!")
        self.run_jobs.append(job)
        if job.submit():
            agent.logger.info("任务投递成功,任务类型%s , ID: %s!" % (mode, job.id))
        else:
            agent.logger.error("任务投递失败!")
            agent.set_error("Tool任务投递失败!", None, "004")

    def get_all_jobs(self):
        """
        获取所有任务对象

        :return: list  Job子类对象列表
        """
        jobs = copy.copy(self.run_jobs)
        jobs.extend(self.queue_jobs)
        return jobs

    def get_unfinish_jobs(self):
        """
        获取未完成的任务对象

        :return: list  Job子类对象列表
        """
        un_done = []
        for job in self.run_jobs:
            if not job.is_end:
                un_done.append(job)
        return un_done

    def get_job(self, agent):
        """
        根据agent对象找出其对应的Job子类，如果没有找到，则返回False

         :return:  Job子类对象
        """
        jobs = copy.copy(self.run_jobs)
        jobs.extend(self.queue_jobs)
        for job in jobs:
            if agent is job.agent:
                return job
        return False

    def get_running_used_cpu(self):
        """
        获取正在运行的Job申请的CPU数
        :return:
        """
        count = 0
        for job in self.run_jobs:
            if not job.is_end:
                count += job.cpu
        return count

    def get_running_used_memory(self):
        """
        获取正在运行的JOB申请的内存数,单位G
        :return:
        """
        count = 0
        for job in self.run_jobs:
            if not job.is_end:
                count += job.memory
        return count/(1024*1024*1024)

    def get_slurm_queue_job_number(self):
        count = 0
        job_list = []
        for job in self.get_unfinish_jobs():
            job_list.append(job.id)
        queue_list = self.get_slurm_queue_list()
        for i in job_list:
            if i in queue_list:
                count += 1
        return count

    def get_slurm_queue_list(self):
        if self._slurm_queue and (datetime.datetime.now() - self._slurm_queue[0]).seconds < 60:
            return self._slurm_queue[1]
        id_list = []
        cmd = "squeue  -h -t PD -u {}".format(self.config.wpm_user)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        text = process.communicate()[0]
        if text:
            for line in text.split("\n"):
                if line:
                    id_list.append(re.split("\s+", line)[1])
        self._slurm_queue = [datetime.datetime.now(), id_list]
        return id_list

    def _watch_waiting_jobs(self):
        """
        监控等待排队的任务，满足任务运行条件时开始运行任务

        :return:
        """
        # while True:
        #     gevent.sleep(30)
        if self.default_mode.lower() == "slurm":
            if self.get_slurm_queue_job_number() > 10:
                return
        jobs = copy.copy(self.queue_jobs)
        for queue_job in jobs:
            if len(self.get_unfinish_jobs()) < self.max_job_number and self.get_running_used_cpu() < self.max_cpu_used\
                    and self.get_running_used_memory() < self.max_memory_used:
                queue_job.agent.is_wait = False
                queue_job.agent.logger.info("开始投递任务!")
                self.run_jobs.append(queue_job)
                mode = self.default_mode.lower()
                if queue_job.agent.mode.lower() != "auto":
                    mode = queue_job.agent.mode.lower()
                if queue_job.submit():
                    queue_job.agent.logger.info("任务投递成功,任务类型%s , ID: %s!" % (mode, queue_job.id))
                    self.queue_jobs.remove(queue_job)
                else:
                    queue_job.agent.logger.error("任务投递失败!")
                    queue_job.agent.set_error("任务投递失败!", None, "005")
                gevent.sleep(0)

    def _watch_processing_jobs(self):
        """
        监视本地进程模式时进程意外终止的问题

        :return:
        """
        if len(self.get_all_jobs()) == 0:
            return
        is_process = False
        for running_job in self.get_unfinish_jobs():
            if hasattr(running_job, "process"):
                is_process = True
                if not running_job.process.is_alive():
                    gevent.sleep(10)
                    if not running_job.is_end:
                        running_job.agent.fire("error", "任务%s进程意外结束，请查看运行日志!" % running_job.agent.id)

        if is_process is False:
            return "exit"

    def remove_all_jobs(self):
        """
        删除所有未完成任务

        :return:
        """
        for job in self.run_jobs:
            if not job.is_end:
                job.delete()


class Job(object):
    """
    Job基类,
    用于扩展各种集群调度平台
    """
    def __init__(self, agent):
        self.agent = agent
        self.id = 0
        self._end = False
        self.submit_time = datetime.datetime.now()
        self.state = ""

    @property
    def cpu(self):
        return self.agent.cpu

    @property
    def memory(self):
        return self.agent.memory

    @property
    def is_end(self):
        """
        返回是否已经完成

        :return: bool
        """
        return self._end

    def submit(self):
        """
        提交任务,需在子类中重写

        :return:
        """
        self.submit_time = datetime.datetime.now()

    def resubmit(self):
        """
        删除任务并重新提交

        :return: None
        """
        self.delete()
        self.id = 0
        self.submit_time = None
        self._end = False
        self.submit()

    def delete(self):
        """
        删除任务,需要在子类中重写此方法

        :return:
        """
        pass

    def check_state(self):
        """
        检测任务状态,需要在子类中重写此方法

        :return: string 返回任务状态代码 如果任务不存在 则返回False
        """
        pass

    def set_end(self):
        """
        将Job状态设置为已完成
        :return:
        """
        self._end = True

    def is_queue(self):
        """
        判断是否正在排队
        :return:
        """
        pass

    def is_running(self):
        """
        判断是否正在运行
        :return:
        """
        pass

    def is_error(self):
        """
        判断是否出现错误
        :return:
        """
        pass

    def is_completed(self):
        """
        判断任务是否完成
        :return:
        """
        pass

    def memory_limit(self):
        pass
