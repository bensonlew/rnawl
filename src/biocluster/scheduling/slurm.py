# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .job import Job
import os
import gevent
import re
import subprocess


class SLURM(Job):
    """
    openPBS任务调度系统,
    用于生产和管理PBS任务
    """
    def __init__(self, agent):
        super(SLURM, self).__init__(agent)
        self.master_ip = agent.config.JOB_MASTER_IP
        self.count = 1

    def create_file(self):
        """
        生成PBS脚本用于投递

        :return:
        """
        file_path = os.path.join(self.agent.work_dir, self.agent.name + ".sbatch")
        script = os.path.abspath(os.path.dirname(__file__) + "/../../../bin/runtool.py")
        cpu, mem = self.agent.get_resource()
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
        with open(file_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH -c {}\n".format(cpu))
            f.write("#SBATCH -D %s\n" % self.agent.work_dir)
            f.write("#SBATCH -n 1\n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -J {}\n".format(self.agent.fullname))
            f.write("#SBATCH -p %s\n" % self.agent.queue)
            f.write("#SBATCH --mem={}\n".format(mem))
            f.write("#SBATCH -o {}/{}_%j.out\n".format(self.agent.work_dir, self.agent.name))
            f.write("#SBATCH -e {}/{}_%j.err\n".format(self.agent.work_dir, self.agent.name))
            f.write("cd {}\n\n".format(self.agent.work_dir))
            f.write("{} {} {}\n".format("python", script, self.agent.name))
            path = os.path.abspath(os.path.dirname(__file__) + "/../../../bin")
            f.write("ssh -i {}/id_rsa root@localhost \"nohup flock -xn /dev/shm/freemem.lock"
                    " -c 'sh {}/freemem.sh' > /dev/null &\"\n".format(path, path))

        return file_path

    def submit(self):
        """
        提交PBS任务,并返回Jobid

        :return: jobid
        """
        super(SLURM, self).submit()
        pbs_file = self.create_file()
        cmd = "sbatch {}".format(pbs_file)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        text = process.communicate()[0]
        if re.match(r'Maximum number', text):
            self.agent.logger.warn("到达最大任务数量，30秒后尝试再次投递!")
            gevent.sleep(30)
            self.submit()
        else:
            m = re.search(r'(\d+)$', text)
            if m:
                self.id = m.group(1)
                return self.id
            else:
                if re.search("invalid\smemory\sconstraint", text):
                    self.agent.logger.error("内存指定出错 {}。任务运行失败".format(text))
                    self.agent.fire("error", "内存指定出错 {}。任务运行失败".format(text))
                elif self.count < 10:
                    self.agent.logger.warn("任务投递系统出现问题，未能获取JobId:%s，等待2分钟检测远程任务是"
                                           "否开始运行!\n" % text)
                    gevent.sleep(120)
                    if self.agent.status == "R" or self.agent.status == "E" and self.agent.job.id:
                        return self.agent.job.id
                    else:
                        self.agent.logger.warn("任务投递系统出现错误:%s，尝试再次投递!\n" % text)
                        self.count += 1
                        return self.submit()
                else:
                    self.agent.logger.error("已重复投递10次任务，终止运行")
                    self.agent.fire("error", "已重复投递10次任务，终止运行")

    def delete(self):
        """
        删除当前任务

        :return: None
        """

        if self.check_state():
            cmd = "scancel {}".format(self.id)
            # process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            # process.communicate()
            os.system(cmd)

    def check_state(self):
        """
        检测任务状态

        :return: string 返回任务状态代码 如果任务不存在 则返回False
        """
        cmd = "scontrol show job {}".format(self.id)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        text = process.communicate()[0]
        m = re.search(r"JobState=(\w+)", text)
        if m:
            self.state = m.group(1)
            return self.state
        else:
            return False

    def is_queue(self):
        # https://slurm.schedmd.com/squeue.html
        if re.match(r"PENDING|CONFIGURING|PD|CF", self.state, re.IGNORECASE):
            self.agent.logger.error("检测任务当前状态为:%s，正在排队" % self.state)
            return True
        return False

    def is_running(self):
        # https://slurm.schedmd.com/squeue.html
        if re.match(r"RUNNING|COMPLETING|R|CG", self.state, re.IGNORECASE):
            self.agent.logger.error("检测任务当前状态为:%s，正在运行中 " % self.state)
            return True
        return False

    def is_error(self):
        # https://slurm.schedmd.com/squeue.html
        if re.match(r"BOOT_FAIL|BF|CANCELLED|CA|FAILED|F|NODE_FAIL|NF|PREEMPTED|PR|STOPPED|ST"
                    r"|SUSPENDED|S|TIMEOUT|TO", self.state, re.IGNORECASE):
            self.agent.logger.error("检测任务当前状态为:%s，运行出现错误 " % self.state)
            return True
        return False

    def is_completed(self):
        # https://slurm.schedmd.com/squeue.html
        if re.match(r"COMPLETED|CD|SPECIAL_EXIT|SE", self.state, re.IGNORECASE):
            self.agent.logger.error("检测任务当前状态为:%s，运行已经完成 " % self.state)
            return True
        return False

    def memory_limit(self):
        slurm_error_path = os.path.join(self.agent.work_dir, "%s_%s.err" % (self.agent.name, self.id))
        if os.path.exists(slurm_error_path):
            with open(slurm_error_path, "r") as f:
                f.seek(0, 2)
                size = os.path.getsize(slurm_error_path)
                point = 5000 if size > 5000 else size
                f.seek(-point, 2)
                lines = f.readlines()
                for line in lines:
                    if re.match(r"^slurmstepd:", line):
                        if re.search(r"memory limit", line):
                            self.agent.logger.info("检测到内存使用超过申请数被系统杀死: %s" % line.strip("\n"))
                            return True
                    elif re.search(r".*/slurm_script: line \d+: \d+ Segmentation fault", line):
                        self.agent.logger.info("检测到内存访问错误被系统杀死: %s" % line.strip("\n"))
                        return True
                    elif re.search(r"^MemoryError", line):
                        self.agent.logger.info("检测到内存访问错误报错: %s" % line.strip("\n"))
                        return True

    def cancelled(self):
        slurm_error_path = os.path.join(self.agent.work_dir, "%s_%s.err" % (self.agent.name, self.id))
        if os.path.exists(slurm_error_path):
            with open(slurm_error_path, "r") as f:
                f.seek(0, 2)
                size = os.path.getsize(slurm_error_path)
                point = 2000 if size > 2000 else size
                f.seek(-point, 2)
                lines = f.readlines()
                for line in lines:
                    if re.match(r"^slurmstepd:", line):
                        if re.search(r"CANCELLED", line):
                            self.agent.logger.info("检测到任务被杀死: %s" % line.strip("\n"))
                            return True
