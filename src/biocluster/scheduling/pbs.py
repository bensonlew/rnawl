# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .job import Job
import os
import gevent
import re


class PBS(Job):
    """
    openPBS任务调度系统,
    用于生产和管理PBS任务
    """
    def __init__(self, agent):
        super(PBS, self).__init__(agent)
        self.master_ip = agent.config.JOB_MASTER_IP

    def create_file(self):
        """
        生成PBS脚本用于投递

        :return:
        """
        file_path = os.path.join(self.agent.work_dir, self.agent.name + ".pbs")
        script = os.path.abspath(os.path.dirname(__file__) + "/../../../bin/runtool.py")
        cpu, mem = self.agent.get_resource()
        with open(file_path, "w") as f:
            f.write("#PBS -N %s\n" % self.agent.fullname)
            f.write("#PBS -l nodes=1:ppn=%s\n" % cpu)
            f.write("#PBS -l mem=%s\n" % mem)
            f.write("#PBS -d %s\n" % self.agent.work_dir)
            f.write("#PBS -q %s\n\n" % self.agent.queue)
            f.write("cd %s\n\n" % self.agent.work_dir)
            f.write("%s %s\n" % (script, self.agent.name))

        return file_path

    def submit(self):
        """
        提交PBS任务,并返回Jobid

        :return: jobid
        """
        super(PBS, self).submit()
        pbs_file = self.create_file()
        output = os.popen('ssh -o GSSAPIAuthentication=no %s "/opt/torque/bin/qsub %s"' % (self.master_ip, pbs_file))
        text = output.read()
        if re.match(r'Maximum number', text):
            self.agent.logger.warn("到达最大任务书，30秒后尝试再次投递!")
            gevent.sleep(30)
            self.submit()
        else:
            m = re.search(r'(\d+)\..*', text)
            if m:
                self.id = m.group(1)
                return self.id
            else:
                self.agent.logger.warn("任务投递系统出现错误:%s，30秒后尝试再次投递!\n" % output)
                gevent.sleep(30)
                self.submit()

    def delete(self):
        """
        删除当前任务

        :return: None
        """

        if self.check_state():
            os.system('ssh -o GSSAPIAuthentication=no root@%s "/opt/torque/bin/qdel -p %s"' % (self.master_ip, self.id))

    def check_state(self):
        """
        检测任务状态

        :return: string 返回任务状态代码 如果任务不存在 则返回False
        """
        output = os.popen('ssh -o GSSAPIAuthentication=no %s "/opt/torque/bin/qstat -f %s"' % (self.master_ip, self.id))
        text = output.read()
        m = re.search(r"job_state = (\w+)", text)
        if m:
            return m.group(1)
        else:
            return False
