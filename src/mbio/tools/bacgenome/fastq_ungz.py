# -*- coding: utf-8 -*-
# __author__ :zhouxuan
# last_modify: 20180119

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import math


class FastqUngzAgent(Agent):
    def __init__(self, parent):
        super(FastqUngzAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "string"},
            {"name": "sample_name", "type": "string"},
            {"name": "direction", "type": "string"},
            {"name": "lib_type", "type": "string"},
            {"name": "result_path", "type": "string"},
            {"name": "nozip", "type": "bool", "default": False}  # 是否需要解压，默认要解压，True的结果提供给fastp质控
        ]
        self.add_option(options)
        self.step.add_steps("fastq_ungz")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fastq_ungz.start()
        self.step.update()

    def stepfinish(self):
        self.step.fastq_ungz.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fastq'):
            raise OptionError("必须输入fastq文件", code="31401601")
        if not self.option('sample_name'):
            raise OptionError("必须输入样品名", code="31401602")
        if not self.option('direction'):
            raise OptionError("必须输入fastq文件的方向：1或2", code="31401603")
        if not self.option('result_path'):
            raise OptionError("必须输入结果存放路径", code="31401604")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        memory = os.path.getsize(self.option('fastq').split(" ")[0])
        n = memory / (1024 * 1024 * 1024)
        n = math.ceil(n)
        self._memory = '{}G'.format(int(n))

    def end(self):
        super(FastqUngzAgent, self).end()


class FastqUngzTool(Tool):
    def __init__(self, config):
        super(FastqUngzTool, self).__init__(config)


    def run(self):
        super(FastqUngzTool, self).run()
        self.run_unzip()
        self.end()

    def run_unzip(self):
        sample = self.option("sample_name")
        fastq = self.option("fastq")
        if not self.option('lib_type'):
            f_name = sample + "." + self.option("direction") + ".fq"
        else:
            f_name = sample + "_" + self.option('lib_type') + "." + self.option("direction") + ".fq"
        new_fastq = os.path.join(self.option("result_path"), f_name)
        unzip_cmd = ''
        if fastq.endswith(".gz") and fastq.endswith(".tar.gz"):
            path = self.work_dir + "/"+ os.path.basename(fastq).split(".tar.gz")[0]
            unzip_cmd = "tar -zxvf " + fastq + " > " + path
            try:
                self.logger.info("start unzip")
                self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info(path)
                if os.path.exists(new_fastq):
                    os.remove(new_fastq)
                os.link(path, new_fastq)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error", code="31401602")
        elif fastq.endswith(".gz") and not fastq.endswith(".tar.gz"):
            unzip_cmd = "zcat " +  fastq + " > " + new_fastq
            try:
                self.logger.info("start unzip")
                self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error", code="31401601")
        elif fastq.endswith("fq") or fastq.endswith("fastq"):  # 允许非压缩的fq或fastq文件作为输入 by GHD @20180329
            unzip_cmd = "cat " + fastq + " > " + new_fastq
            try:
                self.logger.info("start unzip")
                self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error", code="31401601")
        if self.option("nozip"):  # 如果使用fastp软件则不需要解压，只需要简单更名
            file_suffix = fastq.split(".")[-1]
            unzip_cmd = "cat " + fastq + " > " + new_fastq + "." + file_suffix
            try:
                self.logger.info("start unzip")
                self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error", code="31401601")
