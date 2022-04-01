# -*- coding: utf-8 -*-
# __author__ :zhouxuan
# last_modify: 20180119

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import math
import shutil


class FastqUngzAgent(Agent):
    def __init__(self, parent):
        super(FastqUngzAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "string"},
            {"name": "sample_name", "type": "string"},
            {"name": "direction", "type": "string"},
            {"name": "result_path", "type": "string"},
            {"name": "pipeline", "type": "string"} #metaasv
        ]
        self.add_option(options)
        self.step.add_steps("fastq_ungz")
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
            raise OptionError("必须输入fastq文件", code="34001401")
        if not self.option('sample_name'):
            raise OptionError("必须输入样品名", code="34001402")
        if not self.option('direction'):
            raise OptionError("必须输入fastq文件的方向：1或2", code="34001403")
        if not self.option('result_path'):
            raise OptionError("必须输入结果存放路径", code="34001404")
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
        self._version = "v1.0"
        self.sh_path =  self.config.PACKAGE_DIR +"/sequence/scripts/unzip.sh"

    def run(self):
        super(FastqUngzTool, self).run()
        if self.option("pipeline") == "metaasv":
            self.run_unzip_metaasv()
        else:
            self.run_unzip()
        self.end()

    def run_unzip(self):
        sample = self.option("sample_name")
        fastq = self.option("fastq")
        f_name = sample + "." + self.option("direction") + ".fq"
        new_fastq = os.path.join(self.option("result_path"), f_name)
        unzip_cmd = "zcat " +  fastq + " > " + new_fastq
        self.logger.info("start unzip")
        self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
        try:
            subprocess.check_output(unzip_cmd, shell=True)
            self.logger.info("unzip done")
        except subprocess.CalledProcessError:
            self.set_error("%s:unzip error" %(f_name), code="34001401") #添加样本名称到报错信息中  @20190909张清臣
            raise Exception("unzip error")

    def run_unzip_metaasv(self):
        sample = self.option("sample_name")
        fastq = self.option("fastq")
        f_name = sample
        new_fastq = os.path.join(self.option("result_path"), f_name)
        unzip_cmd = "zcat " +  fastq + " > " + new_fastq
        self.logger.info("start unzip")
        self.logger.info("command: %s" % unzip_cmd)  # 添加command log MODIFIED GHD @20180119
        try:
            subprocess.check_output(unzip_cmd, shell=True)
            self.logger.info("unzip done")
        except subprocess.CalledProcessError:
            self.set_error("%s:unzip error" %(f_name), code="34001401") #添加样本名称到报错信息中  @20190909张清臣
            raise Exception("unzip error")