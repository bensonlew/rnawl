#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from mbio.packages.metagbin.common_function import link_dir
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
import re

class UnzipAgent(Agent):
    """
    解压一个文件或者文件夹
    """
    def __init__(self, parent):
        super(UnzipAgent, self).__init__(parent)
        options = [
            {"name": "file_path", "type": "string"},  # 输入的文件名称
            {"name": "file_dir", "type": "string"},  # 输入的文件夹名称
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("file_path") and not self.option("file_dir"):
            raise OptionError("需要输入文件")
        if self.option("file_path") and self.option("file_dir"):
            raise OptionError("不能同时输入文件和文件夹")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(UnzipAgent, self).end()


class UnzipTool(Tool):
    """
    解压
    """
    def __init__(self, config):
        super(UnzipTool, self).__init__(config)
        self.path = ''
        self.file = ''

    def unzip_run(self):
        if self.option('file_path'):
            self.logger.info('处理单个')
            path = self.option('file_path')
            self.path = os.path.dirname(path)
            self.file = os.path.basename(path)
            self.run_unzip()
        else:
            self.logger.info('处理多个')
            self.path = self.option('file_dir')
            if not os.path.isdir(self.path):
                self.set_error("%s并非文件路径", variables=(self.path))
            files = os.listdir(self.path)
            if len(files) == 0:
                self.set_error("%s为空文件路径", variables=(self.path))
            for file in files:
                self.file = file
                self.run_unzip()

    def run_unzip(self):
        """
        解压文件
        :return:
        """
        if os.path.exists(self.work_dir + "/unzip_dir"):
            pass
        else:
            os.mkdir(self.work_dir + "/unzip_dir")
        file_path = os.path.join(self.path, self.file)
        if re.search(r'.tar.gz', self.file):
            file_name = self.work_dir + "/unzip_dir/"+ ".".join(self.file.split('.')[0:-2])
            unzip_cmd = "zcat " +  file_path + " > " + file_name
            try:
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error")
                raise Exception("unzip error")
        elif re.search(r'.gz', self.file):
            file_name = self.work_dir + "/unzip_dir/"+ ".".join(self.file.split('.')[0:-1])
            unzip_cmd = "zcat " +  file_path + " > " + file_name
            try:
                subprocess.check_output(unzip_cmd, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error")
                raise Exception("unzip error")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        link_dir(self.work_dir + "/unzip_dir", self.output_dir)

    def run(self):
        super(UnzipTool, self).run()
        self.unzip_run()
        self.set_output()
        self.end()