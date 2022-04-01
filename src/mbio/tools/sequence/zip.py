# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess,shutil


class ZipAgent(Agent):
    """
    catreads:压缩一个文件，或者压缩一个文件夹
    version 1.0
    author: guhaidong
    last_modify: 2018.01.10
    """

    def __init__(self, parent):
        super(ZipAgent, self).__init__(parent)
        options = [
            {"name": "file_path", "type": "string"},  # 输入的文件名称
            {"name": "file_dir", "type": "string"},  # 输入的文件夹名称
            {"name": "method", "type": "string", "default": "tar"},  # 压缩方法
        ]
        self.add_option(options)
        self.step.add_steps('zips')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.zips.start()
        self.step.update()

    def step_end(self):
        self.step.zips.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("file_path") and not self.option("file_dir"):
            raise OptionError("需要输入文件", code="34003101")
        if self.option("file_path") and self.option("file_dir"):
            raise OptionError("不能同时输入文件和文件夹", code="34003102")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'


class ZipTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(ZipTool, self).__init__(config)
        self.path = ''
        self.file = ''
        self.gz_path = "../../../../../bin/gzip"

    def check(self):
        if not os.path.isdir(self.path):
            self.set_error("%s并非文件路径", variables=(self.path), code="34003101")
        if not os.path.isfile(os.path.join(self.path, self.file)):
            self.set_error("%s无此文件", variables=(os.path.join(self.path, self.file)), code="34003102")

    def zip_run(self):
        if self.option('file_path'):
            self.logger.info('处理单个')
            path = self.option('file_path')
            self.path = os.path.dirname(path)
            self.file = os.path.basename(path)
            self.check()
            output = self.zip_file()
        else:
            self.logger.info('处理多个')
            self.path = self.option('file_dir')
            if not os.path.isdir(self.path):
                self.set_error("%s并非文件路径", variables=(self.path), code="34003103")
            files = os.listdir(self.path)
            if len(files) == 0:
                self.set_error("%s为空文件路径", variables=(self.path), code="34003104")
            for file in files:
                self.file = file
                output = self.zip_file()

    def gz_run(self):
        if self.option('file_path'):
            self.logger.info('处理单个')
            path = self.option('file_path')
            self.file = os.path.basename(path)
            if os.path.exists(self.work_dir + '/' + self.file):
                os.remove(self.work_dir + '/' + self.file)
            shutil.copy2(path, self.work_dir + '/' + self.file)
            if os.path.exists(self.work_dir + '/' + self.file + '.gz'):
                os.remove(self.work_dir + '/' + self.file + '.gz')
            cmd = "{} {}".format(self.gz_path,self.work_dir + '/' + self.file)
            command = self.add_command("run_gz", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行run_gz完成")
            else:
                self.set_error("运行run_gz运行出错!")


    def zip_file(self):
        """
        压缩文件
        :return:
        """
        if os.path.isfile(os.path.join(self.output_dir, self.file + '.tar.gz')):
            os.remove(os.path.join(self.output_dir, self.file + '.tar.gz'))
        cmd = 'cd ' + self.path + '; tar czvf ' + self.file + '.tar.gz ' + self.file + '; cd -; mv ' + self.path + '/' + self.file + '.tar.gz  ' + self.output_dir
        self.logger.info("record command:")
        self.logger.info(cmd)
        try:
            os.system(cmd)
        except:
            self.set_error("解压缩%s失败", variables=(self.file), code="34003105")
        output = os.path.join(self.path, self.file + '.tar.gz')
        return output

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        if self.option('method') == "gz":
            if os.path.exists(self.output_dir + '/' + self.file + '.gz'):
                os.remove(self.output_dir + '/' + self.file + '.gz')
            os.link(self.work_dir + '/' + self.file + '.gz', self.output_dir + '/' + self.file + '.gz')
        self.logger.info("设置结果目录成功")

    def run(self):
        super(ZipTool, self).run()
        if self.option('method') == "tar":
            self.zip_run()
        elif self.option('method') == "gz":
            self.gz_run()
        self.set_output()
        self.end()
