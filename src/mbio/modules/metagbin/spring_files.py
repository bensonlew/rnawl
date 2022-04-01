# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re,shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class SpringFilesModule(Module):
    """
    多个files进行解压或压缩
    author: gaohao
    last_modify: 2020.04.20
    """
    def __init__(self, work_id):
        super(SpringFilesModule, self).__init__(work_id)
        options = [
            {"name": "file_dir", "type": "string"}, #需要压缩的目录或需要解压的目录
            {"name": "type", "type": "string"},  # compress or decompress
        ]
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('type'):
            raise OptionError('必须输入type!')
        else:
            if self.option('type') not in ["compress", "decompress"]:
                raise OptionError('必须输入正确的{}！compress or decompress'.format(self.option('type')))

    def run_amphora(self):
        if self.option('type') in ["compress"]:
            files = os.listdir(self.option("file_dir"))
            for k in files:
                if not re.search("list.txt", k):
                    self.spring = self.add_tool('metagbin.spring_file')
                    opts = {
                        "input": self.option("file_dir") + '/' + k,
                        "type": self.option("type"),
                    }
                    self.spring.set_options(opts)
                    self.modules.append(self.spring)
            self.logger.info(self.modules)
        elif self.option('type') in ["decompress"]:
            files = os.listdir(self.option("file_dir"))
            for k in files:
                if not re.search("list.txt", k):
                    self.spring = self.add_tool('metagbin.spring_file')
                    opts = {
                        "input": self.option("file_dir") + '/' + k,
                        "type": self.option("type"),
                    }
                    self.spring.set_options(opts)
                    self.modules.append(self.spring)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(SpringFilesModule, self).run()
        self.run_amphora()


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            link_dir(module.output_dir, self.output_dir)
        if os.path.exists(self.output_dir + "/list.txt"):
            os.remove(self.output_dir + "/list.txt")
        os.link(self.option("file_dir") + "/list.txt", self.output_dir + "/list.txt")
        self.end()

    def end(self):
        super(SpringFilesModule, self).end()