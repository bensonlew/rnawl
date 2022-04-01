
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.file import download


class DownloadAgent(Agent):
    """
    author: guanqing.zou
    last_modify: 20181017
    """
    def __init__(self, parent):
        super(DownloadAgent, self).__init__(parent)
        options = [
            {"name": "item","type": "string","default": ""}  # "[['s3_path",'./file_name'],]"  如果不下载，则 '[]'
            ]
        self.add_option(options)


    def check_options(self):
        if not self.option("item"):
            raise OptionError("必须设置参数", code="32100601")


    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        super(DownloadAgent, self).end()


class DownloadTool(Tool):
    def __init__(self, config):
        super(DownloadTool, self).__init__(config)
        self.list = eval(self.option('item'))

    def run_download(self):
        if self.list != '[]':
            self.logger.info("开始下载")
            for i,j in self.list:
                download(i,j)
                #self.download_from_s3(i,j)
            self.logger.info("完成下载")


    def set_output(self):
        pass

    def run(self):
        """
        运行
        :return:
        """
        super(DownloadTool, self).run()
        self.run_download()
        self.end()


