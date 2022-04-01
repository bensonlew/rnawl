# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# 测试childend的问题
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import time


class TestQiuAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(TestQiuAgent, self).__init__(parent)
        options = [
            {"name": "test", "type": "string", "default": 'qiuqiu'}]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        """
        self._cpu = 10
        self._memory = ''

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['test.txt', '', '测试文件']
        ])
        super(TestQiuAgent, self).end()


class TestQiuTool(Tool):
    def __init__(self, config):
        super(TestQiuTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(TestQiuTool, self).run()
        with open(self.output_dir + '/test.txt', 'w') as w:
            w.write('this is test info!!!')
        for i in range(3):
            time.sleep(3)
            self.logger.info('test tool running')
        self.end()
