# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# 何胜用于测试不使用rpc服务的tool
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import datetime
import time
import pickle




####


class TestAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(TestAgent, self).__init__(parent)
        options = [
            {"name": "test", "type": "string", "default": 'shenghe'},]
        self.add_option(options)
        self._run_mode = 'local'

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
        print self.get_upload_files()
        super(TestAgent, self).end()

    def run(self):
        super(Agent, self).run()
        config_file = self.save_config()
        # self.job = self._job_manager.add_job(self)
        with open(config_file, "r") as f:
            config = pickle.load(f)
            config.DEBUG = True  # runtool设置了这个值
        self._start_run_time = datetime.datetime.now()
        tool = TestTool(config)
        tool.run()
        self.finish_callback()

    def finish_callback(self):
        """
        收到远程发送回的 :py:class:`biocluster.core.actor.State` end状态时的处理函数，设置当前Agent状态为结束

        :return:
        """
        self.load_output()
        self._status = "E"
        self._end_run_time = datetime.datetime.now()
        secends = (self._end_run_time - self._start_run_time).seconds
        self.logger.info("任务运行结束，运行时间:%ss" % secends)
        # self.job.set_end()
        self.end()


class TestTool(Tool):
    def __init__(self, config):
        super(TestTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(TestTool, self).run()
        with open(self.output_dir + '/test.txt', 'w') as w:
            w.write('this is test info!!!')
        for i in range(3):
            time.sleep(1)
            self.logger.info('test tool running')
        self.end()
