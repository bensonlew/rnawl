# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
# from biocluster.config import config
from collections import namedtuple, defaultdict
import unittest
import datetime
import subprocess
import random
import re
import os
import sys
import shutil

class ZipresultAgent(Agent):
    def __init__(self, parent):
        super(ZipresultAgent, self).__init__(parent)
        options = [
            {"name":"result_file","type":"string"},
            {"name":"output_name","type":'string'},
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen
    
    def check_option(self):
        '''
        参数检查
        '''
        if not self.option("result_file"):
            raise OptionError("没有输入文件")

    def set_resource(self):
        self._cpu = 2 
        self._memory = '20G'

    def end(self):
        super(ZipresultAgent, self).end()


class ZipresultTool(Tool):
    def __init__(self, config):
        super(ZipresultTool, self).__init__(config)
        self._version = "v.10"
        self.software_dir = self.config.SOFTWARE_DIR
        self.result_file = self.option('result_file')
        self.output_name = self.option('output_name')

    def run(self):
        super(ZipresultTool, self).run()
        self.run_zip()
        self.set_output()
        self.end()

    def run_zip(self,):
        name = os.path.basename(self.result_file)
        self.path = os.path.abspath(self.result_file).split(name)[0]
        os.chdir(self.path)
        cmd = "zip -q -r '{}.zip' '{}'".format(self.output_name,name)
        # command = self.add_command("zip", cmd, ignore_error = True)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("运行压缩结果文件完成")
        # else:
        #     self.set_error("运行压缩结果文件完成")
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "zipresult_{}.sh".format(now_time)
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/tool_lab/Yoogene/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成结果压缩文件")
        self.logger.info(shell)
        command1 =self.add_command("zipresult", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行zip完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if len(self.output_dir)  > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.path,"{}.zip".format(self.output_name)),os.path.join(self.output_dir, "{}.zip".format(self.output_name)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "zipresult_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.zipresult",
            "options": dict(
                result_file="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_yfullresult_353/Yfullresult/output",
  # tree="/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/20200801_tree_fix.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
    
if __name__ == '__main__':
    unittest.main()
