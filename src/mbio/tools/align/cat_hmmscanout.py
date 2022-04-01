# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest


class CatHmmscanoutAgent(Agent):
    """
    CatBlastout:往同一文件中合并分开的结果
    version 1.0
    author: zhouxuan
    last_modify: 20170524
    """

    def __init__(self, parent):
        super(CatHmmscanoutAgent, self).__init__(parent)
        options = [
            {"name": "hmmscan_out", "type": "infile", "format": "paternity_test.data_dir"},
            {"name": "hmmscan_result", "type": "outfile", "format": "paternity_test.tab"},
        ]
        self.add_option(options)
        self.step.add_steps('cat_hmmscan')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cat_hmmscan.start()
        self.step.update()

    def step_end(self):
        self.step.cat_hmmscan.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("hmmscan_out").is_set:
            raise OptionError("请传入hmmscan输出结果文件", code="31101001")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""],
        ])
        super(CatHmmscanoutAgent, self).end()

class CatHmmscanoutTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatHmmscanoutTool, self).__init__(config)
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def cat_files(self):
        """
        合并结果文件
        """
        file_name = os.listdir(self.option('hmmscan_out').prop['path'])
        self.result_name = os.path.join(self.output_dir, "hmmscan_result")
        command_prefix = "cat"
        command_times = 0
        for f in file_name:
            command_times += 1
            file_path = os.path.join(self.option('hmmscan_out').prop['path'], f)
            cmd = '{} {} {}'.format(self.sh_path, file_path, self.result_name)
            self.logger.info("start cat {}".format(f))
            command_name = command_prefix + "_%s" % command_times
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(f))
            else:
                self.set_error("cat %s error", variables=(f), code="31101001")
                # raise Exception("cat {} error".format(f))

    def run(self):
        super(CatHmmscanoutTool, self).run()
        self.cat_files()
        self.set_output()
        self.end()

    def set_output(self):
        try:
            self.option("hmmscan_result", self.result_name)
        except Exception as e:
            self.set_error("结果设置错误: %s", variables=(e), code="31101002")
            raise Exception("结果设置错误{}".format(e))