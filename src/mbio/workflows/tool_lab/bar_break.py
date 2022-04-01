# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200416

import os
import re
import math
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class BarBreakWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BarBreakWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "set_group", "type": "string", "default": True},
            {"name": "ishape", "type": "string", "default": "sd"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "low_point", "type": "float"},  # 下断点值
            {"name": "high_point", "type": "float"},  # 上断点值
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.bar_break = self.add_tool("tool_lab.bar_break")

    def check_options(self):
        if not self.option("bar_table"):
            raise OptionError("必须输入snp_table文件")
        if not self.option("low_point"):
            raise OptionError("必须输入下断点值")
        if not self.option("high_point"):
            raise OptionError("必须输入上断点值")
        if not self.option("ishape"):
            raise OptionError("必须输入ishape取值")
        if not self.option("set_group"):
            raise OptionError("必须输入set_group取值")

    def run_bar_break(self):
        options = {
            "bar_table": self.option('bar_table'),
            "group_table": self.option('group_table'),
            "low_point": self.option('low_point'),
            'high_point': self.option('high_point'),
            "main_id": self.option('main_id'),
            "ishape":self.option('ishape'),
        }
        self.bar_break.set_options(options)
        self.bar_break.on("end", self.set_output, "bar_break")
        self.bar_break.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bar_break':
            self.linkdir(obj.output_dir, 'bar_break')

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))
        time.sleep(1)
        self.end()

    def run(self):
        self.run_bar_break()
        super(BarBreakWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BarBreakWorkflow, self).end()
