# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20180612

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class CoverageWindowModule(Module):
    """
     用于coverage_window 接口的module，用于分配各个tool。
    """
    def __init__(self, work_id):
        super(CoverageWindowModule, self).__init__(work_id)
        options = [
            {"name": "bam_dir", "type": 'infile', "format": "dna_evolution.bam_dir"},  # 这里depth_dir虽然写的是depth，其实是包含bam的文件夹，本module的第一个tool就是把bam转化为第二个tool可以使用的文件
            {"name": "step_num", "type": "int", "default": 200},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "project_sn", "type": "string"}
        ]
        self.add_option(options)
        self.coverage_window_tools = []
        self.samtool_depth_tools = []

    def check_options(self):
        if not self.option("bam_dir"):
            raise OptionError("请输入深度信息文件路径", code="24901001")
        if not self.option("step_num") in [5, 10, 50, 100, 200, 500]:
            raise OptionError("step_num 只能是5K，10K，50K,100K，200K，500K", code="24901002")
        return True

    def coverage_window_run(self):
        for f in os.listdir(self.option("bam_dir").prop["path"]):
            if re.match(".*bai$", f):   # 筛选去掉bai文件
                pass
            else:
                depth_file = os.path.join(self.option("bam_dir").prop["path"], f)
                coverage_window = self.add_tool("dna_evolution.coverage_window")
                coverage_window.set_options({
                    "bam_file": depth_file,
                    "step_num": self.option("step_num")
                })
                self.coverage_window_tools.append(coverage_window)
        for i in range(len(self.coverage_window_tools)):
            self.coverage_window_tools[i].on("end", self.set_output, 'coverage_window_dir')
        if self.coverage_window_tools:
            if len(self.coverage_window_tools) > 1:
                self.on_rely(self.coverage_window_tools, self.set_db)
            elif len(self.coverage_window_tools) == 1:
                self.coverage_window_tools[0].on('end', self.set_db)
        else:
            raise Exception("coverage_window_tools列表为空！")
        for tool in self.coverage_window_tools:
            gevent.sleep(1)
            tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'coverage_window_dir':
            self.linkdir(obj.output_dir, 'coverage_window_dir')
        else:
            pass

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

        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        """
        保存结果到mongo
        """
        if self.option("main_id"):
            api_coverage = self.api.api("dna_evolution.coverage_window")
            api_coverage.get_sample(self.output_dir + "/coverage_window_dir", self.option("main_id"))
        self.end()

    def run(self):
        super(CoverageWindowModule, self).run()
        self.coverage_window_run()

    def end(self):
        super(CoverageWindowModule, self).end()
