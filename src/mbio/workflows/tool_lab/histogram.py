# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import subprocess
import re


class HistogramWorkflow(Workflow):
    """
    频率直方图工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HistogramWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_table", "type": "infile", "format": "tool_lab.histogram_table"},
            {"name": "fq", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.histogram = self.add_tool("tool_lab.histogram")

    def check_options(self):
        """
        检查参数
        """
        if not self.option("input_table"):
            raise OptionError('必须设置作图数据')
        if not self.option("fq"):
            raise OptionError('必须设置频率')

    def run_histogram(self):
        options = {
            "input_table": self.option('input_table'),
            "fq": self.option('fq'),
        }
        self.histogram.set_options(options)
        self.histogram.on("end", self.set_output, "histogram")
        self.histogram.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'histogram':
            self.linkdir(obj.output_dir, 'histogram')
        self.set_db()

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

    def set_db(self):
        """
        保存结果到Mongo库
        """
        self.logger.info("保存结果到Mongo")
        main_id = self.option("main_id")
        histogram_api = self.api.api("tool_lab.histogram_api")
        histogram_table = self.output_dir + '/histogram/histogram.csv'
        histogram_api.add_histogram_detail(main_id, histogram_table)
        self.logger.info("保存结果到Mongo结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_histogram()
        super(HistogramWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "频率直方图数据"],
            ["./histogram.csv", "csv", "频率直方图数据"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(HistogramWorkflow, self).end()
