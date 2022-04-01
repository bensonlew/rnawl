# -*- coding: utf-8 -*-
# __author__ = 'wuqin'
# from biocluster.agent import Agent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class BubbleWorkflow(Workflow):
    """
    气泡图小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BubbleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "x_axis", "type": "string"},
            {"name": "y_axis", "type": "string"},
            {"name": "size", "type": "string"},
            {"name": "group", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.bubble = self.add_tool("tool_lab.bubble")

    def check_options(self):
        """
        参数检查
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 2:
            raise OptionError('列数少于2，不可进行分析', code="32702904")
        if self.option('x_axis') not in self.option('tooltable').prop['col_sample']:
            raise OptionError('x轴数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('y_axis') not in self.option('tooltable').prop['col_sample']:
            raise OptionError('y轴数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('size'):
            if self.option('size') not in self.option('tooltable').prop['col_sample']:
                raise OptionError('气泡大小数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('group'):
            if self.option('group') not in self.option('tooltable').prop['col_sample']:
                raise OptionError('分组数据不在输入表格中，查看是否数据选择错误', code="32702909")
        return True

    def run_bubble(self):
        options = {
            "tooltable": self.option('tooltable'),
            "x_axis": self.option('x_axis'),
            "y_axis": self.option('y_axis'),
            "size": self.option('size'),
            "group": self.option('group')
        }
        self.bubble.set_options(options)
        self.bubble.on("end", self.set_output, "bubble")
        self.bubble.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bubble':
            self.linkdir(obj.output_dir, 'bubble')
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
        self.logger.info("开始导表")
        api_bubble = self.api.api("tool_lab.bubble")
        api_bubble.add_bubble_detail(self.option('main_id'), self.output_dir + '/bubble/bubble.xls')
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_bubble()
        super(BubbleWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "气泡图数据"],
            ["./bubble.xls", "xls", "气泡图数据"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BubbleWorkflow, self).end()
