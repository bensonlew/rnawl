# -*- coding: utf-8 -*-
# __author__ = 'wuqin'
# from biocluster.agent import Agent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class ViolinWorkflow(Workflow):
    """
    小提琴图小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ViolinWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "grouptable", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "log", "type": "string", "default": "none"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.violin = self.add_tool("tool_lab.violin")

    def check_options(self):
        """
        参数检查
        """
        if not self.option('tooltable'):
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('log'):
            if self.option('log') not in ['log2', 'log10', 'no']:
                raise OptionError('数据处理方式不存在', code="32702909")
        if self.option('tooltable').prop['sample_num'] < 2:
            raise OptionError('列数少于2，不可进行分析', code="32702904")
        return True

    def run_violin(self):
        if self.option('log') == "no":
            log = "none"
        else:
            log = self.option('log')
        options = {
            "tooltable": self.option('tooltable'),
            "grouptable": self.option('grouptable'),
            "log": log
        }
        self.violin.set_options(options)
        self.violin.on("end", self.set_output, "violin")
        self.violin.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'violin':
            self.linkdir(obj.output_dir, 'violin')
        self.set_db()
        # self.end()

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
        api_violin = self.api.api("tool_lab.violin")
        api_violin.add_violin_detail(self.option('main_id'), self.output_dir + '/violin')
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_violin()
        super(ViolinWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option('grouptable').is_set:
            result_dir.add_relpath_rules([
                [".", "", "小提琴图数据"],
                ["./violin_sum.xls", "xls", "小提琴图求和数据"],
                ["./violin_mean.xls", "xls", "小提琴图平均数数据"],
                ["./violin_median.xls", "xls", "小提琴图中位数数据"],
                ["./violin.xls", "xls", "小提琴图画图数据"],
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "小提琴图数据"],
                ["./violin.xls", "xls", "小提琴图画图数据"]
            ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(ViolinWorkflow, self).end()
