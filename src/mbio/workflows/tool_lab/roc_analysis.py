# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200706

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class RocAnalysisWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RocAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "roc_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "smooth", "type": "bool", "default": False},
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.roc_analysis = self.add_tool("tool_lab.roc_analysis")

    def check_options(self):
        if not self.option("roc_table"):
            raise OptionError("必须输入roc_table文件")
        return True

    def run_regression(self):
        options = {
            "roc_table": self.option('roc_table'),
            "main_id": self.option('main_id'),
            "smooth": self.option('smooth'),
        }
        self.roc_analysis.set_options(options)
        self.roc_analysis.on("end", self.set_output, "roc_analysis")
        self.roc_analysis.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'roc_analysis':
            self.linkdir(obj.output_dir, 'roc_analysis')
        self.end()

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

    def run(self):
        self.run_regression()
        super(RocAnalysisWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(RocAnalysisWorkflow, self).end()
