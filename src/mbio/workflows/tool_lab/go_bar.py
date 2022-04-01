# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200706

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class GoBarWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GoBarWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "go_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.go_bar = self.add_tool("tool_lab.go_bar")

    def check_options(self):
        if not self.option("go_table"):
            raise OptionError("必须输入snp_table文件")
        return True

    def run_regression(self):
        options = {
            "go_table": self.option('go_table'),
            "main_id": self.option('main_id'),
        }
        self.go_bar.set_options(options)
        self.go_bar.on("end", self.set_output, "go_bar")
        self.go_bar.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'go_bar':
            self.linkdir(obj.output_dir, 'go_bar')
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
        super(GoBarWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GoBarWorkflow, self).end()
