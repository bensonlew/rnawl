# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow


class FileCheckWorkflow(Workflow):
    def __init__(self, wsheet_object):
        super(FileCheckWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "up", "type": "string"},
            {"name": "down", "type": "string"},
            {"name": "date", "type": "string"},
            {"name": "batch", "type": "string"},
            {"name": "flowcell", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.file_check = self.add_tool("medical.paternity_test_v2.file_check_V4")

    def run(self):
        # self.start_listener()
        # self.fire("start")
        options = {
            'up': self.option('up'),
            'down': self.option('down'),
            'date': self.option('date'),
            'batch': self.option('batch'),
            'flowcell': self.option("flowcell"),
            'main_id': self.option("main_id")
        }
        self.file_check.set_options(options)
        self.file_check.on("end", self.end)  # modified by hongdong 20180102 24-39 修改workflow不能正常运行结束
        self.file_check.run()
        # self.output_dir = self.file_check.output_dir
        super(FileCheckWorkflow, self).run()
        # self.end()

    def end(self):
        super(FileCheckWorkflow, self).end()
