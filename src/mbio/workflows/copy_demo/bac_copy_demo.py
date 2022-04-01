# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.workflow import Workflow
from mbio.packages.bacgenome.copy_demo import CopyDemo

class BacCopyDemoWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BacCopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''},
            {"name": "target_member_type", "type": 'int', "default": ''},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire("start")
        copy_task = CopyDemo(self.option("task_id"), self.option("target_task_id"), self.option("target_project_sn"),
                             self.option("target_member_id"), self.option("target_member_type"))
        copy_task.run()
        self.end()

    def end(self):
        super(BacCopyDemoWorkflow, self).end()