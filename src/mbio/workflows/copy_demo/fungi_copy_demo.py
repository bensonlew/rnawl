# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.workflow import Workflow
from mbio.packages.fungi_genome.copy_demo import CopyDemo
from biocluster.core.exceptions import OptionError

class FungiCopyDemoWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FungiCopyDemoWorkflow, self).__init__(wsheet_object)
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
        if not self.option("task_id"):
            raise OptionError("缺少demo的task_id,请检查")
        if not self.option("target_task_id"):
            raise OptionError("缺少demo的新task_id,请检查")
        if not self.option("target_project_sn"):
            raise OptionError("缺少demo的新project_sn,请检查")
        if not self.option("target_member_id"):
            raise OptionError("缺少demo的新member_id,请检查")

    def run(self):
        self.start_listener()
        self.fire("start")
        copy_task = CopyDemo(self.option("task_id"), self.option("target_task_id"), self.option("target_project_sn"),
                             self.option("target_member_id"), self.option("target_member_type"))
        copy_task.run()
        self.end()

    def end(self):
        super(FungiCopyDemoWorkflow, self).end()