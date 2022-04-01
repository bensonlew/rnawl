# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# last modified add target_member_type @ 20171201
from biocluster.workflow import Workflow


class MetagenomicCopyDemoWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetagenomicCopyDemoWorkflow, self).__init__(wsheet_object)
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
        from mbio.packages.metagenomic.copy_demo import CopyDemo
        copy_task = CopyDemo(self.option("task_id"), self.option("target_task_id"), self.option("target_project_sn"),
                             self.option("target_member_id"), self.option("target_member_type"))  # 不要db参数,增加target_member_type
        copy_task.run()
        self.end()

    def end(self):
        super(MetagenomicCopyDemoWorkflow, self).end()
