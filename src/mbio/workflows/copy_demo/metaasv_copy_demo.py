# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mbio.packages.metaasv.copy_demo import CopyDemo


class MetaasvCopyDemoWorkflow(Workflow):
    """
    metaasv demo的拷贝
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self._project_type = "metaasv"
        super(MetaasvCopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''},
            {"name": "project_type", "type": 'string', "default": ''},
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
        super(MetaasvCopyDemoWorkflow, self).end()
