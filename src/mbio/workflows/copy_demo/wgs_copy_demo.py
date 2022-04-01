# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0508

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.wgs.wgs_copy_demo import WgsCopyDemo


class WgsCopyDemoWorkflow(Workflow):
    """
    拉取WGS demo数据
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgsCopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
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
        copy_task = WgsCopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),\
                    new_project_sn=self.option("target_project_sn"), new_member_id=self.option("target_member_id"))
        copy_task.run()
        self.end()

    def end(self):
        super(WgsCopyDemoWorkflow, self).end()
