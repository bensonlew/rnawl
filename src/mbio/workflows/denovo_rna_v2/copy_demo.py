# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mbio.packages.denovo_rna_v2.copy_demo import CopyDemoMongo


class CopyDemoWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''},
            {"name": "project_type", "type": 'string', "default": ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire("start")
        copy_demo = CopyDemoMongo(self.option("task_id"),
                                  self.option("target_task_id"),
                                  self.option("target_project_sn"),
                                  self.option("target_member_id"),
                                  self.option("project_type"),)
        copy_demo.run()
        self.end()

    def end(self):
        # self.step.finish()
        # self.step.update()
        # self.logger.info("运行结束!")
        # self._update("end")
        # self.set_end()
        # self.fire('end')
        super(CopyDemoWorkflow, self).end()
