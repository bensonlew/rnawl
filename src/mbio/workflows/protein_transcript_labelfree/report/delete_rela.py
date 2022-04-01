# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mbio.packages.protein_transcript_labelfree.delete_rela import DeleteRelaMongo


class DeleteRelaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeleteRelaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "project_type", "type": 'string', "default": ''},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire("start")
        delete_demo = DeleteRelaMongo(self.option("task_id"), self.option("project_type"))
        delete_demo.run()
        self.end()

    def end(self):
        super(DeleteRelaWorkflow, self).end()
