# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
from bson.objectid import ObjectId


class DeleteDemoV2Workflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeleteDemoV2Workflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "project_type", "type": 'string', "default": ''},
            {"name": "db_version", "type": 'int', "default": ''},
            {"name": "main_id", "type": 'string', "default": ''},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.delete_demo = self.api.api("project_demo.delete_demo")

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire("start")
        delete_demo = DeleteDemoMongo(self.option("task_id"), self.option("project_type"), db_version=int(self.option("db_version")))
        main_id = ObjectId(self.option("main_id"))
        try:
            delete_demo.run()
        except:
            self.delete_demo.update_main_table('sg_task_delete', main_id, 'failed')
        else:
            self.delete_demo.update_main_table('sg_task_delete', main_id, 'end')
        self.end()

    def end(self):
        super(DeleteDemoV2Workflow, self).end()
