# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

import os

from biocluster.workflow import Workflow

from mbio.packages.project_demo.delete_records import DeleteRecordsMongo


class DeleteRecordsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeleteRecordsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'task_id', 'type': 'string', 'default': ''},
            {'name': 'project_type', 'type': 'string', 'default': ''},
            {'name': 'submit_location', 'type': 'string', 'default': ''},
            {'name': 'status', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire('start')
        self.logger.debug('start initializing the instance of class <DeleteRecordsMongo>')
        delete_demo = DeleteRecordsMongo(self.option('task_id'),
                                         self.option('project_type'),
                                         self.option('submit_location').lower(),
                                         self.option('status').lower())
        delete_demo.logger = self.logger
        delete_demo.run()
        self.logger.debug('succeed in calling <run> method of class <DeleteRecordsMongo>')
        self.end()

    def end(self):
        # open(os.path.join(self.output_dir, 'done')).close()
        super(DeleteRecordsWorkflow, self).end()
