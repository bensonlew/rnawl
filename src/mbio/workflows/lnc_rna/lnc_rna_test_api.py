# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow

class LncRnaTestApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LncRnaTestApiWorkflow, self).__init__(wsheet_object)
        options = list()
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        pass

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_output(self, event):
        obj = event['bind_object']

    def run(self):
        super(LncRnaTestApiWorkflow, self).run()

    def end(self):
        super(LncRnaTestApiWorkflow, self).end()
