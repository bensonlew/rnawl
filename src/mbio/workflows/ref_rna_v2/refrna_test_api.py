# -*- coding: utf-8 -*-

"""空workflow用于测试导表"""

from biocluster.workflow import Workflow

class RefrnaTestApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RefrnaTestApiWorkflow, self).__init__(wsheet_object)
        options = []
        self.task_id = self.sheet.id
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        """
        检查参数设置
        """

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()


    def set_output(self, event):
        obj = event["bind_object"]

    def run(self):
        super(RefrnaTestApiWorkflow, self).run()

    def end(self):
        super(RefrnaTestApiWorkflow, self).end()
