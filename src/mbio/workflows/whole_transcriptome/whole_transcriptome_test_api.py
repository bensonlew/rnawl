# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow


class WholeTranscriptomeTestApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WholeTranscriptomeTestApiWorkflow, self).__init__(wsheet_object)
        self.add_option(list())
        self.set_options(self._sheet.options())

    def check_options(self):
        for name in dir(self.sheet):
            pass
            # if name in ['__delattr__', '__getattribute__', '__hash__', '__repr__', '__setattr__', '__str__']:
            #     continue
            # self.logger.info('{} -> {}'.format(name, getattr(self._sheet, name)))
        else:
            self.logger.debug(self.sheet.data)
            self.logger.debug('Whole Transcriptome Test Api Workflow')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_output(self, event):
        obj = event['bind_object']

    def run(self):
        super(WholeTranscriptomeTestApiWorkflow, self).run()

    def end(self):
        super(WholeTranscriptomeTestApiWorkflow, self).end()
