# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class DownloadAgent(Agent):
    '''
    last_modify: 2019.03.07
    '''
    def __init__(self, parent):
        super(DownloadAgent, self).__init__(parent)
        options = [
            {'name': 'file_in', 'type': 'infile', 'format':'lnc_rna.common'},
            {'name': 'file_name', 'type': 'string', 'default': None},
            {'name': 'file_out', 'type': 'outfile', 'format':'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('download')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.download.start()
        self.step.update()

    def step_end(self):
        self.step.download.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('file_in').is_set:
            self.logger.debug('{} = {}'.format('file_in', self.option('file_in').path))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

    def end(self):
        super(DownloadAgent, self).end()

class DownloadTool(Tool):
    def __init__(self, config):
        super(DownloadTool, self).__init__(config)

    def run(self):
        super(DownloadTool, self).run()
        self.set_output()
        self.end()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.option('file_in').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('file_out').set_path(link_name)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'download_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.download',
            'instant': False,
            'options': {
                'file_in': 's3://refrnav2/files/m_188/188_5bea79e58826b/tsg_33530/workflow_results/Sequence_database/refrna_seqs.db',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()