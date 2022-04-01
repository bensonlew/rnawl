# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.functions import toolfuncdeco


class DownloadAgent(Agent):
    """
    last_modify: 2019.06.13
    """

    def __init__(self, parent):
        super(DownloadAgent, self).__init__(parent)
        options = [
            {'name': 'ifile', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'basename', 'type': 'string', 'default': None},
            {'name': 'ofile', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
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

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

    @toolfuncdeco
    def end(self):
        super(DownloadAgent, self).end()


class DownloadTool(Tool):
    def __init__(self, config):
        super(DownloadTool, self).__init__(config)

    @toolfuncdeco
    def run(self):
        super(DownloadTool, self).run()
        self.set_output()
        self.end()

    @toolfuncdeco
    def set_output(self):
        source = self.option('ifile').path
        if self.option('basename'):
            link_name = os.path.join(self.output_dir, self.option('basename'))
        else:
            link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('ofile').set_path(link_name)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'download_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.download',
            'instant': False,
            'options': {
                'ifile': 's3nb://refrnav2/files/m_34394/34394_5d80a9c6283a7/i-sanger_204579/intermediate_results'
                         '/Align/AlignBam/',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
