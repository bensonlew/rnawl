# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class PfamAnnotAgent(Agent):
    '''
    last_modify: 2019.02.14
    '''
    def __init__(self, parent):
        super(PfamAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'file_in', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'file_out', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {"name": "pfam_version", "type": "string", "default": "32"},
        ]
        self.add_option(options)
        self.step.add_steps('pfam_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pfam_annot.start()
        self.step.update()

    def step_end(self):
        self.step.pfam_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('file_in').is_set:
            self.logger.debug('{} = {}'.format('file_in', self.option('file_in').prop['path']))
            self.infile_size = os.path.getsize(self.option('file_in').prop['path'])
        else:
            raise OptionError('input PFAM DOMAIN must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 4 + 4))

    def end(self):
        super(PfamAnnotAgent, self).end()

class PfamAnnotTool(Tool):
    def __init__(self, config):
        super(PfamAnnotTool, self).__init__(config)

    def run(self):
        super(PfamAnnotTool, self).run()
        self.set_output()
        self.end()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.option('file_in').prop['path']
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.option('file_out', link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
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
            'id': 'pfam_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.pfam_annot',
            'instant': False,
            'options': {
                'file_in': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/pfam/pfam.filter.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
