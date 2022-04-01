# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import unittest

class MirnaPrecursorModule(Module):
    '''
    last_modify: 2019.03.26
    '''
    def __init__(self, work_id):
        super(MirnaPrecursorModule, self).__init__(work_id)
        options = [
            {'name': 'lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'known_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 't2g', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(MirnaPrecursorModule, self).run()
        self.run_precursor_candidate()

    def run_precursor_candidate(self):
        self.step.add_steps('precursor_candidate')
        self.precursor_candidate = self.add_tool('whole_transcriptome.longrna.precursor_candidate')
        options = {
            'file_in': self.option('lncrna_fa'),
            'species': self.option('species')
        }
        self.precursor_candidate.set_options(options)
        self.precursor_candidate.on('start', self.set_step, {'start': self.step.precursor_candidate})
        self.precursor_candidate.on('end', self.set_step, {'end': self.step.precursor_candidate})
        self.precursor_candidate.on('end', self.run_precursor_stat)
        self.precursor_candidate.run()

    def run_precursor_stat(self):
        self.step.add_steps('precursor_stat')
        self.precursor_stat = self.add_tool('whole_transcriptome.longrna.precursor_stat')
        options = {
            'tblout': self.precursor_candidate.option('tblout'),
            'species': self.option('species'),
            'known_list': self.option('known_list'),
            'novel_list': self.option('novel_list'),
            't2g': self.option('t2g'),
        }
        self.precursor_stat.set_options(options)
        self.precursor_stat.on('start', self.set_step, {'start': self.step.precursor_stat})
        self.precursor_stat.on('end', self.set_step, {'end': self.step.precursor_stat})
        self.precursor_stat.on('end', self.set_output)
        self.precursor_stat.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.precursor_stat.option('tabular').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('tabular').set_path(link_name)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(MirnaPrecursorModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'whole_transcriptome.longrna.mirna_precursor_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'whole_transcriptome.longrna.mirna_precursor',
            'instant': False,
            'options': {
                'lncrna_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/all_lncrna.fa',
                'species': 'hsa',
                'known_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/known_lncrna_ids.list',
                'novel_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/novel_lncrna_ids.list',
                't2g': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/trans2gene',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()