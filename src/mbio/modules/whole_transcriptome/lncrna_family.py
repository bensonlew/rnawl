# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import glob
import os
import unittest

class LncrnaFamilyModule(Module):
    '''
    last_modify: 2019.04.19
    '''
    def __init__(self, work_id):
        super(LncrnaFamilyModule, self).__init__(work_id)
        options = [
            {'name': 'lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'number', 'type': 'int', 'default': 10},
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
        super(LncrnaFamilyModule, self).run()
        self.run_family_prep()

    def run_family_prep(self):
        self.step.add_steps('family_prep')
        self.family_prep = self.add_tool('lnc_rna.structure.family_prep')
        options = {
            'lncrna_fa': self.option('lncrna_fa'),
            'number': self.option('number')
        }
        self.family_prep.set_options(options)
        self.family_prep.on('start', self.set_step, {'start': self.step.family_prep})
        self.family_prep.on('end', self.set_step, {'end': self.step.family_prep})
        self.family_prep.on('end', self.run_cmscan)
        self.family_prep.run()

    def run_cmscan(self):
        for n, seqfile in enumerate(glob.glob(os.path.join(self.family_prep.output_dir, '*'))):
            self.step.add_steps('cmscan_{}'.format(n))
            cmscan = self.add_tool('lnc_rna.structure.cmscan')
            cmscan.set_options({'seqfile': seqfile})
            cmscan.on('start', self.set_step, {'start': getattr(self.step, 'cmscan_{}'.format(n))})
            cmscan.on('end', self.set_step, {'start': getattr(self.step, 'cmscan_{}'.format(n))})
            self.tools.append(cmscan)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_family_stat)
        else:
            self.on_rely(self.tools, self.run_family_stat)
        for tool in self.tools:
            tool.run()

    def run_family_stat(self):
        tsv_list = os.path.join(self.work_dir, 'tsv.list')
        open(tsv_list, 'w').writelines('{}\n'.format(tool.option('tsv').path) for tool in self.tools)
        self.step.add_steps('family_stat')
        self.family_stat = self.add_tool('lnc_rna.structure.family_stat')
        options = {
            'tsv_list': tsv_list,
            'known_list': self.option('known_list'),
            'novel_list': self.option('novel_list'),
            't2g': self.option('t2g'),
        }
        self.family_stat.set_options(options)
        self.family_stat.on('start', self.set_step, {'start': self.step.family_stat})
        self.family_stat.on('end', self.set_step, {'end': self.step.family_stat})
        self.family_stat.on('end', self.set_output)
        self.family_stat.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.family_stat.option('tabular').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('tabular').set_path(link_name)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(LncrnaFamilyModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'lncrna_family_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.lncrna_family',
            'instant': False,
            'options': {
                'lncrna_fa': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33834/FilterByExpress/output/filtered_file/all_lncrna.fa',
                'known_list': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33834/FilterByExpress/output/filtered_file/known_lncrna_ids.list',
                'novel_list': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33834/FilterByExpress/output/filtered_file/novel_lncrna_ids.list',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33834/Assemble/output/NewTranscripts/trans2gene',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
