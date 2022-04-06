# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class PfamAnnotAgent(Agent):
    '''
    last_modify: 2019.04.17
    '''
    def __init__(self, parent):
        super(PfamAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'pfam'},
            {'name': 'txpt_pfam_domain', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_pfam_domain', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
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
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        infile_size = os.path.getsize(self.option('pfam_domain').prop['path'])
        self._memory = '{}G'.format(int(float(infile_size) / 1024 ** 3 * 4 + 8))

    def end(self):
        super(PfamAnnotAgent, self).end()

class PfamAnnotTool(Tool):
    def __init__(self, config):
        super(PfamAnnotTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.pfam_gene_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/pfam_gene.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_venn_list.py')
        self.prefix = os.path.basename(self.option('pfam_domain').path)[:-4]
        self.txpt_tsv = os.path.join(self.output_dir, '{}.T.tsv'.format(self.option('database')))
        self.gene_tsv = os.path.join(self.output_dir, '{}.G.tsv'.format(self.option('database')))
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))

    def run(self):
        super(PfamAnnotTool, self).run()
        self.run_pfam_gene()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.set_output()
        self.end()

    def run_pfam_gene(self):
        shutil.copy(self.option('pfam_domain').path, self.txpt_tsv)
        print os.environ
        cmd = '{} {}'.format(self.python, self.pfam_gene_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --domain {}'.format(self.txpt_tsv)
        cmd += ' --output {}'.format(self.gene_tsv)
        cmd_name = 'run_pfam_gene'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.txpt_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.gene_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.gene_list)
        cmd_name = 'run_get_venn_list_g'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('txpt_pfam_domain').set_path(self.txpt_tsv)
        self.option('gene_pfam_domain').set_path(self.gene_tsv)
        self.option('txpt_list').set_path(self.txpt_list)
        self.option('gene_list').set_path(self.gene_list)
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
            'name': 'ref_rna_v2.annotation.pfam_annot',
            'instant': False,
            'options': {
                'pfam_domain': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotFilter/output/pfam/pfam.filter.xls',
                'longest_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotClass/AnnotFile/result/longest.t2g.tsv',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
