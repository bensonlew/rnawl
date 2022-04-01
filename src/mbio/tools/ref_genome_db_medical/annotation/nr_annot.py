# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng, liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class NrAnnotAgent(Agent):
    '''
    last_modify: 2019.04.22
    '''
    def __init__(self, parent):
        super(NrAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'nr_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'nr'},
            {'name': 'txpt_table', 'type': 'outfile', 'format': 'ref_rna_v2.blast_table'},
            {'name': 'gene_table', 'type': 'outfile', 'format': 'ref_rna_v2.blast_table'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('nr_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.nr_annot.start()
        self.step.update()

    def step_end(self):
        self.step.nr_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(NrAnnotAgent, self).end()

class NrAnnotTool(Tool):
    def __init__(self, config):
        super(NrAnnotTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        # self.xml2table_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/xml2table.py')

        # self.txml2gxml_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/txml2gxml.py')
        self.tids2gids = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/tids2gids.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/get_venn_list.py')
        self.prefix = "annot_nr"
        self.txpt_ids = os.path.join(self.work_dir, '{}.T.ids'.format(self.prefix))
        self.gene_ids = os.path.join(self.work_dir, '{}.G.ids'.format(self.prefix))
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))

    def run(self):
        super(NrAnnotTool, self).run()
        # self.run_txml2gxml()
        # self.run_xml2table_t()
        self.run_idst2g()
        # self.run_xml2table_g()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.set_output()
        self.end()

    def run_idst2g(self):
        cmd = '{} {}'.format(self.python, self.tids2gids)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --tids {}'.format(self.option('nr_ids').path)
        cmd += ' --output_t {}'.format(self.txpt_ids)
        cmd += ' --output_g {}'.format(self.gene_ids)
        cmd += ' --output_format {}'.format("{nr}({description})")
        cmd_name = 'run_ids2g'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.txpt_ids)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.gene_ids)
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
                        self.set_error('fail to run %s, abord', variables=(n), code="33711202")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))

        self.option('txpt_table').set_path(self.txpt_ids)
        self.option('gene_table').set_path(self.gene_ids)
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
            'id': 'nr_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.nr_annot',
            'instant': False,
            'options': {
                'nr_ids': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/nr_annot.tsv',
                'longest_t2g': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/begin/longest.t2g.tsv',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
