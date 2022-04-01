# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class DisgenetAnnotAgent(Agent):
    '''
    last_modify: 2019.04.16
    '''
    def __init__(self, parent):
        super(DisgenetAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'disgenet_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'disgenet'},
            {'name': 'disgenet_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "disgenet_version", "type": "string", "default": "72"},
        ]
        self.add_option(options)
        self.step.add_steps('cog_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cog_annot.start()
        self.step.update()

    def step_end(self):
        self.step.cog_annot.finish()
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
        super(DisgenetAnnotAgent, self).end()

class DisgenetAnnotTool(Tool):
    def __init__(self, config):
        super(DisgenetAnnotTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        # self.xml2table_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/xml2table.py')
        self.txml2gxml_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/txml2gxml.py')
        self.disgenet_transition_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/disgenet_transition.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/get_venn_list.py')
        self.prefix = "disgenet_annot"
        self.txpt_xml = os.path.join(self.work_dir, '{}.T.xml'.format(self.prefix))
        self.txpt_tsv = os.path.join(self.work_dir, '{}.T.tsv'.format(self.prefix))
        self.disgenet_tsv = os.path.join(self.output_dir, '{}.G.tsv'.format(self.option('database')))
        self.disgenet_detail = os.path.join(self.output_dir, '{}.G.detail.tsv'.format(self.option('database')))

        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))

    def run(self):
        super(DisgenetAnnotTool, self).run()
        # self.run_xml2table()
        self.run_disgenet_transition()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.set_output()
        self.end()

    def run_disgenet_transition(self):
        cmd = '{} {}'.format(self.python, self.disgenet_transition_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --disgenet_id {}'.format(self.option("disgenet_ids").path)
        cmd += ' --output {}'.format(self.disgenet_tsv)
        cmd_name = 'run_disgenet_genes'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.option("disgenet_ids").path)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.disgenet_tsv)
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
                        self.set_error('fail to run %s, abord', variables=(n), code="33710702")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('disgenet_table').set_path(self.disgenet_tsv)
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
            'id': 'disgenet_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.disgenet_annot',
            'instant': False,
            'options': {
                'disgenet_ids': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/disgenet_annot.tsv',
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
