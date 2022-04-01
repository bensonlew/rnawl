# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class EggnogAnnotAgent(Agent):
    '''
    last_modify: 2019.04.16
    '''
    def __init__(self, parent):
        super(EggnogAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'cog'},
            {'name': 'cog_table', 'type': 'outfile', 'format': 'ref_rna_v2.cog_table'},
            {'name': 'cog_gene_table', 'type': 'outfile', 'format': 'ref_rna_v2.cog_table'},
            {'name': 'txpt_summary', 'type': 'outfile', 'format': 'ref_rna_v2.cog_summary'},
            {'name': 'gene_summary', 'type': 'outfile', 'format': 'ref_rna_v2.cog_summary'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'eggnog_version', 'type': 'string', 'default': None},

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
        infile_size = os.path.getsize(self.option('blast_xml').prop['path'])
        self._memory = '{}G'.format(int(float(infile_size) / 1024 ** 3 * 2 + 4))

    def end(self):
        super(EggnogAnnotAgent, self).end()

class EggnogAnnotTool(Tool):
    def __init__(self, config):
        super(EggnogAnnotTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.xml2table_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/xml2table.py')
        self.txml2gxml_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/txml2gxml.py')
        self.cog_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/cog_annotation.py')
        self.cog_transition_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/cog_transition.py')
        self.cog_summary_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/cog_summary.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_venn_list.py')
        self.prefix = os.path.basename(self.option('blast_xml').path)[:-4]
        self.txpt_xml = os.path.join(self.work_dir, '{}.T.xml'.format(self.prefix))
        self.txpt_tsv = os.path.join(self.work_dir, '{}.T.tsv'.format(self.prefix))
        self.cog_tsv = os.path.join(self.output_dir, '{}.tsv'.format(self.option('database')))
        self.cog_gene_tsv = os.path.join(self.output_dir, '{}.G.tsv'.format(self.option('database')))
        self.txpt_summary_tsv = os.path.join(self.output_dir, 'summary.T.tsv')
        self.gene_summary_tsv = os.path.join(self.output_dir, 'summary.G.tsv')
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))

    def run(self):
        super(EggnogAnnotTool, self).run()
        self.run_xml2table()
        self.run_cog_annotation()
        self.run_cog_transition()
        self.run_cog_summary()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.set_output()
        self.end()

    def run_xml2table(self):
        shutil.copy(self.option('blast_xml').path, self.txpt_xml)
        cmd = '{} {} -a'.format(self.python, self.xml2table_py)
        cmd += ' -i {}'.format(self.txpt_xml)
        cmd += ' -o {}'.format(self.txpt_tsv)
        cmd_name = 'run_xml2table'
        self.run_code(cmd_name, cmd)

    def run_cog_annotation(self):
        cmd = '{} {}'.format(self.python, self.cog_annotation_py)
        cmd += ' -i {}'.format(self.txpt_tsv)
        cmd += ' -o {}'.format(self.cog_tsv)
        cmd_name = 'run_cog_annotation'
        self.run_code(cmd_name, cmd)

    def run_cog_transition(self):
        cmd = '{} {}'.format(self.python, self.cog_transition_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --cog {}'.format(self.cog_tsv)
        cmd += ' --output {}'.format(self.cog_gene_tsv)
        cmd_name = 'run_cog_genes'
        self.run_code(cmd_name, cmd)

    def run_cog_summary(self):
        cmd = '{} {}'.format(self.python, self.cog_summary_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --cog {}'.format(self.cog_tsv)
        cmd += ' --txpt {}'.format(self.txpt_summary_tsv)
        cmd += ' --gene {}'.format(self.gene_summary_tsv)
        cmd_name = 'run_cog_summary'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.txpt_summary_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.gene_summary_tsv)
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
        self.option('cog_table').set_path(self.cog_tsv)
        self.option('cog_gene_table').set_path(self.cog_tsv)
        self.option('txpt_summary').set_path(self.txpt_summary_tsv)
        self.option('gene_summary').set_path(self.gene_summary_tsv)
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
            'id': 'eggnog_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.eggnog_annot',
            'instant': False,
            'options': {
                'blast_xml': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotFilter/output/eggnog/blast.filter.xml',
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
