# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class GoAnnotAgent(Agent):
    '''
    last_modify: 2019.04.16
    '''
    def __init__(self, parent):
        super(GoAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'go'},
            {'name': 'txpt_id2terms', 'type': 'outfile', 'format': 'ref_rna_v2.go_list'},
            {'name': 'gene_id2terms', 'type': 'outfile', 'format': 'ref_rna_v2.go_list'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "pir_version", "type": "string", "default": "2019"}, #pir database version
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version

        ]
        self.add_option(options)
        self.step.add_steps('go_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST2GO'

    def step_start(self):
        self.step.go_annot.start()
        self.step.update()

    def step_end(self):
        self.step.go_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        infile_size = os.path.getsize(self.option('blast2go_annot').prop['path'])
        self._memory = '{}G'.format(int(float(infile_size) / 1024 ** 3 * 8 + 10))

    def end(self):
        super(GoAnnotAgent, self).end()

class GoAnnotTool(Tool):
    def __init__(self, config):
        super(GoAnnotTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.go_transition_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/go_transition.py')
        self.go_relation_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/go_relation.py')
        self.go_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/go_annotation2.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_venn_list.py')
        self.prefix = os.path.basename(self.option('blast2go_annot').path)[:-4]
        self.b2g_txpt_tsv = os.path.join(self.work_dir, '{}.T.tsv'.format(self.prefix))
        self.b2g_gene_tsv = os.path.join(self.work_dir, '{}.G.tsv'.format(self.prefix))
        self.id2terms_txpt_tsv = os.path.join(self.output_dir, 'id2terms.T.tsv')
        self.id2terms_gene_tsv = os.path.join(self.output_dir, 'id2terms.G.tsv')
        self.b2g_host = 'localhost'
        self.b2g_user = 'biocluster102'
        self.b2g_passwd = 'sanger-dev-123'
        self.b2g_db = 'b2gdb'
        self.txpt_dir = os.path.join(self.output_dir, 'T')
        self.gene_dir = os.path.join(self.output_dir, 'G')
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']

    def run(self):
        super(GoAnnotTool, self).run()
        self.run_go_transition()
        self.run_go_relation_t()
        self.run_go_relation_g()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.run_go_annotation_t()
        self.run_go_annotation_g()
        self.set_output()
        self.end()

    def run_go_transition(self):
        shutil.copy(self.option('blast2go_annot').path, self.b2g_txpt_tsv)
        cmd = '{} {}'.format(self.python, self.go_transition_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --blast2go {}'.format(self.b2g_txpt_tsv)
        cmd += ' --output {}'.format(self.b2g_gene_tsv)
        cmd_name = 'run_go_genes'
        self.run_code(cmd_name, cmd)

    def run_go_relation_t(self):
        cmd = '{} {}'.format(self.python, self.go_relation_py)
        cmd += ' -i {}'.format(self.b2g_txpt_tsv)
        cmd += ' -o {}'.format(self.id2terms_txpt_tsv)
        cmd_name = 'run_go_merge_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_go_relation_g(self):
        cmd = '{} {}'.format(self.python, self.go_relation_py)
        cmd += ' -i {}'.format(self.b2g_gene_tsv)
        cmd += ' -o {}'.format(self.id2terms_gene_tsv)
        cmd_name = 'run_go_merge_g'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.id2terms_txpt_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.id2terms_gene_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.gene_list)
        cmd_name = 'run_get_venn_list_g'
        self.run_code(cmd_name, cmd, block=False)

    def run_go_annotation_t(self):
        if os.path.isdir(self.txpt_dir):
            shutil.rmtree(self.txpt_dir)
        os.mkdir(self.txpt_dir)
        '''
        cmd = '{} {} {} {} {} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_txpt_tsv,
            self.b2g_host,
            self.b2g_user,
            self.b2g_passwd,
            self.b2g_db,
            self.txpt_dir
        )
        '''
        cmd = '{} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_txpt_tsv,
            self.txpt_dir,
            self.go_obo
        )
        cmd_name = 'run_go_annot_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_go_annotation_g(self):
        if os.path.isdir(self.gene_dir):
            shutil.rmtree(self.gene_dir)
        os.mkdir(self.gene_dir)
        '''
        cmd = '{} {} {} {} {} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_gene_tsv,
            self.b2g_host,
            self.b2g_user,
            self.b2g_passwd,
            self.b2g_db,
            self.gene_dir
        )
        '''
        cmd = '{} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_gene_tsv,
            self.gene_dir,
            self.go_obo
        )
        cmd_name = 'run_go_annot_g'
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
                        self.set_error('fail to run %s, abord', variables=(n), code="33710402")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('txpt_id2terms').set_path(self.id2terms_txpt_tsv)
        self.option('gene_id2terms').set_path(self.id2terms_gene_tsv)
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
            'id': 'go_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.go_annot',
            'instant': False,
            'options': {
                'blast2go_annot': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotFilter/output/go/blast2go.filter.xls',
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
