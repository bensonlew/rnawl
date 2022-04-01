# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import shutil
import unittest

class LncrnaOrthologAgent(Agent):
    '''
    last_modify: 2019.04.23
    '''
    def __init__(self, parent):
        super(LncrnaOrthologAgent, self).__init__(parent)
        options = [
            {'name': 'lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'target_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'dbtype', 'type': 'string', 'default': 'nucl'},
            {'name': 'outfmt', 'type': 'int', 'default': 6},
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'identity', 'type': 'float', 'default': 50.0},
            {'name': 'num_threads', 'type': 'int', 'default': 2},
            {'name': 'type_tsv', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('lncrna_ortholog')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.lncrna_ortholog.start()
        self.step.update()

    def stepfinish(self):
        self.step.lncrna_ortholog.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '4G'

    def end(self):
        super(LncrnaOrthologAgent, self).end()

class LncrnaOrthologTool(Tool):
    def __init__(self, config):
        super(LncrnaOrthologTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.makeblastdb = 'bioinfo/align/ncbi-blast-2.3.0+/bin/makeblastdb'
        self.blastn = 'bioinfo/align/ncbi-blast-2.3.0+/bin/blastn'
        self.filter_fasta_by_id_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_fasta_by_id.py')
        self.ortholog_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/ortholog_statistics.py')
        self.test_fa = os.path.join(self.work_dir, 'test.fa')
        self.ctrl_fa = os.path.join(self.work_dir, 'ctrl.fa')
        self.forward_tblout = os.path.join(self.work_dir, 'forward.tblout')
        self.reverse_tblout = os.path.join(self.work_dir, 'reverse.tblout')
        self.tabular = os.path.join(self.output_dir, 'lncRNA_ortholog.tabular')

    def run(self):
        super(LncrnaOrthologTool, self).run()
        self.run_filter_fasta_by_id()
        self.run_makeblastdb_test()
        self.run_makeblastdb_ctrl()
        self.run_blastn_forward()
        self.run_blastn_reverse()
        self.run_ortholog_statistics()
        self.set_output()
        self.end()

    def run_filter_fasta_by_id(self):
        cmd = '{} {}'.format(self.python, self.filter_fasta_by_id_py)
        cmd += ' -i {}'.format(self.option('lncrna_fa').path)
        cmd += ' -l {}'.format(self.option('keep_list').path)
        cmd += ' -o {}'.format(self.test_fa)
        cmd_name = 'run_filter_fasta_by_id'
        self.run_code(cmd_name, cmd)

    def run_makeblastdb_test(self):
        cmd = '{} -in {} -dbtype {}'.format(self.makeblastdb, self.test_fa, self.option('dbtype'))
        cmd_name = 'run_makeblastdb_test'
        self.run_code(cmd_name, cmd)

    def run_makeblastdb_ctrl(self):
        shutil.copy(self.option('target_fa').path, self.ctrl_fa)
        cmd = '{} -in {} -dbtype {}'.format(self.makeblastdb, self.ctrl_fa, self.option('dbtype'))
        cmd_name = 'run_makeblastdb_ctrl'
        self.run_code(cmd_name, cmd)

    def run_blastn_forward(self):
        cmd = '{} -outfmt {}'.format(self.blastn, self.option('outfmt'))
        cmd += ' -query {}'.format(self.test_fa)
        cmd += ' -db {}'.format(self.ctrl_fa)
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -out {}'.format(self.forward_tblout)
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastn_forward'
        self.run_code(cmd_name, cmd)

    def run_blastn_reverse(self):
        cmd = '{} -outfmt {}'.format(self.blastn, self.option('outfmt'))
        cmd += ' -query {}'.format(self.ctrl_fa)
        cmd += ' -db {}'.format(self.test_fa)
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -out {}'.format(self.reverse_tblout)
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastn_reverse'
        self.run_code(cmd_name, cmd)

    def run_ortholog_statistics(self):
        cmd = '{} {}'.format(self.python, self.ortholog_statistics_py)
        cmd += ' -f {}'.format(self.forward_tblout)
        cmd += ' -r {}'.format(self.reverse_tblout)
        cmd += ' -i {}'.format(self.option('identity'))
        cmd += ' -t {}'.format(self.option('type_tsv').path)
        cmd += ' -o {}'.format(self.tabular)
        cmd_name = 'run_ortholog_statistics'
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
        self.option('tabular').set_path(self.tabular)
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
            'id': 'lncrna_ortholog_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.lncrna_ortholog',
            'instant': False,
            'options': {
                'lncrna_fa': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33833/upload_results/Filtered_result/filtered_file/all_lncrna.fa',
                'keep_list': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33833/upload_results/Filtered_result/filtered_file/known_lncrna_ids.list',
                'target_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/lncrna/lncrna.fa',
                'evalue': 1e-5,
                'identity': 50,
                'type_tsv': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33833/upload_results/Filtered_result/filtered_file/trans_type.xls'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
