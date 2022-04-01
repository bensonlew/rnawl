# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import shutil
import unittest

class PrecursorCandidateAgent(Agent):
    '''
    last_modify: 2019.03.27
    '''
    def __init__(self, parent):
        super(PrecursorCandidateAgent, self).__init__(parent)
        options = [
            {'name': 'file_in', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'outfmt', 'type': 'int', 'default': 6},
            {'name': 'evalue', 'type': 'int', 'default': 1},
            {'name': 'num_threads', 'type': 'int', 'default': 2},
            {'name': 'tblout', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('precursor_candidate')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.precursor_candidate.start()
        self.step.update()

    def stepfinish(self):
        self.step.precursor_candidate.finish()
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
        super(PrecursorCandidateAgent, self).end()

class PrecursorCandidateTool(Tool):
    def __init__(self, config):
        super(PrecursorCandidateTool, self).__init__(config)
        self.makeblastdb = 'bioinfo/align/ncbi-blast-2.3.0+/bin/makeblastdb'
        self.blastn = 'bioinfo/align/ncbi-blast-2.3.0+/bin/blastn'
        self.python = 'program/Python/bin/python'
        self.file_in = self.option('file_in').prop['path']
        self.hairpin_fa = os.path.join(
            self.config.SOFTWARE_DIR,
            'database/mirbase/species/{}/{}.hairpin.fa'.format(self.option('species'), self.option('species'))
        )
        self.precursor_candidate_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/precursor_candidate.py')
        self.database = os.path.join(self.work_dir, os.path.basename(self.hairpin_fa))
        self.blastn_tblout = os.path.join(self.work_dir, 'blastn.tblout')

    def run(self):
        super(PrecursorCandidateTool, self).run()
        self.run_makeblastdb()
        self.run_blastn()
        self.set_output()
        self.end()

    def run_makeblastdb(self):
        shutil.copy(self.hairpin_fa, self.database)
        cmd = '{} -in {} -dbtype nucl'.format(self.makeblastdb, self.database)
        cmd_name = 'run_makeblastdb'
        self.run_code(cmd_name, cmd)

    def run_blastn(self):
        cmd = '{} -outfmt {}'.format(self.blastn, self.option('outfmt'))
        cmd += ' -query {}'.format(self.file_in)
        cmd += ' -db {}'.format(self.database)
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -out {}'.format(self.blastn_tblout)
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastn'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.blastn_tblout
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('tblout').set_path(link_name)
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
            'id': 'precursor_candidate_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.precursor_candidate',
            'instant': False,
            'options': {
                'species': 'hsa',
                'file_in': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/mirna_precursor/all_lncrna.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()