# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import shutil
import unittest

class BuildIdxAgent(Agent):
    '''
    last_modify: 2019.03.12
    '''
    def __init__(self, parent):
        super(BuildIdxAgent, self).__init__(parent)
        options = [
            {'name': 'ref_fa', 'type': 'infile', 'format':'lnc_rna.fasta'},
            {'name': 'reference_fasta', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('build_idx')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.build_idx.start()
        self.step.update()

    def step_end(self):
        self.step.build_idx.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('ref_fa').is_set:
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').prop['path']))
            self.infile_size = os.path.getsize(self.option('ref_fa').prop['path'])
        else:
            raise OptionError('reference FASTA file must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 4 + 12))

    def end(self):
        super(BuildIdxAgent, self).end()

class BuildIdxTool(Tool):
    def __init__(self, config):
        super(BuildIdxTool, self).__init__(config)
        self.samtools = 'bioinfo/align/samtools-1.6/samtools'
        self.java = 'program/sun_jdk1.8.0/bin/java'
        self.picard_jar = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/gene-structure/picard.jar')
        self.reference_fasta = os.path.join(self.output_dir, os.path.basename(
                '{}.fasta'.format(self.option('ref_fa').path[:self.option('ref_fa').path.rfind('.')]
        )))
        shutil.copy(self.option('ref_fa').path, self.reference_fasta)
        self.reference_dict = '{}.dict'.format(self.reference_fasta[:-6])
        self.reference_fasta_fai = '{}.fai'.format(self.reference_fasta)

    def run(self):
        super(BuildIdxTool, self).run()
        self.run_samtools_faidx()
        self.run_picard_create_sequence_dictionary()
        self.set_output()
        self.end()

    def run_samtools_faidx(self):
        cmd = '{} faidx {}'.format(self.samtools, self.reference_fasta)
        cmd_name = 'run_samtools_faidx'
        self.run_code(cmd_name, cmd)

    def run_picard_create_sequence_dictionary(self):
        cmd = '{} -jar {} CreateSequenceDictionary'.format(self.java, self.picard_jar)
        cmd += ' R={}'.format(self.reference_fasta)
        cmd += ' O={}'.format(self.reference_dict)
        cmd_name = 'run_picard_create_sequence_dictionary'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
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
        self.option('reference_fasta').set_path(self.reference_fasta)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'build_idx_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.build_idx',
            'instant': False,
            'options': {
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()