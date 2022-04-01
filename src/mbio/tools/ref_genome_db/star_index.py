# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import unittest

class StarIndexAgent(Agent):
    def __init__(self, parent):
        super(StarIndexAgent, self).__init__(parent)
        options = [
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'genome_dir', 'type': 'outfile', 'format': 'ref_genome_db.common_dir'},
            {'name': 'genome_fasta_files', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            {'name': 'sjdb_gtf_file', 'type': 'infile', 'format': 'ref_genome_db.gtf'},
            {'name': 'sjdb_overhang', 'type': 'int', 'default': 149},
        ]
        self.add_option(options)
        self._memory_increase_step = 60

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('threads', self.option('threads')))
        self.logger.debug('{} - {}'.format('genome_fasta_files', self.option('genome_fasta_files').prop['path']))
        self.logger.debug('{} - {}'.format('sjdb_gtf_file', self.option('sjdb_gtf_file').prop['path']))
        self.logger.debug('{} - {}'.format('sjdb_overhang', self.option('sjdb_overhang')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        self._memory = '{}G'.format(os.path.getsize(self.option('genome_fasta_files').prop['path'])/1024**3+120)

    def end(self):
        super(StarIndexAgent, self).end()

class StarIndexTool(Tool):
    def __init__(self, config):
        super(StarIndexTool, self).__init__(config)
        self.star = 'bioinfo/rna/star-2.5/bin/Linux_x86_64/STAR'
        self.star27 = 'bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/STAR'
        self.star_version = "2.5"

    def run(self):
        super(StarIndexTool, self).run()
        self.run_star_index()
        self.set_output()
        self.end()

    def run_star_index(self):
        cmd = '{} '.format(self.star)
        cmd += '--runThreadN {} '.format(self.option('threads'))
        cmd += '--runMode genomeGenerate '
        cmd += '--genomeDir {} '.format(self.output_dir)
        cmd += '--genomeFastaFiles {} '.format(self.option('genome_fasta_files').prop['path'])
        cmd += '--sjdbGTFfile {} '.format(self.option('sjdb_gtf_file').prop['path'])
        cmd += '--sjdbOverhang {} '.format(self.option('sjdb_overhang'))
        cmd += '--limitGenomeGenerateRAM {} '.format(1200000000000)
        cmd_name = 'star_index'
        self.run_code(cmd_name, cmd)

    def run_star_version27(self):
        self.star_version = "2.7"
        cmd = '{} '.format(self.star27)
        cmd += '--runThreadN {} '.format(self.option('threads'))
        cmd += '--runMode genomeGenerate '
        cmd += '--genomeDir {} '.format(self.output_dir)
        cmd += '--genomeFastaFiles {} '.format(self.option('genome_fasta_files').prop['path'])
        cmd += '--sjdbGTFfile {} '.format(self.option('sjdb_gtf_file').prop['path'])
        cmd += '--sjdbOverhang {} '.format(self.option('sjdb_overhang'))
        cmd += '--limitGenomeGenerateRAM {} '.format(1200000000000)
        cmd_name = 'star_index27'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd, ignore_error=True)
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
        elif command.return_code in [105] and self.star_version == "2.5":
            self.run_star_version27()
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('genome_dir', self.output_dir)
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
            'id': 'star_index_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.star_index',
            'instant': False,
            'options': {
                'genome_fasta_files': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/Homo_sapiens.GRCh38.94.fa',
                'sjdb_gtf_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/Homo_sapiens.GRCh38.94.gtf',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
