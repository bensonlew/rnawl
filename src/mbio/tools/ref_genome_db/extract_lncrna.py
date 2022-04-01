# -*- coding: utf-8 -*-
# __author__ = 'scp'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import subprocess
import unittest

class ExtractLncrnaAgent(Agent):
    def __init__(self, parent):
        super(ExtractLncrnaAgent, self).__init__(parent)
        options = [
            {'name': 'reference_in', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            {"name": "gtf_in", "type": "infile", "format": "ref_genome_db.gtf"},
            {"name": "gene_biotype", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "trans_biotype", "type": "infile", "format": "ref_genome_db.common"},
            {'name': 'fasta_out', 'type': 'outfile', 'format': 'ref_genome_db.fasta'},
            {'name': 'gtf_out', 'type': 'outfile', 'format': 'ref_genome_db.gtf'},
            {"name": "lncrna_list", "type": "outfile", "format": "ref_genome_db.common"},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('reference_in', self.option('reference_in').prop['path']))
        self.logger.debug('{} - {}'.format('gtf_in', self.option('gtf_in').prop['path']))
        self.logger.debug('{} - {}'.format('gene_biotype', self.option('gene_biotype').prop['path']))
        self.logger.debug('{} - {}'.format('trans_biotype', self.option('trans_biotype').prop['path']))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(ExtractLncrnaAgent, self).end()

class ExtractLncrnaTool(Tool):
    def __init__(self, config):
        super(ExtractLncrnaTool, self).__init__(config)
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.perl = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.extract_lncrna_gtf = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db/extract_lncrna_gtf.pl')

    def run(self):
        super(ExtractLncrnaTool, self).run()
        self.run_extract_lncrna_gtf()
        self.run_extract_lncrna_fa()
        self.set_output()
        self.end()

    def run_extract_lncrna_gtf(self):
        self.gtf_out = os.path.join(self.work_dir, "lncrna.gtf")
        self.lncrna_list = os.path.join(self.work_dir, "lncrna.txt")
        cmd = '{} {} {} {} {} {} {}'.format(
            self.perl,
            self.extract_lncrna_gtf,
            self.option('gtf_in').prop['path'],
            self.option('gene_biotype').prop['path'],
            self.option('trans_biotype').prop['path'],
            self.gtf_out,
            self.lncrna_list
        )
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}运行出错".format(cmd)
        cmd_name = 'extract_lncrna_gtf'
        self.run_code(cmd_name, cmd)

    def run_extract_lncrna_fa(self):
        self.fa_out = os.path.join(self.work_dir, "lncrna.fa")
        cmd = self.cufflinks_path + "gffread %s -g %s -w %s" % (
            self.gtf_out,
            self.option('reference_in').prop['path'],
            self.fa_out)
        cmd_name = 'extract_lncrna_fa'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
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
        self.option("fasta_out", self.fa_out)
        self.option("gtf_out", self.gtf_out)
        self.option("lncrna_list", self.lncrna_list)
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
            'id': 'extract_lncrna_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.extract_lncrna',
            'instant': False,
            'options': {
                'reference_in': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Homo_sapiens/ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'gtf_in': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Homo_sapiens/ensembl/Homo_sapiens.GRCh38.96.gtf',
                'gene_biotype': '/mnt/ilustre/users/sanger-dev/workspace/20190530/Single_Extract_6634/ExtractBiotype/gene_biotype.txt',
                'trans_biotype': '/mnt/ilustre/users/sanger-dev/workspace/20190530/Single_Extract_6634/ExtractBiotype/trans_biotype.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()