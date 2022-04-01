# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun-20190927'

import unittest
import os
from biocluster.module import Module


class FindcircModule(Module):
    def __init__(self, work_id):
        super(FindcircModule, self).__init__(work_id)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'fq1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'fq2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}


        ]
        self.add_option(options)


    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(FindcircModule, self).run()
        self.run_bwa()

    def run_bwa(self):
        self.bwa = self.add_tool('whole_transcriptome.circrna.bwamem')
        if self.option('fq2').is_set:
            opts = {
                'sample': self.option('sample'),
                'genome': self.option('genome'),
                'fq1': self.option('fq1'),
                'fq2': self.option('fq2'),
                # 'rnasam': self.option('rnasam')
            }
        else:
            opts = {
                'sample': self.option('sample'),
                'genome': self.option('genome'),
                'fq1': self.option('fq1'),
                # 'rnasam': self.option('rnasam')
            }
        self.bwa.set_options(opts)
        self.bwa.on('end', self.run_findcircbwa)
        self.bwa.run()

    def run_findcircbwa(self):
        self.findcircbwa = self.add_tool('whole_transcriptome.circrna.findcircbwa')
        opts = {
            'sample' : self.option('sample'),
            'genome' : self.option('genome'),
            'rnasam' : self.bwa.option('rnasam'),
            'junction_reads': self.option('junction_reads'),
            'circrna_length': self.option('circrna_length')
        }
        self.findcircbwa.set_options(opts)
        self.findcircbwa.on('end',self.run_signal)
        self.findcircbwa.run()


    def run_signal(self):
        self.signal = self.add_tool('whole_transcriptome.circrna.signal')
        opts = {
            'sample': self.option('sample'),
            'genome': self.option('genome'),
            'circRNA_merge': self.findcircbwa.option('find_circ_filter'),
            # 'bed': self.option('bed'),
            # 'fasta': self.option('fasta'),
            # 'signal': self.option('signal')
        }
        self.signal.set_options(opts)
        self.signal.on('end', self.run_circtype)
        self.signal.run()

    def run_circtype(self):
        self.circtype = self.add_tool('whole_transcriptome.circrna.circtype')
        opts = {
            'sample': self.option('sample'),
            'annotate': self.option('annotate'),
            # 'type': self.option('type'),
            'signal': self.signal.option('signal')
        }
        self.circtype.set_options(opts)
        self.circtype.on('end', self.run_rpm)
        self.circtype.run()


    def run_rpm(self):
        self.rpm = self.add_tool('whole_transcriptome.circrna.rpm')
        opts = {
            'sample': self.option('sample'),
            'BSJ': self.circtype.option('type'),
            'rnasam': self.bwa.option('rnasam'),

        }
        self.rpm.set_options(opts)
        self.rpm.on('end',self.set_output)
        self.rpm.run()

    def set_output(self):
        p = self.rpm.option('RPM').path
        link_names = os.path.join(self.output_dir, os.path.basename(p))
        if os.path.exists(link_names):
            os.remove(link_names)
        os.link(p, link_names)
        self.rpm.option('RPM').set_path(link_names)
        self.end()

    def end(self):
        super(FindcircModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'findcirc_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.circrna.findcirc',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'fq1': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_l.fastq',
                'fq2': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_r.fastq',
                'annotate': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf',
                'sample': 'zjx',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
