# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.module import Module


class CircrnaModule(Module):
    def __init__(self, work_id):
        super(CircrnaModule, self).__init__(work_id)
        CIRC_METHOD = ('ciri2,find_circ', 'circ2', 'find_circ', 'circ_finder', 'circexplorer')
        options = [
            {'name': 'circ_method', 'type': 'string', 'default': CIRC_METHOD[0]},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
        ]
        self.add_option(options)
        self.modules = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(CircrnaModule, self).run()
        if self.option('circ_method') == 'ciri2,find_circ':
            self.run_circaddfindcirc()

    def run_circaddfindcirc(self):
        is_se, sp2fqs_dict = read_fastq_dir(self.option('fastq_dir').path)
        for sample, fastq_list in sp2fqs_dict.items():
            ciriaddfindcirc = self.add_module('whole_transcriptome.circrna.ciriaddfindcirc')
            opts = {'genome': self.option('genome').path, 'annotate': self.option('annotate').path, 'sample': sample}
            if is_se:
                opts.update({'fq1': fastq_list[0]})
            else:
                opts.update({'fq1': fastq_list[0], 'fq2': fastq_list[1]})
            ciriaddfindcirc.set_options(opts)
            self.modules.append(ciriaddfindcirc)
        else:
            self.on_rely(self.modules, self.set_output)
        for module in self.modules:
            module.run()

    def set_output(self):
        for file_name in os.listdir(self.trans_build.output_dir):
            source = os.path.join(self.trans_build.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        else:
            self.end()

    def end(self):
        super(CircrnaModule, self).end()

def read_fastq_dir(fastq_dir):
    is_se = False
    sp2fqs_dict = dict()
    for line in open(os.path.join(fastq_dir, 'list.txt')):
        eles = line.strip().split('\t')
        fastq = os.path.join(fastq_dir, eles[0])
        sample, mate_type = eles[1:]
        if sample in sp2fqs_dict:
            if mate_type == 'l':
                sp2fqs_dict[sample].insert(0, fastq)
            elif mate_type == 'r':
                sp2fqs_dict[sample].append(fastq)
        else:
            sp2fqs_dict[sample] = [fastq]
    else:
        pe_sample_count = len(filter(lambda item: len(item[1]) > 1, sp2fqs_dict.items()))
        if pe_sample_count:
            if pe_sample_count == len(sp2fqs_dict):
                is_se = False
            else:
                raise Exception('mix mate type found in {} -> {}'.format(fastq_dir, sp2fqs_dict))
        else:
            is_se = True
        return is_se, sp2fqs_dict


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.circrna.circrna',
            'instant': False,
            'options': {
                'circ_method': 'ciri2,find_circ',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191010/Longrna_workflow_2617_9288/FastpRna/output/fastq',
                'annotate': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
