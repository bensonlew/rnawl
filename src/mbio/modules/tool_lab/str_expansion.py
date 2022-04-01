# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import shutil
import unittest

from biocluster.module import Module


class StrExpansionModule(Module):
    def __init__(self, work_id):
        super(StrExpansionModule, self).__init__(work_id)
        options = [
            {'name': 'bamlist', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'variant-catalog', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'sample_list', 'type': 'infile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.modules = list()
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(StrExpansionModule, self).run()
        self.run_str_expansion()

    def run_str_expansion(self):
        with open(self.option('bamlist').path, 'r') as b:
            for bam in b.readlines():
                bam_name = os.path.basename(bam.strip()).split('.bam')[0]
                self.str = self.add_tool('tool_lab.str_expansionhunter')
                opts = {'bam': bam.strip(),
                        'ref': self.option('ref'),
                        'bam_name': bam_name,
                        'variant-catalog': self.option('variant-catalog')
                        }
                self.str.set_options(opts)
                self.tools.append(self.str)
            else:
                self.on_rely(self.tools, self.end)
            for tool in self.tools:
                tool.run()

    def set_output(self):
        for file_name in os.listdir(self.str_merge.output_dir):
            source = os.path.join(self.str_merge.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.end()

    def end(self):
        super(StrExpansionModule, self).end()




class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_ath(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Str_expansion_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'tool_lab.str_expansion',
            'instant': False,
            'options': {
                'bamlist': '/mnt/ilustre/users/sanger-dev/workspace/20210202/Str_STR7071/RnaseqMapping/output/bamlist',
                'ref': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/SL4.0_ITAG4.0/dna/S_lycopersicum_chromosomes.4.00.fa',
                'variant-catalog': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/variant_catalog/variant_catalog_ssr_exon.json'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
