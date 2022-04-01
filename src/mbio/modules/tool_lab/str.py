# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import shutil
import unittest

from biocluster.module import Module


class StrModule(Module):
    def __init__(self, work_id):
        super(StrModule, self).__init__(work_id)
        options = [
            {'name': 'bamlist', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
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
        super(StrModule, self).run()
        self.run_str_predict()

    def run_str_predict(self):
        with open(self.option('bamlist').path, 'r') as b:
            for bam in b.readlines():
                bam_name = os.path.basename(bam.strip()).split('.bam')[0]
                self.str = self.add_tool('tool_lab.str')
                opts = {'bam': bam.strip(),
                        'ref': self.option('ref_fa'),
                        'bam_name': bam_name
                        }
                self.str.set_options(opts)
                self.tools.append(self.str)
            else:
                self.on_rely(self.tools, self.run_str_merge)
            for tool in self.tools:
                tool.run()

    def run_str_merge(self):
        self.manifest = os.path.join(self.work_dir, 'manifest.tsv')
        sample_dict = dict()
        with open(self.option('sample_list').path, 'r') as s:
            for lines in s.readlines():
                sample, sample_type = lines.strip().split('\t')
                sample_dict[sample] = sample_type
        with open(self.manifest, 'w') as m:
            for tool in self.tools:
                m.write(tool.option('bam_name') + '\t' + sample_dict[tool.option('bam_name')] + '\t' + tool.option('str_json').path + '\n')
        self.str_merge = self.add_tool('tool_lab.str_merge')
        opts = {'manifest': self.manifest,
                'ref': self.option('ref_fa'),
                }
        self.str_merge.set_options(opts)
        self.str_merge.on('end', self.run_str_case_control)
        self.str_merge.run()

    def run_str_case_control(self):
        self.str_case_control = self.add_tool('tool_lab.str_case_control')
        opts = {
            'manifest': self.manifest,
            'merge_json': self.str_merge.option('merge_json')
        }
        self.str_case_control.set_options(opts)
        self.str_case_control.on('end', self.set_output)
        self.str_case_control.run()




    def set_output(self):
        for file_name in os.listdir(self.str_merge.output_dir):
            source = os.path.join(self.str_merge.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.end()

    def end(self):
        super(StrModule, self).end()




class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_ath(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Str_predict_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'tool_lab.str',
            'instant': False,
            'options': {
                'bamlist': '/mnt/ilustre/users/sanger-dev/workspace/20210114/Str_STR2926/RnaseqMapping/output/bamlist',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/SL4.0_ITAG4.0/dna/S_lycopersicum_chromosomes.4.00.fa',
                'sample_list': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/sample_list.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
