# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import shutil
import unittest

from biocluster.module import Module


class AsprofileModule(Module):
    def __init__(self, work_id):
        super(AsprofileModule, self).__init__(work_id)
        options = [
            {'name': 'gtf_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'group_table', 'type': 'infile', 'format': 'ref_rna_v2.common'}

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
        super(AsprofileModule, self).run()
        self.run_count_ref()


    def run_count_ref(self):
        self.count_ref = self.add_tool('ref_rna_v3.ASprofile.count_ref')
        self.count_ref.set_options({'ref_fa': self.option('ref_fa')})
        self.count_ref.on('end', self.run_gtf_compare)
        self.count_ref.run()

    def run_gtf_compare(self):
        sp2fqs_dict = read_gtf_dir(self.option('gtf_dir').path)
        for sample, gtf_list in sp2fqs_dict.items():
            self.gffcompare = self.add_tool('ref_rna_v3.ASprofile.gffcompare')
            opts = {'merged_gtf': gtf_list,
                    'ref_gtf': self.option('ref_gtf'),
                    'sample': sample,
                    'seq_path': self.option('ref_fa')

                    }

            self.gffcompare.set_options(opts)
            self.tools.append(self.gffcompare)
        else:
            self.on_rely(self.tools, self.run_asprofile)
        for tool in self.tools:
            tool.run()

    def run_asprofile(self):
        gtf_list = os.path.join(self.work_dir, 'gtf.list')
        sp2fqs_dict = dict()
        for tool in self.tools:
            sp2fqs_dict[tool.option('sample')] = [tool.option('old_gtf').path, tool.option('new_gtf').path]

        for sample, gtf_list in sp2fqs_dict.items():
            self.asprofile = self.add_tool('ref_rna_v3.ASprofile.asprofile')
            opts = {'transcripts_old': gtf_list[0],
                    'transcripts_new': gtf_list[1],
                    'hdrs': self.count_ref.option('hdrs'),
                    'sample': sample,

                    }

            self.asprofile.set_options(opts)
            self.modules.append(self.asprofile)
        else:
            self.on_rely(self.modules, self.run_merge)
        for module in self.modules:
            module.run()

    def run_merge(self):

        self.merge = self.add_tool('ref_rna_v3.ASprofile.merge')
        result_list = os.path.join(self.work_dir, 'result.list')
        statistics_list = os.path.join(self.work_dir, 'statistics_list')
        group_dict = dict()
        for line in open(self.option('group_table').path):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                group_dict[sample] = group
        with open(result_list, 'w') as handle:
            for module in self.modules:
                handle.write('{}\t{}\t{}\n'.format(module.option('sample'), group_dict[module.option('sample')], module.option('as_result').path))
        with open(statistics_list, 'w') as handle1:
            for module in self.modules:
                handle1.write('{}\t{}\n'.format(module.option('sample'), module.option('as_statistics').path))
        opts = {
            'result_list': result_list,
            'statistics_list': statistics_list

        }

        self.merge.set_options(opts)
        self.merge.on('end', self.set_output)
        self.merge.run()

    def set_output(self):
        for file_name in os.listdir(self.merge.output_dir):
            source = os.path.join(self.merge.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.end()

    def end(self):
        super(AsprofileModule, self).end()


def read_gtf_dir(gtf_dir):

    sp2fqs_dict = dict()
    with open(os.path.join(gtf_dir, 'list'), 'r') as path:
        for line in path.readlines():
            gtf, sample = line.strip().split('\t')
            gtf_path = os.path.join(gtf_dir, gtf)
            sp2fqs_dict[sample]=gtf_path
    return sp2fqs_dict



class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_ath(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'asprofile_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'tool_lab.asprofile',
            'instant': False,
            'options': {
                'gtf_dir': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/gtf/',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
