# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,shicaiping,qinjincheng'

import glob
import os
import shutil
import unittest
import collections
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
import random
import matplotlib.pyplot as plt


class CoverageRandomModule(Module):
    def __init__(self, work_id):
        super(CoverageRandomModule, self).__init__(work_id)
        options = [
            {'name': 'bed', 'type': 'infile', 'format': 'ref_rna_v2.bed'},  # bed格式文件
            {'name': 'bam', 'type': 'infile', 'format': 'align.bwa.bam'},  # bam格式文件,排序过的
            {'name': 'min_len', 'type': 'int', 'default': 100}  # Minimum mRNA length (bp).
        ]
        self.add_option(options)
        self.tools = []
        self.files = []

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def split_bed(self):
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)
        self.bed_list = split_file + '/bedlist'

        bed_file = pd.read_table(self.option('bed').path, header=None, index_col=None, sep='\t')
        lines_number = bed_file.shape[0]
        need_number = lines_number//10
        with open(self.bed_list, 'w+') as bl:
            split_name = split_file + '/' + os.path.basename(self.option('bed').path) + 'filter'
            need_list = random.sample(range(0,lines_number), need_number)
            need_bed = bed_file.iloc[need_list]
            need_bed.to_csv(split_name, header=None, index=None, sep='\t')
            bl.write(split_name + '\n')

    def coverage_new_run(self):
        with open(self.bed_list, 'r') as bed_list:
            for b in bed_list.readlines():
                self.coverage = self.add_tool('ref_rna_v3.mapping.coverage')
                self.coverage.set_options({
                    'bam': self.option('bam').path,
                    'bed': b.strip(),
                    'min_len': self.option('min_len')
                })
                self.tools.append(self.coverage)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
            for t in self.tools:
                t.run()
        else:
            self.tools[0].on('end', self.set_output)
            self.tools[0].run()
    # def set_output(self):
    #     coverage_dict = collections.defaultdict(int)
    #     reads_all = list()
    #     for tool in self.tools:
    #         geneBodyCoverage = glob.glob(tool.output_dir + '/*geneBodyCoverage.txt')[0]
    #         self.file_name = os.path.basename(geneBodyCoverage)
    #         self.sample_name = self.file_name.strip().split('geneBodyCoverage.txt')[0]
    #         with open(geneBodyCoverage, 'r') as c:
    #             line = c.readlines()[1]
    #             name = line.strip().split('\t')[0]
    #             reads = line.strip().split('\t')[1:]
    #             reads_all.append(reads)
    #     for i in reads_all:
    #         for j in range(0, 100):
    #             coverage_dict[j] += int(float(i[j]))
    #     with open(os.path.join(self.output_dir, self.file_name), 'w+') as o:
    #         o.write("Percentile\t" + '\t'.join([str(i) for i in range(1,101)]))
    #         o.write('\n' + name + '\t' + '\t'.join([str(coverage_dict[i]) for i in range(0,100)]))
    #     self.end()

    def set_output(self):
        for tool in self.tools:
            geneBodyCoverage = glob.glob(tool.output_dir + '/*geneBodyCoverage.txt')[0]
            coverage = pd.read_table(geneBodyCoverage, header=None, index_col=0, sep='\t')
            y = coverage.iloc[1].tolist()
            x = coverage.iloc[0].tolist()
            plt.plot(x, y)
        plt.savefig(os.path.join(self.work_dir, 'test.jpg'))
        self.end()


    def run(self):
        super(CoverageRandomModule, self).run()
        self.split_bed()
        self.samtool_index()


    def end(self):
        super(CoverageRandomModule, self).end()


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'coverage_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'ref_rna_v3.coverage_random',
            'instant': False,
            'options': dict(
                bam='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/neg_2.bam',
                bed='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/Homo_sapiens.GRCh38.101.gtf.filter.bed',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
