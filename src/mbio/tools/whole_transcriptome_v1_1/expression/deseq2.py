# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class Deseq2Agent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(Deseq2Agent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'together', 'type': 'bool', 'default': True},
            {'name': 'is_batch', 'type': 'bool', 'default': False},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(Deseq2Agent, self).end()


class Deseq2Tool(Tool):
    def __init__(self, config):
        super(Deseq2Tool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'deseq2': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome_v1_1/expression/deseq2.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'design': os.path.join(self.work_dir, 'design.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'contrast': os.path.join(self.work_dir, 'contrast.txt')
        }

    def run(self):
        super(Deseq2Tool, self).run()
        self.pre_deseq2()
        self.run_deseq2()
        self.set_output()
        self.end()

    def pre_deseq2(self):
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('control_table').path, self.file['contrast'])
        samples = list()
        groups = list()
        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)
        if self.option('is_batch') == True:
            batchs_dict = dict()
            for line in open(self.option('batch_matrix').path):
                if line.strip() and line[0] != '#':
                    sample, batch = line.strip().split('\t')
                    batchs_dict[sample] = batch
            batchs = [batchs_dict[i] for i in samples]
        contrasts = list()
        for line in open(self.file['contrast']):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')
                contrasts.append({'ctrl': ctrl, 'case': case})
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df = df.astype(int)
        df.to_csv(self.file['count'], sep='\t')
        opts = {
            'count': self.file['count'],
            'samples': samples,
            'groups': groups,
            'contrasts': contrasts,
            'together': self.option('together'),
            'output': self.output_dir,
            'is_batch': self.option('is_batch')
        }
        if self.option('is_batch') == True:
            opts.update({
                'batchs': batchs
            })
        json.dump(opts, open(self.file['json'], 'w'), indent=4)

    def run_deseq2(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['deseq2'], self.file['json'])
        runcmd(self, 'run_deseq2', cmd)

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'deseq2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.deseq2',
            'instant': False,
            'options': {
                'count_matrix': PARAMS['count_matrix'],
                'group_table': PARAMS['group_table'],
                'control_table': PARAMS['control_table'],
                'together': PARAMS['together'],
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run DESeq2 on biocluster')
    parser.add_argument('json', action='store', help='setting file in JSON format')

    args = parser.parse_args()

    if hasattr(args, 'json'):
        global PARAMS
        PARAMS = json.load(open(args.json))
        dirname = os.path.dirname(os.path.abspath(args.json))
        for k in ('count_matrix', 'group_table', 'control_table'):
            PARAMS[k] = os.path.join(dirname, PARAMS[k])
        suite = unittest.TestSuite()
        suite.addTests([TestFunction('test')])
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        parser.print_help()
        sys.exit(-1)
