# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.module import Module
from mbio.packages.ref_rna_v3.functions import modlfuncdeco
import os
import unittest


class ExtractSeqModule(Module):
    '''
    last_modify: 20200708
    '''
    def __init__(self, work_id):
        super(ExtractSeqModule, self).__init__(work_id)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "extract_type", "type": "string", "default": 'mapped'},
            {"name": "seq_type", "type": "string", "default": 'PE'},
        ]
        self.add_option(options)

    @modlfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @modlfuncdeco
    def run(self):
        super(ExtractSeqModule, self).run()
        self.run_tool()

    @modlfuncdeco
    def run_tool(self):
        self.tools = [self.set_tool(n, line) for n, line in enumerate(open(self.option('input_file').path))]
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    @modlfuncdeco
    def set_tool(self, n, line):
        items = line.strip().split('\t')
        options = {
            'input_file': items[0],
            'seq_type': self.option('seq_type'),
            'extract_type': self.option('extract_type'),
        }
        if len(items) == 2:
            options.update({'basename': items[1]})
        else:
            options.update({'basename': items[0].split("/")[-1].replace('.bam', '')})
        self.step.add_steps('extract_{}'.format(n))
        extract = self.add_tool('ref_rna_v3.extract_seq')
        extract.set_options(options)
        extract.on('start', self.set_step, {'start': getattr(self.step, 'extract_{}'.format(n))})
        extract.on('end', self.set_step, {'end': getattr(self.step, 'extract_{}'.format(n))})
        return extract

    def set_output(self, event):
        for tool in self.tools:
            for result in os.listdir(tool.output_dir):
                source = os.path.join(tool.output_dir, result)
                target = os.path.join(self.output_dir, result)
                if os.path.exists(target):
                    os.remove(target)
                os.link(os.path.join(tool.output_dir, result), os.path.join(self.output_dir, result))
            self.logger.info('succeed in linking {} to {}'.format(source, target))
        self.end()

    @modlfuncdeco
    def end(self):
        super(ExtractSeqModule, self).end()