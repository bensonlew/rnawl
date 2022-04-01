# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.module import Module
from biocluster.option import OptionError


class ArghubModule(Module):
    def __init__(self, parent):
        super(ArghubModule, self).__init__(parent)
        options = [
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'read1', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read2', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'aligner', 'type': 'string', 'default': 'blast'},
            {'name': 'db_type', 'type': 'string', 'default': 'core'},
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'identity', 'type': 'float', 'default': 50},
            {'name': 'nthread', 'type': 'int', 'default': 8},
            {'name': 'output', 'type': 'string'}
        ]
        self.add_option(options)
        self.tools = []

    def check_options(self):
        a = [
            self.option('prot').is_set, self.option('rrna').is_set,
            self.option('read1').is_set, self.option('read2').is_set
        ]
        if not any(a):
            raise OptionError('需要输入蛋白序列、rRNA序列或双端reads')
        elif any(a[2:]) and (not all(a[2:])):
            raise OptionError('reads输入必须为双端reads')

    def run(self):
        super(ArghubModule, self).run()
        self.set_run()
        self.on_rely(self.tools, self.set_output)
        [tool.run() for tool in self.tools]
    
    def set_run(self):
        if self.option('read1').is_set:
            ariba = self.add_tool('arghub.ariba')
            self.tools.append(ariba)
            self.set_ariba(ariba)
        else:
            if self.option('prot').is_set:
                align1 = self.add_tool('arghub.align')
                self.tools.append(align1)
                self.set_align(align1, self.option('prot').prop['path'], 'blastp')
            if self.option('rrna').is_set:
                align2 = self.add_tool('arghub.align')
                self.tools.append(align2)
                self.set_align(align2, self.option('rrna').prop['path'], 'blastn')

    def set_ariba(self, ariba):
        opts = {
            'read1': self.option('read1').path,
            'read2': self.option('read2').path,
            'nthread': self.option('nthread'),
            'db_type': self.option('db_type'),
        }
        ariba.set_options(opts)

    def set_align(self, tool, input_seq, blast):
        opts = {
            'input': input_seq,
            'blast': blast,
            'aligner': self.option('aligner'),
            'evalue': self.option('evalue'),
            'identity': self.option('identity'),
            'nthread': self.option('nthread'),
            'db_type': self.option('db_type'),
            'out_prefix': os.path.basename(input_seq)
        }
        tool.set_options(opts)

    def set_output(self):
        o = []
        for t in self.tools:
            self.link(t.option('output').path, cover=True)
            o.append(t.option('output').path)
        self.option('output', '\n'.join(o))
        self.end()
