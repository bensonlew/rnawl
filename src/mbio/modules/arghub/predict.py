# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.module import Module
from biocluster.option import OptionError


class PredictModule(Module):
    def __init__(self, parent):
        super(PredictModule, self).__init__(parent)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'mode', 'type': 'string', 'default': 'single'},
            {'name': 'gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'cds', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'r_gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'trna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 't_gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
        ]
        self.add_option(options)
        self.gene = self.add_tool('arghub.prodigal')
        self.rrna = self.add_tool('arghub.barrnap')
        self.trna = self.add_tool('predict.trnascanse')
        self.tools = []

    def check_options(self):
        if not self.option('genome').is_set:
            raise OptionError('必须设置输入 genome 参数')
        self.genome = self.option('genome').prop['path']

    def run(self):
        super(PredictModule, self).run()
        self.run_set(self.gene, opts={'input': self.genome, 'mode': self.option('mode')})
        self.run_set(self.trna, opts={'input_genome': self.genome})
        self.run_set(self.rrna, opts={'input_genome': self.genome})
        self.on_rely(self.tools, self.set_output)

    def run_set(self, tool, opts={}):
        self.tools.append(tool)
        tool.set_options(opts)
        tool.run()

    def set_output(self):
        self.logger.info('start set output {}'.format(self.tools))
        for tool in self.tools:
            self.logger.info('tool {}'.format(tool))
            for f in os.listdir(tool.output_dir):
                self.link(os.path.join(tool.output_dir, f), cover=True)
                if 'prot.faa' in f:
                    self.option('prot', os.path.join(tool.output_dir, f))
                if 'rRNA.fna' in f:
                    self.option('rrna', os.path.join(tool.output_dir, f))
                if 'tRNA.fnn' in f and os.path.getsize(os.path.join(tool.output_dir, f)) > 0:
                    self.option('trna', os.path.join(tool.output_dir, f))
                if 'cds.fna' in f:
                    self.option('cds', os.path.join(tool.output_dir, f))
                if '.gff' in f and os.path.getsize(os.path.join(tool.output_dir, f)) > 0:
                    if 'rRNA.gff' in f:
                        k = 'r_gff'
                    elif 'tRNA.gff' in f:
                        k = 't_gff'
                    elif 'tRNA.gbk.gff' in f:
                        continue
                    else:
                        k = 'gff'
                    self.logger.info('#### {} {} ###'.format(k, f))
                    self.option(k, os.path.join(tool.output_dir, f))
        self.end()
