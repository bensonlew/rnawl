# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.module import Module
from biocluster.option import OptionError
from biocluster.iofile import FileBase


class OutfmtModule(Module):
    def __init__(self, parent):
        super(OutfmtModule, self).__init__(parent)
        options = [
            {'name': 'argout_list', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'arg_type', 'type': 'string'},
            {'name': 'sample', 'type': 'string'},
            {'name': 'rename', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'prot_list', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'trna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rnt', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'mge_result', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'mge_elem', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'elem_trna', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'elem_gene', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'mge_n', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'arg_n', 'type': 'outfile', 'format': 'sequence.profile_table'},
        ]
        self.add_option(options)
        self.argout = self.add_tool('arghub.argout')
        self.combine = self.add_tool('arghub.combine')

    def check_options(self):
        if self.option("prot_list").is_set:
            prot_path = os.path.join(self.work_dir, "all_prot.faa")
            if os.path.exists(prot_path):
                os.remove(prot_path)
            for file_path in self.option("prot_list").file_list:
                cmd = "cat {} >> {}".format(file_path, prot_path)
                os.system(cmd)
            # self.option("prot").set_path(prot_path)
        if not (self.option("prot").is_set or self.option("rrna").is_set):
            raise OptionError("prot 或 rrna  必须设置")

    def run(self):
        super(OutfmtModule, self).run()
        self.on_rely([self.argout, self.combine], self.set_output)
        self.argout.on('end', self.run_combine)
        self.run_argout()

    def run_argout(self):
        opts = {
            'file_list': self.option('argout_list'),
            'file_type': self.option('arg_type'),
        }
        if self.option('rename').is_set:
            opts['rename'] = self.option('rename')
        self.argout.set_options(opts)
        self.argout.run()

    def run_combine(self):
        opts = {
            'argout': self.argout.option('output'),
        }
        if self.option('arg_type') == 'cds':
            opts['prot'] = self.option('prot')
            opts['cds'] = self.option('cds')
        elif self.option('arg_type') == 'prot':
            opts['prot'] = self.option('prot')
        elif self.option('arg_type') == 'rrna':
            opts['rrna'] = self.option('rrna')
        elif self.option('arg_type') in ['single', 'meta']:
            opts.update({
                'prot': self.option('prot'),
                'sample': self.option('sample'),
                'cds': self.option('cds'),
                'trna': self.option('trna'),
                'rrna': self.option('rrna'),
                'rnt': self.option('rnt'),
                'mge': self.option('mge_result'),
                'mge_elem': self.option('mge_elem'),
            })
        self.logger.info('{} {}'.format(self.option('arg_type'), opts))
        self.logger.info("self.option('prot') {} {}".format(self.option('prot'), isinstance(self.option('prot'), FileBase)))
        self.combine.set_options(opts)
        self.combine.run()

    def set_output(self):
        self.option('arg_n', self.combine.option('arg_n').path)
        if self.option('arg_type') in ['single', 'meta']:
            self.option('mge_n', self.combine.option('mge_n').path)
            self.option('elem_trna', self.combine.option('elem_trna').path)
            self.option('elem_gene', self.combine.option('elem_gene').path)
        self.end()
