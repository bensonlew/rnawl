# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.module import Module
from biocluster.option import OptionError


class ArgMgeModule(Module):
    def __init__(self, parent):
        super(ArgMgeModule, self).__init__(parent)
        options = [
            {'name': 'file_for_mge', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'genome', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'genome_type', 'type': 'string', 'default': 'single'},
            {'name': 'name_file_list', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'cds', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'trna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'rnt', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'rename', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'mge_result', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'mge_elem', 'type': 'outfile', 'format': 'sequence.profile_table'},
        ]
        self.add_option(options)
        self.mge_file = self.add_tool('arghub.cat_file')
        self.rename = self.add_tool('arghub.rename')
        self.gff = ''
        self.cds = ''
        self.prot = ''
        self.trna = ''

    def check_options(self):
        if not self.option('genome').is_set:
            raise OptionError('缺少参数genome')
        self.logger.info('seq_num {}'.format(self.option('genome').prop["seq_number"]))
        n = 0
        with open(self.option('genome').path, 'r') as r:
            for l in r:
                if l.startswith('>'):
                    n+=1
                if n > 50:
                    break
        if n <= 50:
            self.mge = self.add_module('mobile_genetic_elements.mobile_elements')  # 可移动元件模块
        else:
            self.mge = self.add_module('mobile_genetic_elements.meta_mobile_elements')  # 可移动元件模块

    def run(self):
        super(ArgMgeModule, self).run()
        self.sample = self.option('genome').path.split('/')[-1].rpartition('.')[0]
        if self.option('name_file_list').is_set:
            self.mge_file.on('end', self.run_rename)
            self.rename.on('end', self.run_mge)
        else:
            self.mge_file.on('end', self.run_mge)
        self.mge.on('end', self.set_output)
        #self.rename.on('end', self.set_output)

        self.run_mge_file()

    def run_mge_file(self):
        opts = {
            'file_list': self.option('file_for_mge').path,
            'sample_name': self.sample
        }
        self.mge_file.set_options(opts)
        self.mge_file.run()

    def run_rename(self):
        opts = {
            'cds': self.mge_file.option('cds').path,
            'prot': self.mge_file.option('prot').path,
            'gff': self.mge_file.option('gff').path,
            'file_list': self.option('name_file_list').path
        }
        self.rename.set_options(opts)
        self.rename.run()

    def run_mge(self):
        opts = {
            'sample': self.sample,
            'genome_fa': self.option('genome').path,
            'rnt': self.mge_file.option('rnt').path
        }
        if self.option('name_file_list').is_set:
            self.gff = self.rename.option('rename_gff').path
            self.prot = self.rename.option('rename_prot').path
            self.cds = self.rename.option('rename_cds').path
        else:
            self.gff = self.mge_file.option('gff').path
            self.prot = self.mge_file.option('prot').path
            self.cds = self.mge_file.option('cds').path

        opts.update({
            'ptt': self.gff,
            'gene_fna': self.cds,
            'faa': self.prot,
        })
        self.logger.info(opts)
        self.mge.set_options(opts)
        self.mge.run()

    def set_output(self):
        self.option('cds', self.cds)
        self.option('prot', self.prot)
        self.option('gff', self.gff)
        self.option('trna', self.mge_file.option('trna').path)
        self.option('rnt', self.mge_file.option('rnt').path)
        self.option('rrna', self.mge_file.option('rrna').path)
        if os.path.exists(os.path.join(self.mge.output_dir, self.sample + '.mge.xls')):
            self.option('mge_result', os.path.join(self.mge.output_dir, self.sample + '.mge.xls'))
        if os.path.exists(os.path.join(self.mge.output_dir, self.sample + '.element.xls')):
            self.option('mge_elem', os.path.join(self.mge.output_dir, self.sample + '.element.xls'))
        if self.option('name_file_list').is_set:
            self.option('rename', os.path.join(self.rename.work_dir, 'names.txt'))
        self.end()
