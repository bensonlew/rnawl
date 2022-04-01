# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os


class InputPrepareAgent(Agent):
    def __init__(self, parent):
        super(InputPrepareAgent, self).__init__(parent)
        options = [
            {'name': 'gff', 'type': 'infile', 'format': 'gene_structure.gff3',},
            {'name': 'genome', 'type': 'infile', 'format': 'sequence.fasta',},
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta',},
            {'name': 'cds_o', 'type': 'outfile', 'format': 'sequence.fasta',},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'trna', 'type': 'outfile', 'format': 'sequence.fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info([self.option('gff').is_set, self.option('cds').is_set, self.option('genome').is_set])
        if not any([self.option('gff').is_set, self.option('cds').is_set, self.option('genome').is_set]):
            raise OptionError('请设置输入 cds 文件或 gff + genome 文件')
        if not self.option('cds').is_set and (not all([self.option('gff').is_set, self.option('genome').is_set])):
            raise OptionError('必须同时设置gff 文件和 genome文件')
        return True

    def set_resource(self):
        self._memory = '2G'
        self._cpu = 2


class InputPrepareTool(Tool):
    def __init__(self, config):
        super(InputPrepareTool, self).__init__(config)
        self.trans = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.gff_read = "bioinfo/WGS/gffread/gffread-master/gffread"

    def run(self):
        super(InputPrepareTool, self).run()
        if self.option('cds').is_set:
            self.cds = self.option('cds').path
            self.prefix = os.path.splitext(os.path.basename(self.cds))[0]
            self.run_trans()
        else:
            self.gff = self.option('gff').prop['path']
            self.genome = self.option('genome').prop['path']
            self.prefix = os.path.splitext(os.path.basename(self.genome))[0]
            self.run_extract_trans()
        self.set_output()
        self.end()

    def run_trans(self):
        cmd = '{} -trim -table 11 -sequence {} -outseq {}.prot'
        cmd = cmd.format(self.trans, self.cds, self.prefix)
        command = self.add_command('trans_cds', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('trans cds to prot done!')
            cmd2 = self.config.SOFTWARE_DIR + "/program/perl-5.24.0/bin/perl -i -lpe 's/\_\d+$//' " + self.prefix + '.prot'
            command2 = self.add_command('transed_cds', cmd2, shell=True).run()
            self.wait(command2)
        else:
            self.set_error('wrong in trans cds to prot')

    def run_extract_trans(self):
        cmd1 = '{} {} -g {} -x {}.cds -o tmp.txt'.format(self.gff_read, self.gff,
                                                             self.genome,
                                                             self.prefix)
        command = self.add_command('extract_cds', cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('extract cds through gff done!')
            self.cds = self.work_dir + '/' + self.prefix + '.cds'
            self.run_trans()
        else:
            self.set_error('wrong in extract cds')
        # 提取非编码序列
        cmd2 = '{} {} -g {} -F --nc -w {}.rna -o tmp.txt'.format(self.gff_read, self.gff,
                                                                     self.genome,
                                                                     self.prefix)
        command2 = self.add_command('extract_non_coding', cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info('extract non-coding sequence done')
            self.r_t_rna()
        else:
            self.set_error('wrong in extracting non-coding sequence')

    def r_t_rna(self):
        with open(self.prefix + '.rna', 'r') as r, open(self.prefix + '.trna', 'w') as trna, open(self.prefix + '.rrna', 'w') as rrna:
            for l in r:
                if l.startswith('>') and 'trna' in l.lower():
                    tag = 'trna'
                elif l.startswith('>'):
                    tag = ''
                if tag == 'trna':
                    trna.write(l)
                else:
                    rrna.write(l)

    def set_output(self):
        for f in ['prot', 'cds', 'rrna', 'trna']:
            file_path = os.path.join(self.work_dir, self.prefix + '.' + f)
            if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
                print self.option(f).path
                if f == 'cds':
                    f += '_o'
                self.option(f).set_path(file_path)
                print self.option(f).path
        if self.option('gff').is_set and not self.option('cds_o').is_set:
            self.set_error('设置gff文件参数后, 必须在结果中配置 cds')

