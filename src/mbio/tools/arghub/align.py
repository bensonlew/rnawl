# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import pandas as pd
import os


class AlignAgent(Agent):
    def __init__(self, parent):
        super(AlignAgent, self).__init__(parent)
        options = [
            {'name': 'aligner', 'type': 'string', 'default': 'blast'}, # blast, hmmscan, diamond
            {'name': 'blast', 'type': 'string', 'default': 'blastp'},
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'identity', 'type': 'float', 'default': 50},
            {'name': 'nthread', 'type': 'int', 'default': 4},
            {'name': 'input', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'db_type', 'type': 'string', 'default': 'core'}, # core, plus
            {'name': 'out_prefix', 'type': 'string', 'default': 'out'},
            {'name': 'output', 'type': 'outfile', 'format': 'sequence.profile_table'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('input').is_set:
            raise OptionError('缺少输入序列文件')

    def set_resource(self):
        self._memory = '5G'
        self._cpu = self.option('nthread')


class AlignTool(Tool):
    def __init__(self, config):
        super(AlignTool, self).__init__(config)
        self.db_path = self.config.SOFTWARE_DIR + '/database/ArgHub/'
        self.db_path += self.option('aligner') + '/' + self.option('db_type')
        self.self_score = self.config.SOFTWARE_DIR +\
            '/database/ArgHub/self_score.' + self.option('db_type')
        self.blast_bin = '/bioinfo/align/ncbi-blast-2.3.0+/bin/'
        self.diamond = '/bioinfo/align/diamond-0.9.11/diamond'
        self.hmmscan = '/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan'
        self.hmmdb = self.config.SOFTWARE_DIR + '/database/ArgHub/hmmscan/arghub.hmm'

    def run(self):
        super(AlignTool, self).run()
        if self.option('aligner') == 'hmmscan':
            self.run_hmmscan()
        else:
            self.run_align()
        self.set_output()
        self.end()

    def run_hmmscan(self):
        # hmmscan
        cmd = self.hmmscan + ' --tblout {0}.tblout --domtblout {0}.domtblout -o {0}.stdout' +\
            ' -E 1e-5 --domE 1e-3 --cpu {1} {2} {3}'
        cmd = cmd.format(self.option('out_prefix'),
                         self.option('nthread'), self.hmmdb,
                         self.option('input').path)
        if os.path.exists('hmmscan.o'):
            self.hmm_out()
        else:
            command = self.add_command('hmmscan', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info('hmmscan done!')
                self.hmm_out()
            else:
                self.set_error('wrong in hmmscan')

    def hmm_out(self):
        out = self.option('out_prefix') + '.domtblout'
        cmd = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin/perl '
        cmd += self.config.PACKAGE_DIR + '/arghub/hmm_out.pl {2}/{0} > {2}/{1}.hmmout'
        cmd = cmd.format(out, self.option('out_prefix'), self.work_dir)
        command = self.add_command('hmm_out', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            if os.system("wc -l {}/{}.hmmout".format(self.work_dir, self.option('out_prefix'))) == 1:
                self.set_error("You input have no result!")
            self.logger.info('hmm out done')
        else:
            self.set_error('wrong hmm_out')

    def run_align(self):
        # blast/diamond
        outfmt = '6 qseqid sseqid pident length mismatch ' +\
            'gapopen qstart qend sstart send evalue bitscore qlen slen gaps'
        if self.option('blast') == 'blastp':
            db = self.db_path + '_prot'
        else:
            db = self.db_path + '_rrna'
            db = db.replace('plus', 'core')
        if self.option('aligner').lower() == 'blast':
            cmd = self.blast_bin + "{} -outfmt '{}' -max_target_seqs 1 -num_threads {} " +\
                '-query {} -db {} -evalue {} -out {}.{}'
        else:
            cmd = self.diamond + ' {}  -f {} -p {} -k 1 ' +\
                '-q {} -d {} -e {} -o {}.{}'
        cmd = cmd.format(self.option('blast'), outfmt, self.option('nthread'),
                         self.option('input').prop['path'], db, self.option('evalue'),
                         self.option('out_prefix'), self.option('blast'))
        command = self.add_command('align', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            if os.path.getsize('{}.{}'.format(self.option('out_prefix'), self.option('blast'))) == 0:
                self.set_error("You input have no result!")
            self.logger.info('alignment done! now, evaluate the result...')
            self.out_filter('{}.{}'.format(self.option('out_prefix'), self.option('blast')))
            try:
                pass #self.out_filter('{}.{}'.format(self.option('out_prefix'), self.option('blast')))
            except Exception as e:
                self.set_error('wrong in evaluation: {}'.format(e))
            self.logger.info('evaluation done')
        else:
            self.set_error('wrong in alignment')

    def out_filter(self, blast_out):
        scores = pd.read_csv(self.self_score, sep='\t')
        scores = dict(zip(scores['arghub_id'], scores['ref_score']))
        cols = 'qseqid sseqid pident length mismatch ' +\
            'gapopen qstart qend sstart send evalue bitscore qlen slen gaps'
        table = pd.read_csv(blast_out, header=None, sep='\t', names=cols.split())
        table = table[table['pident'] > self.option('identity')]

        table['evaluation'] = table.apply(self.evaluation, axis=1, args=(scores,))
        os.rename(blast_out, blast_out + '_ori')
        table = table.rename(columns={'sseqid': 'arghub_id', 'qseqid': 'gene_id',
                                      'pident': 'identity', 'bitscore': 'score'})
        table.drop_duplicates('gene_id', inplace=True)
        table.to_csv(blast_out, index=False, sep='\t')

    def evaluation(self, row, scores):
        evaluation = 'questionable'
        cov = float(row['length'] - row['gaps']) / row['slen']
        ref_score = scores[row['sseqid']]
        print(cov)
        print(row['gaps'])
        print(row['pident'])
        if cov == 1 and row['pident'] == 100 and row['gaps'] == 0:
            evaluation = 'perfect'
        elif row['bitscore'] >= 0.9 * ref_score:
            evaluation = 'high'
        elif row['pident'] >= 80 and cov >= 0.8:
            evaluation = 'moderate'
        elif row['pident'] >= 80 or cov >= 0.8:
            evaluation = 'low'
        return evaluation

    def set_output(self):
        self.logger.info('start set output')
        if self.option('aligner') == 'hmmscan':
            out = self.option('out_prefix') + '.hmmout'
        else:
            out = '{}.{}'.format(self.option('out_prefix'), self.option('blast'))
        out_file = os.path.join(self.output_dir, out)
        if os.path.exists(out_file):
            os.remove(out_file)
        os.link(out, out_file)
        self.option('output', out_file)

