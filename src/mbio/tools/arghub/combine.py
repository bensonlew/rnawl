# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd
from collections import defaultdict


class CombineAgent(Agent):
    def __init__(self, parent):
        super(CombineAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': 'sample'},
            {'name': 'argout', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'trna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rnt', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'mge', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'mge_elem', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'elem_trna', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'elem_gene', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'mge_n', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'arg_n', 'type': 'outfile', 'format': 'sequence.profile_table'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('argout').is_set:
            raise OptionError('找不到输入文件')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'


class CombineTool(Tool):
    '''
    对arghub分析产生的结果文件结合数据注释文件进行格式化
    '''
    def __init__(self, config):
        super(CombineTool, self).__init__(config)
        self.id_loc = {}
        self.mge_args = {}  # 抗性基因和可移动元件的位置关系

    def run(self):
        super(CombineTool, self).run()
        try:
            self.out_for_arg()
        except Exception as e:
            self.set_error('wrong in out format {}'.format(e))
        if self.option('mge').is_set:
            self.combine_mge()
            self.combine_elem()
        self.end()

    def out_for_arg(self):
        arg_result = pd.read_csv(self.option('argout').path, sep='\t')
        self.ids = list(arg_result['gene_id'])
        self.logger.info('ids: {}'.format(self.ids))
        self.id_loc = {}
        nuc_dict = {}
        if self.option('cds').is_set:
            cds = self.get_seq(self.option('cds').path)
            nuc_dict.update(cds)
            #arg_result['nuc_seq'] = arg_result['gene_id'].apply(lambda x: cds[x])
        if self.option('rrna').is_set:
            if self.option('mge').is_set:
                rrna, id_loc = self.get_seq(self.option('rrna').path, True)
                self.id_loc.update(id_loc)
            else:
                rrna = self.get_seq(self.option('rrna').path)
            nuc_dict.update(rrna)
            #arg_result['nuc_seq'] = arg_result['gene_id'].apply(lambda x: rrna[x])
        if nuc_dict:
            arg_result['nuc_seq'] = arg_result['gene_id'].apply(lambda x: nuc_dict[x])

        if self.option('mge').is_set:
            faa, id_loc = self.get_seq(self.option('prot').path, True)
            self.id_loc.update(id_loc)
            arg_result['prot_seq'] = arg_result['gene_id'].apply(lambda x: faa[x] if x in faa else '-')
            arg_result = arg_result.rename(columns={"length": "align_len"})
            arg_result['location'] = arg_result['gene_id'].apply(lambda x: self.id_loc[x][0])
            arg_result['start'] = arg_result['gene_id'].apply(lambda x: self.id_loc[x][1])
            arg_result['end'] = arg_result['gene_id'].apply(lambda x: self.id_loc[x][2])
            arg_result['strand'] = arg_result['gene_id'].apply(lambda x: self.id_loc[x][3])
            arg_result['length'] = arg_result['end'] - arg_result['start'] + 1
            arg_result.apply(lambda x: self.id_loc[x['gene_id']].extend(
                [x['anti_class'], x['db_type']]), axis=1)
        elif self.option('prot').is_set:
            faa = self.get_seq(self.option('prot').path)
            arg_result['prot_seq'] = arg_result['gene_id'].apply(lambda x: faa[x] if x in faa else '-')
        arg_result.to_csv('arg_results.txt', sep='\t', index=False)
        self.option('arg_n', self.work_dir + '/arg_results.txt')

    def get_seq(self, seq, loc=False):
        self.logger.info('get seq:{} #'.format(seq))
        pre_name = ''
        pre_sequence = ''
        f = open(seq, 'r')
        seq_dict = {}
        seq_loc = {}
        name_p = re.compile(r'>(\S+)')
        while 1:
            l = f.readline().strip()
            if not l:
                break
            elif l.startswith('>'):
                this_name = name_p.search(l).groups()[0]
                self.logger.info('this_name:{}#'.format(this_name))
                if loc:
                    _, location, start, end, strand = l.split(' # ')
                    if int(start) > int(end):
                        start, end = end, start
                    if int(strand) == 1:
                        strand = '+'
                    elif int(strand) == -1:
                        strand = '-'
                self.logger.info('pre_name in self.ids:{}#'.format(pre_name in self.ids))
                if pre_name in self.ids:
                    seq_dict[pre_name] = pre_sequence
                    seq_loc[pre_name] = loc and [
                        location, int(start),
                        int(end), strand
                    ]
                pre_name = this_name
                pre_sequence = ''
            else:
                pre_sequence += l.strip()

        seq_dict[pre_name] = pre_sequence
        self.logger.info('seq_dict:{}#'.format(seq_dict))
        if loc:
            return seq_dict, seq_loc
        else:
            return seq_dict

    def combine_mge(self):
        mge = pd.read_csv(self.option('mge').path, sep='\t')
        mge_arglist = []
        for line in mge.itertuples():
            arglist = []
            for gene in self.id_loc:
                info = self.id_loc[gene]
                if gene in self.mge_args or line.location != info[0]:
                    continue
                if int(line.start) > int(line.end):
                    line.start, line.end = line.end, line.start
                gene_len = abs(info[1] - info[2])
                a = [int(line.start), int(line.end), info[1], info[2]]
                if max(a) - min(a) < line.length + gene_len:
                    self.mge_args[gene] = 'in'
                    arglist.append(gene)
            if arglist:
                mge_arglist.append(','.join(arglist))
            else:
                mge_arglist.append("--")
        mge['mge_arglist'] = mge_arglist
        mge.to_csv('mge.xls', sep='\t', index=False)
        self.option('mge_n', self.work_dir + '/mge.xls')

    def combine_elem(self):
        sample = self.option('sample')
        with open('elem_gene.txt', 'w') as elem_g:
            # sample	data_type	type	elem_name	location	start	end	strand	length
            elem_g.write('sample\ttype\tmelem_type\tmelem_name\tanti_type'
                         '\targ_loc\tlocation\tstart\tend\tstrand\tlength\n')
            for gene in self.id_loc:
                line = '{}\telem\tresistance gene\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                arg_loc = self.mge_args[gene] if gene in self.mge_args else 'out'
                print(gene)
                print(self.id_loc[gene])
                location, start, end, strand, anti_type, _ = self.id_loc[gene]
                elem_g.write(
                    line.format(sample, gene, anti_type, arg_loc, location,
                                start, end, strand,
                                abs(start - end) + 1))
        self.option('elem_gene', self.work_dir + '/elem_gene.txt')

        trna = defaultdict(str)
        with open(self.option('trna').path, 'r') as tr:
            s_id = ''
            for l in tr:
                l = l.strip()
                if l.startswith('>'):
                    s_id = l.split(' ')[0].split('>')[-1]
                else:
                    trna[s_id] += l

        with open(self.option('rnt').path,
                  'r') as rnt, open('elem_trna.txt', 'w') as elem_t:
            elem_t.write(
                'sample\ttype\tmelem_type\tmelem_name\tlocation\tstart\tend\tstrand\tlength\tseq\n'
            )
            line = sample + '\telem\ttRNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
            rnt.readline()
            for l in rnt:
                if 'tRNA' not in l:
                    continue
                l = l.strip().split('\t')
                start, end = l[0].split('..')
                elem_t.write(
                    line.format(l[4], l[3], start, end, l[1], l[2],
                                trna[l[4]]))
        self.option('elem_trna', self.work_dir + '/elem_trna.txt')
