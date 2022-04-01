# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import re
from Bio import SeqIO


class UniFormatAgent(Agent):
    def __init__(self, parent):
        super(UniFormatAgent, self).__init__(parent)
        opts = [
            {'name': 'gff', 'type': 'infile', 'format': 'gene_structure.gff3'},
            {'name': 'gff_o', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 't_gff', 'type': 'infile', 'format': 'gene_structure.gff3'},
            {'name': 'r_gff', 'type': 'infile', 'format': 'gene_structure.gff3'},
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'prot_o', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'trna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'trna_o', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'rrna_o', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'cds_o', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rnt', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'predict', 'type': 'int', 'default': 0},  # 1 or 0, 1表示预测结果
            {'name': 'name_file', 'type': 'outfile', 'format': 'sequence.profile_table' },  # 保存有基因名称的列表用于改名用
        ]
        self.add_option(opts)

    def check_options(self):
        if len(self.get_option_object()) < 2:
            raise OptionError('指定输入文件')

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'


class UniFormatTool(Tool):
    def __init__(self, config):
        super(UniFormatTool, self).__init__(config)
        os.environ['OPENBLAS_NUM_THREADS'] = '1'

    def run(self):
        super(UniFormatTool, self).run()
        if self.option('t_gff').is_set:
            self.logger.info("###")
            self.logger.info(self.option('t_gff').path)
            self.logger.info(os.path.getsize(self.option('t_gff').path))
            self.logger.info("###")
            self.rm_overlap()
        id2index = self.re_gff()
        if self.option('predict') == 1:
            self.name_file(id2index)
        self.uf_coding()
        self.uf_rna()
        self.logger.info(self.option('name_file'))
        for f in self.get_option_object():
            print f
            if hasattr(self.option(f), 'path'):
                print self.option(f).path
        self.end()

    def rm_overlap(self):
        '''
        去除 tRNA和gene中有重复的部分
        '''
        import pandas as pd
        bedtool = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        gene_gff = pd.read_csv(self.option('gff').path, comment='#', sep='\t', header=None)
        print gene_gff.head()
        self.logger.info("###2")
        self.logger.info(self.option('t_gff').path)
        self.logger.info(os.path.getsize(self.option('t_gff').path))
        self.logger.info("###2")
        with open(self.option('t_gff').path, 'r') as tt, open('trna.gff', 'w') as temp_t:
            line = '{}\t{}\ttRNA\t{}\t{}\t.\t{}\t0\tID={};Product={};\n'
            tt.readline()
            for l in tt:
                l = l.strip().split('\t')
                strand = '+'
                scaf = l[0].rpartition('_')[0]
                scaf_1 = l[1].rpartition('_')[0]
                if scaf:
                    scaf += '_' + scaf_1
                else:
                    scaf = scaf_1
                if int(l[2]) > int(l[3]):
                    strand = '-'
                    l[2], l[3] = l[3], l[2]
                filter_1 = gene_gff[0] == scaf
                filter_2 = gene_gff[3] < int(l[2])
                filter_3 = gene_gff[4] > int(l[2])
                filter_4 = gene_gff[3] < int(l[3])
                filter_5 = gene_gff[4] > int(l[3])
                f1 = filter_1 & filter_2 & filter_3
                f2 = filter_1 & filter_4 & filter_5
                if not (any(f1) or any(f2)):
                    temp_t.write(
                        line.format(scaf, l[1], l[2], l[3], strand, l[1], l[4])
                    )
        self.option('t_gff', os.path.join(self.work_dir, 'trna.gff'))

    def name_file(self, id2index):
        '''
        用于生成改名相关的记录格式 sample\told_name\tnew_name
        '''
        sample = os.path.basename(self.option('gff').path).rpartition('.gff')[0]
        with open('name_list.txt', 'w') as w:
            [w.write('{}\n'.format(k)) for k in id2index]
        self.option('name_file', os.path.join(self.work_dir, 'name_list.txt'))

    def re_gff(self):
        '''
        修改gff文件格式，用于可移动元件分析
        '''
        file_name = os.path.basename(self.option('gff').path)
        out_file = os.path.join(self.output_dir, file_name)
        id2index = {}
        if self.option('predict') == 1:
            partner = re.compile(r'^(\S+).*ID\=\d+\_([^\;]+)')  #用于匹配id
        else:
            partner = re.compile(r'^\S+.*ID\=([^\;]+)')  #用于匹配id
        sub1 = re.compile(r'ID\=[^\;]+')  # 用于替换编号 else:
        with open(out_file, 'w') as out, open(self.option('gff').path, 'r') as in_f:
            i = 0
            for li in in_f:
                match = partner.search(li)
                if not match:
                    out.write(li)
                    continue
                g_id = '_'.join(match.groups())
                id2index[g_id] = g_id
                li = sub1.sub('ID=' + g_id, li).split('\t')
                li[1] = g_id
                out.write('\t'.join(li))
        self.option('gff_o', out_file)
        self.option('gff').set_path(out_file)
        if self.option("rrna").is_set:
            self.option('rrna_o', self.option("rrna").path)
        return id2index

    def uf_coding(self):
        '''
        每个序列id后面加上 起始 终止 链格式为 >id # scaf  # start # end # strand
        '''
        seq_header = {}
        partner = re.compile(r'ID\=([^\;]+)')
        with open(self.option('gff').path, 'r') as gff:
            for l in gff:
                s_id = partner.search(l)
                if not s_id:
                    continue
                else:
                    s_id = s_id.groups()[0]
                l = l.strip().split('\t')
                scaf, start, end = l[0], l[3], l[4]
                if l[6] == '+':
                    strand = 1
                else:
                    strand = -1
                if start > end:
                    start, end = end, start
                seq_header[s_id] = '{} # {} # {} # {} # {}'.format(s_id, scaf, start, end, strand)
        for f in ['cds', 'prot']:
            f_path = self.option(f).path
            f_fasta = SeqIO.parse(f_path, 'fasta')
            f_name = os.path.basename(f_path)
            f_out = os.path.join(self.output_dir, f_name)
            with open(f_out, 'w') as out:
                for o in f_fasta:
                    if o.id in seq_header:
                        out.write('>{}\n{}\n'.format(seq_header[o.id], o.seq))
            self.option(f + '_o',  f_out)

    def uf_rna(self):
        '''
        构建meg分析需要用到的rna相关信息文件
        '''
        out_header = '\t'.join([
            'Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product'
        ])
        r_out = open(os.path.join(self.output_dir, 'rna.rnt'), 'w')
        r_out.write(out_header + '\n')
        if self.option('t_gff').is_set:
            # 输出 .rnt 格式
            p_p = re.compile(r'Product\=(\S+)\;')
            with open(self.option('t_gff').prop['path'], 'r') as t_gff:
                for li in t_gff:
                    li = li.strip().split('\t')
                    scaf, s_id, start, end, strand = li[0], li[1], li[3], li[4], li[6]
                    product = p_p.search(li[-1]).groups()[0]
                    r_out.write("{}..{}\t{}\t{}\t{}\t{}\t-\t-\t-\ttRNA {}\n".format(
                        start, end, strand, int(end) - int(start) + 1, scaf, s_id, product
                    ))
            # 修改trna序列文件id，跟gff文件一致
            trna_out = open(os.path.join(self.output_dir, 'trna.fnn'), 'w')
            with open(self.option('trna').path, 'r') as trna:
                for t in trna:
                    if t.startswith('>'):
                        trna_out.write('>{}\n'.format(t.split()[3]))
                    else:
                        trna_out.write(t)
            trna_out.close()
            self.option('trna_o', os.path.join(self.output_dir, 'trna.fnn'))
        if self.option('r_gff').is_set or self.option('gff'):
            # 输出.rnt格式
            if self.option('r_gff').is_set:
                gff_path = self.option('r_gff').prop['path'] 
            elif self.option('gff').is_set:
                gff_path = self.option('gff').prop['path'] 
            head_dict = {}
            with open(gff_path, 'r') as gff:
                for li in gff:
                    if li.startswith('#'):
                        continue
                    li = li.strip().split('\t')
                    if li[2] not in ['rRNA', 'tRNA'] or li[3].isalpha():
                        continue
                    scaf, start, end, strand = li[0], li[3], li[4], li[6]
                    if int(start) > int(end):
                        start, end = end, start
                    product = re.search(r'product\=([^\;]+)', li[-1]).groups()[0]
                    product = re.sub(r'-', ' ', product)
                    s_id = re.search(r'(?:Name|ID)\=([^\;]+)', li[-1]).groups()[0]
                    r_out.write("{}..{}\t{}\t{}\t{}\t{}\t-\t-\t-\t{}\n".format(
                        start, end, strand, int(end) - int(start) + 1, scaf, s_id, product
                    ))
                    head_dict[s_id] = [">{} # {} # {} # {} # ".format(s_id, scaf, start, end), strand]
            if self.option('rrna').is_set:
                with open(self.option('rrna').path, 'r') as r, open(os.path.join(self.output_dir, 'rrna.fnn'), 'w') as rrna_out:
                    for l in r:
                        if l.startswith('>'):
                            s_id = l.strip().split()[0].split('>')[-1]
                            header, strand = head_dict[s_id]
                            if strand == '+':
                                rrna_out.write("{}1\n".format(header))
                            else:
                                rrna_out.write("{}-1\n".format(header))
                        else:
                            rrna_out.write(l)
                self.option('rrna_o', os.path.join(self.output_dir, 'rrna.fnn'))
        self.option('rnt', os.path.join(self.output_dir, 'rna.rnt'))
