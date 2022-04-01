# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys
import os


class DNAInfos(object):
    def __init__(self, fasta_file):
        self.seqs = SeqIO.parse(fasta_file, 'fasta', generic_dna)
        self._gc_skew = []
        self._gc_content = []
        self._seq_len = []

    @property
    def gc_skew(self):
        return self._gc_skew

    @property
    def gc_content(self):
        return self._gc_content

    @property
    def seq_len(self):
        return self._seq_len

    def parse_gc(self, window, step):
        for seq in self.seqs:
            self._seq_len.append([seq.id, str(len(seq))])
            chucks = list(range(0, len(seq), step))
            for i in chucks[:-1]:
                window_seq = seq[i:i+window].upper()
                self.__gc_infos(window_seq, i)
            wds = seq[chucks[-1]:].upper()
            self.__gc_infos(wds, chucks[-1])

    def __gc_infos(self, seq, start):
        g_count = seq.seq.count('G')
        c_count = seq.seq.count('C')
        gc_content = float(g_count + c_count) / len(seq)
        if int(g_count + c_count) ==0:
            gc_skew = 0
        else:
            gc_skew = float(g_count-c_count) / (g_count + c_count)
        self._gc_skew.append([seq.id, str(start + 1), str(start + len(seq)), str(gc_skew)])
        self._gc_content.append([seq.id, str(start + 1), str(start + len(seq)), str(gc_content)])


def output_list(list, out_name):
    with open(out_name, 'w') as w:
        [w.write('\t'.join(l) + '\n') for l in list]


if __name__ == "__main__":
    '''
    给定窗口大小和步长计算基因组的gc含量和gc skew值
    '''
    if len(sys.argv) == 1:
        exit('参数设置：\n python gc_infos.py DNAfasta [window] [step]')
    fasta_file = sys.argv[1]
    window = sys.argv[2] if len(sys.argv) > 2 else 10000
    step = sys.argv[3] if len(sys.argv) == 4 else 10000

    seqs = DNAInfos(fasta_file)
    seqs.parse_gc(int(window), int(step))
    gc_skew = seqs.gc_skew
    gc_content = seqs.gc_content
    seq_len = seqs.seq_len

    prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    output_list(gc_skew, prefix + '.gc_skew.xls')
    output_list(gc_content, prefix + '.gc_content.xls')
    output_list(seq_len, prefix + '.seq_len.xls')
