# -*- coding: utf-8 -*-
# __author__ = "guhaidong"

import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SeqIO
import re

class Pocp(object):
    """
    进行pocp计算
    """
    def __init__(self):
        super(Pocp, self).__init__()
        self.seq_name = set()
        self.seqid2len = dict()
        self.pocp_count = 0
        self.hit_seq_name = set()
        self.identity = 40
        self.query_percent = 50

    def parse_fa(self, input):
        for seq in SeqIO.parse(input, "fasta"):
            self.seq_name.add(seq.id)
            self.seqid2len[seq.id] = len(seq)

    def parse_table(self, table, table_type):
        if table_type == "m5":
            self.parse_m5(table)
        elif table_type == "m6":
            self.parse_m6(table)

    def parse_m5(self, table):
        records = NCBIXML.parse(open(table))
        for rec in records:
            query = re.split(' ', rec.query, maxsplit=1)[0]
            query_len = rec.query_length
            for align in rec.alignments:
                if query in self.hit_seq_name:
                    break
                for hsp in align.hsps:
                    if query in self.hit_seq_name:
                        break
                    if self.seqid2len:
                        query_align_percent = float(hsp.align_length) / self.seqid2len[query] * 100
                    else:
                        query_align_percent = float(hsp.align_length) / query_len * 100
                    if float(hsp.identities) >= self.identity and query_align_percent >= self.query_percent:
                        self.pocp_count += 1
                        self.hit_seq_name.add(query)

    def parse_m6(self, table):
        pass  # 先处理m5

    def get_pocp_value(self):
        pocp_value = float(self.pocp_count) / len(self.seq_name) * 100
        return pocp_value

    # def export(self, output_table):
    #     pocp_result = float(self.pocp_count) / len(self.seq_name)
    #     return pocp_result