#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/27 14:57
@file    : gsea_exp_matrix.py
"""
import pandas as pd


class GseaExpMatrix(object):
    def __init__(self, exp_matrix_file, gene_detail_dic, out_path, geneset_type = "msigdb"):
        self._exp_matrix = exp_matrix_file
        self._gene_detail = gene_detail_dic
        self._out_path = out_path
        self.geneset_type = geneset_type

    def run(self):
        gene_desc = {}
        df = pd.read_table(self._exp_matrix, sep='\t', index_col=None)
        indexes = df['seq_id']


        gsymbols = []

        for g in indexes:
            print(g)
            if self.geneset_type == "msigdb":
                sub_dic = self._gene_detail.get(g)
                if not sub_dic:
                    gsymbols.append('__non__')
                    continue
                g_symbol = sub_dic.get('gene_name', '__non__').upper()
                g_desc = sub_dic.get('description', 'non')
                gsymbols.append(g_symbol)
                gene_desc[g_symbol] = g_desc
            else:
                g_symbol = g
                gsymbols.append(g_symbol)
                gene_desc[g_symbol] = g



        df['seq_id'] = gsymbols
        df = df.groupby('seq_id')
        df = df.max()
        df.index.name = 'geneid'
        if self.geneset_type == "msigdb":
            if '__non__' in list(df.columns):
                df = df.drop('__non__')
        samples = list(df.columns)
        df.insert(0, 'description', [gene_desc.get(i, 'non') for i in df.index])
        df.to_csv(self._out_path, sep='\t')
        return samples
