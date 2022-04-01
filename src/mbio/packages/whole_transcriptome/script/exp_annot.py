# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging

import numpy as np
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

t_tpm_matrix = 'T.tpm.txt'
g_tpm_matrix = 'G.tpm.txt'
annot_matrix = 'all_annot.xls'
relation_map = 'trans_type.xls'


t_tpm_df = pd.read_table(t_tpm_matrix, index_col=0)
g_tpm_df = pd.read_table(g_tpm_matrix, index_col=0)
annot_df = pd.read_table(annot_matrix)
relat_df = pd.read_table(relation_map, names=['transcript_id', 'gene_id', 'category', 'kind'], index_col=0)

t_annot_df = annot_df.set_index('transcript_id')
g_annot_df = annot_df.set_index('gene_id', 'is_gene')

logging.info('succeed in loading data')

def get_annot_dict(pd_obj, level='T'):
    if level == 'T':
        assert isinstance(pd_obj, pd.core.series.Series)
        annot_dict = pd_obj.rename(str.lower).to_dict()
    elif level == 'G':
        assert isinstance(pd_obj, pd.core.series.Series)
        annot_dict = pd_obj.rename(str.lower).to_dict()
    return annot_dict

mrna_record_dict_list = list()
lncrna_record_dict_list = list()
mrna_gene_record_dict_list = list()
lncrna_gene_record_dict_list = list()

count = 0
sample_list = list(t_tpm_df.columns[1:-2])

for idx, row in t_tpm_df.iterrows():
    record_dict = row.rename({'rna_type': 'category'}).to_dict()
    record_dict['transcript_id'] = idx
    record_dict['gene_id'] = relat_df.loc[idx]['gene_id']
    record_dict['kind'] = 'new' if record_dict.pop('is_new') else 'ref'
    record_dict['level'] = 'T'
    if record_dict['category'] == 'mRNA':
        pd_obj = t_annot_df.loc[idx]
        annot_dict = get_annot_dict(pd_obj)
        record_dict.update(annot_dict)
        mrna_record_dict_list.append(record_dict)
    elif record_dict['category'] == 'lncRNA':
        lncrna_record_dict_list.append(record_dict)
    count += 1
    if not count % 10000:
        logging.debug('succeed in processing {} lines'.format(count))
else:
    logging.info('succeed in classing linear RNA')

basic_column_list = ['transcript_id', 'gene_id'] + sample_list + ['level', 'category', 'kind']
lncrna_df = pd.DataFrame(lncrna_record_dict_list)
lncrna_df = lncrna_df.reindex(basic_column_list, axis=1)
lncrna_df.to_csv('lncRNA.txt', sep='\t', index=False)

whole_column_list = basic_column_list + ['gene_name', 'description', 'length', 'entrez', 'go', 'ko_id', 'ko_name', 'paths', 'cog', 'cog_description', 'nr', 'swissprot', 'pfam']
mrna_df = pd.DataFrame(mrna_record_dict_list)
mrna_df = mrna_df.reindex(whole_column_list, axis=1)
mrna_df.to_csv('mRNA.txt', sep='\t', index=False)

for idx, row in g_tpm_df.iterrows():
    record_dict = row.rename({'rna_type': 'category'}).to_dict()
    record_dict['gene_id'] = idx
    record_dict['kind'] = 'new' if record_dict.pop('is_new') else 'ref'
    record_dict['level'] = 'G'
    if record_dict['category'] == 'mRNA':
        pd_obj = g_annot_df.loc[(idx, 'yes')]
        annot_dict = get_annot_dict(pd_obj, level='G')
        record_dict.update(annot_dict)
        mrna_record_dict_list.append(record_dict)
    elif record_dict['category'] == 'lncRNA':
        lncrna_record_dict_list.append(record_dict)
    count += 1
    if not count % 10000:
        logging.debug('succeed in processing {} lines'.format(count))
else:
    logging.info('succeed in classing linear RNA')






















logging.info('succeed in exporting results')




