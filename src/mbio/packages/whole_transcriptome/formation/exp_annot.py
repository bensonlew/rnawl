# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    t_exp_df = pd.read_table(args.t_exp, index_col=0)
    g_exp_df = pd.read_table(args.g_exp, index_col=0)
    need_cols = ['gene_id', 'gene_name', 'gene_description']
    lnc_des = pd.read_table(args.biomart, sep="\t").loc[:, need_cols].drop_duplicates()
    lnc_des.rename(columns={'gene_description': 'description'}, inplace=True)
    lnc_des = lnc_des.set_index('gene_id')
    query_obj = QueryBuilder(args.annot)
    relat_df = pd.read_table(args.relate, names=['transcript_id', 'gene_id', 'category', 'kind'], index_col=0)
    output_dir = args.output
    logging.info('succeed in loading requried data')
    export_transcript_result(t_exp_df, query_obj, relat_df, output_dir, lnc_des)
    logging.info('succeed in exporting transcript results in {}'.format(output_dir))
    export_gene_result(g_exp_df, query_obj, output_dir, lnc_des)
    logging.info('succeed in exporting gene results in {}'.format(output_dir))


def export_transcript_result(t_exp_df, query_obj, relat_df, output_dir, lnc_des):
    mrna_record_dict_list = list()
    lncrna_record_dict_list = list()
    for idx, row in t_exp_df.iterrows():
        record_dict = row.rename({'rna_type': 'category'}).to_dict()
        record_dict['transcript_id'] = idx
        record_dict['gene_id'] = relat_df.loc[idx]['gene_id']
        record_dict['kind'] = 'new' if record_dict.pop('is_new') else 'ref'
        record_dict['level'] = 'T'
        if record_dict['category'] == 'mRNA':
            annot_dict = query_obj.get_annot_dict(idx, level='T')
            record_dict.update(annot_dict)
            mrna_record_dict_list.append(record_dict)
        elif record_dict['category'] == 'lncRNA':
            if record_dict['gene_id'] in lnc_des.index.values:
                desc = lnc_des.loc[record_dict['gene_id']].to_dict()
            else:
                desc = {'gene_name': '-', 'description': '-'}
            record_dict.update(desc)
            lncrna_record_dict_list.append(record_dict)
    sample_list = list(t_exp_df.columns[:-2])
    basic_column_list = ['transcript_id', 'gene_id'] + sample_list + ['level', 'category', 'kind']
    lncrna_df = pd.DataFrame(lncrna_record_dict_list)
    lncrna_df = lncrna_df.reindex(basic_column_list + ['gene_name', 'description'], axis=1)
    lncrna_df.to_csv(os.path.join(output_dir, 'T.lncRNA.txt'), sep='\t', index=False)
    whole_column_list = basic_column_list + ['gene_name', 'description', 'length', 'entrez', 'go', 'ko_id', 'ko_name',
                                             'paths', 'cog', 'cog_description', 'nr', 'swissprot', 'pfam',
                                             'pathways_class1', 'pathways_class2']
    mrna_df = pd.DataFrame(mrna_record_dict_list)
    mrna_df = mrna_df.reindex(whole_column_list, axis=1)
    mrna_df.to_csv(os.path.join(output_dir, 'T.mRNA.txt'), sep='\t', index=False)


def export_gene_result(g_exp_df, query_obj, output_dir, lnc_des):
    mrna_gene_record_dict_list = list()
    lncrna_gene_record_dict_list = list()
    for idx, row in g_exp_df.iterrows():
        record_dict = row.rename({'rna_type': 'category'}).to_dict()
        record_dict['gene_id'] = idx
        record_dict['kind'] = 'new' if record_dict.pop('is_new') else 'ref'
        record_dict['level'] = 'G'
        if record_dict['category'] == 'mRNA':
            annot_dict = query_obj.get_annot_dict(idx, level='G')
            record_dict.update(annot_dict)
            mrna_gene_record_dict_list.append(record_dict)
        elif record_dict['category'] == 'lncRNA':
            if idx in lnc_des.index.values:
                desc = lnc_des.loc[idx].to_dict()
            else:
                desc = {'gene_name': '-', 'description': '-'}
            record_dict.update(desc)
            lncrna_gene_record_dict_list.append(record_dict)
    sample_list = list(g_exp_df.columns[:-2])
    basic_column_list = ['gene_id'] + sample_list + ['level', 'category', 'kind']
    lncrna_df = pd.DataFrame(lncrna_gene_record_dict_list)
    lncrna_df = lncrna_df.reindex(basic_column_list + ['gene_name', 'description'], axis=1)
    lncrna_df.to_csv(os.path.join(output_dir, 'G.lncRNA.txt'), sep='\t', index=False)
    whole_column_list = basic_column_list + ['gene_name', 'description', 'transcript_id', 'length', 'entrez', 'go',
                                             'ko_id', 'ko_name', 'paths', 'cog', 'cog_description', 'nr', 'swissprot',
                                             'pfam', 'pathways_class1', 'pathways_class2']
    mrna_df = pd.DataFrame(mrna_gene_record_dict_list)
    mrna_df = mrna_df.reindex(whole_column_list, axis=1)
    mrna_df.to_csv(os.path.join(output_dir, 'G.mRNA.txt'), sep='\t', index=False)


class QueryBuilder:
    def __init__(self, annot_file):
        self._df = pd.read_table(annot_file)
        self.t_df = self._df.set_index('transcript_id')
        self.g_df = self._df.set_index(['gene_id', 'is_gene'])

    def get_annot_dict(self, idx, level):
        if level == 'T':
            pd_obj = self.t_df.loc[idx]
        elif level == 'G':
            pd_obj = self.g_df.loc[idx].loc['yes']
        assert isinstance(pd_obj, pd.core.series.Series)
        annot_dict = pd_obj.rename(str.lower).to_dict()
        return annot_dict


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate files for dumping data into mongo')
    parser.add_argument('-t', action='store', required=True,
                        help='transcript expression matrix file', metavar='<FILE>', dest='t_exp')
    parser.add_argument('-g', action='store', required=True,
                        help='gene expression matrix file', metavar='<FILE>', dest='g_exp')
    parser.add_argument('-a', action='store', required=True,
                        help='annotation matrix file', metavar='<FILE>', dest='annot')
    parser.add_argument('-r', action='store', required=True,
                        help='relationship map file', metavar='<FILE>', dest='relate')
    parser.add_argument('-b', action='store', required=True,
                        help='desc file', metavar='<FILE>', dest='biomart')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
