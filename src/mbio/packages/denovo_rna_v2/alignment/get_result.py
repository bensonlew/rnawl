# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging

import pandas as pd

from mbio.packages.denovo_rna_v2.functions import watcher

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


@watcher
def main(args):
    {'genome': parse_genome_hits, 'database': parse_database_hits, 'upload': parse_upload_hits}[args.subparsers](args)


@watcher
def parse_upload_hits(args):
    m8df = pd.read_table(args.m8_out, names=['query_id', 'subject_id', 'identity', 'length', 'query_start', 'query_end',
                                             'subject_start', 'subject_end', 'evalue', 'score'],
                         usecols=[0, 1, 2, 3, 6, 7, 8, 9, 10, 11])
    m8df.to_csv(args.output, sep='\t', index=False)


@watcher
def parse_database_hits(args):
    m8df = pd.read_table(args.m8_out,
                         names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                'sstart', 'send', 'evalue', 'bitscore'])
    bmdf = pd.read_table(args.bm_file, header=None)
    mrdf = pd.read_table(args.g2t2p, names=['gid', 'tid', 'pid']).fillna(method='pad')

    flag = False
    if args.db_type in ('transcript', 'cds'):
        lit_set = set(m8df['sseqid'].tolist())
        lar_set = set(mrdf['tid'].tolist())
        uni_set = lit_set & lar_set
        if float(len(uni_set)) / len(lit_set) > 0.5:
            flag = True
    elif args.db_type == 'peptide':
        lit_set = set(m8df['sseqid'].tolist())
        lar_set = set(mrdf['pid'].tolist())
        uni_set = lit_set & lar_set
        if float(len(uni_set)) / len(lit_set) > 0.5:
            flag = True

    def dropdot(cell):
        if pd.isna(cell) or args.bm_type == 'type3':
            ret = cell
        else:
            if flag:
                ret = cell
            else:
                if cell.rfind('.') == -1:
                    ret = cell
                else:
                    ret = cell[:cell.rfind('.')]
        return ret

    bmdf[1] = bmdf[1].apply(dropdot)
    bmdf = bmdf.set_index(1)
    mrdf['tid'] = mrdf['tid'].apply(dropdot)
    mrdf['pid'] = mrdf['pid'].apply(dropdot)
    mrdf = mrdf.set_index('pid')

    def row2tuple(i, row):
        query_id = row['qseqid']
        query_start = row['qstart']
        query_end = row['qend']
        subject_id = row['sseqid']
        if args.db_type == 'transcript':
            traid = row['sseqid']
            if flag:
                retid = traid
            else:
                retid = traid if traid.rfind('.') == -1 else traid[:traid.rfind('.')]
            idx = retid if args.bm_type in ('type1', 'type2') else retid
        elif args.db_type == 'cds':
            cdsid = row['sseqid']
            if flag:
                retid = cdsid
            else:
                retid = cdsid if cdsid.rfind('.') == -1 else cdsid[:cdsid.rfind('.')]
            idx = retid if args.bm_type in ('type1', 'type2') else cdsid
        elif args.db_type == 'peptide':
            pepid = row['sseqid']
            if flag:
                tmpid = pepid
            else:
                tmpid = pepid if pepid.rfind('.') == -1 else pepid[:pepid.rfind('.')]
            mridx = tmpid if args.bm_type in ('type1', 'type2') else pepid
            idx = mrdf.loc[mridx, 'tid']
        if isinstance(idx, pd.core.series.Series):
            idx = list(idx)[0]
        gene_id = bmdf.loc[idx][0]
        gene_name = '-' if args.bm_type == 'type3' else '-' if pd.isna(bmdf.loc[idx][2]) else bmdf.loc[idx][2]
        gene_description = bmdf.loc[idx][{'type1': 7, 'type2': 5, 'type3': 3}[args.bm_type]]
        gene_description = '-' if pd.isna(gene_description) else gene_description
        if row['sstart'] < row['send']:
            orientation = '+'
            subject_start = row['sstart']
            subject_end = row['send']
        else:
            orientation = '-'
            subject_start = row['send']
            subject_end = row['sstart']
        length = row['length']
        score = row['bitscore']
        evalue = row['evalue']
        identity = row['pident']
        if i and not i % 10000:
            logging.debug('succeed in parsing {} lines'.format(i))
        return (
            query_id, query_start, query_end, subject_id, gene_id, gene_name, gene_description, orientation,
            subject_start,
            subject_end, length, score, evalue, identity)

    data = [row2tuple(i, row) for i, row in m8df.iterrows()]
    redf = pd.DataFrame(data, columns=['query_id', 'query_start', 'query_end', 'subject_id', 'gene_id', 'gene_name',
                                       'gene_description', 'orientation', 'subject_start', 'subject_end', 'length',
                                       'score', 'evalue', 'identity'])
    redf.to_csv(args.output, sep='\t', index=False)


@watcher
def parse_genome_hits(args):
    m8df = pd.read_table(args.m8_out,
                         names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                'sstart', 'send', 'evalue', 'bitscore'])
    btdf = pd.read_table(args.bt_out, header=None)
    bmdf = pd.read_table(args.bm_file, header=None, index_col=1)

    def get_biomart_gene_info(txptid, bmdf, bmtype):
        gene_id = bmdf.loc[txptid][0]
        if bmtype == 'type1':
            gene_name = bmdf.loc[txptid][2]
            gene_description = bmdf.loc[txptid][7]
        elif bmtype == 'type2':
            gene_name = bmdf.loc[txptid][2]
            gene_description = bmdf.loc[txptid][5]
        elif bmtype == 'type3':
            gene_name = '-'
            gene_description = bmdf.loc[txptid][3]
        return gene_id, gene_name, gene_description

    def row2tuple(row):
        gene_id, gene_name, gene_description = get_biomart_gene_info(row[9], bmdf, args.bm_type)
        return (
            row[3], m8df.loc[row[4], 'qstart'], m8df.loc[row[4], 'qend'], '{}:{}-{}'.format(row[0], row[1] + 1, row[2]),
            gene_id, gene_name, gene_description, row[5], m8df.loc[row[4], 'length'], m8df.loc[row[4], 'bitscore'],
            m8df.loc[row[4], 'evalue'], m8df.loc[row[4], 'pident'])

    data = [row2tuple(row) for i, row in btdf.iterrows()]
    redf = pd.DataFrame(data,
                        columns=['query_id', 'query_start', 'query_end', 'genomic_location', 'overlapping_gene_id',
                                 'overlapping_gene_name', 'overlapping_gene_description', 'orientation', 'length',
                                 'score', 'evalue', 'identity'])
    redf = redf.drop_duplicates()
    redf['overlapping_gene_description'] = redf['overlapping_gene_description'].fillna('-')
    redf.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Obtain result table file')
    subparsers = parser.add_subparsers(dest='subparsers', help='sub-command help')

    parser_g = subparsers.add_parser('genome', help='parse outputs when algining with genome')
    parser_g.add_argument('--m8_out', action='store', required=True,
                          help='blast m8 tabular', metavar='<FILE>', dest='m8_out')
    parser_g.add_argument('--bt_out', action='store', required=True,
                          help='bedtools output', metavar='<FILE>', dest='bt_out')
    parser_g.add_argument('--bm_file', action='store', required=True,
                          help='biomart file', metavar='<FILE>', dest='bm_file')
    parser_g.add_argument('--bm_type', action='store', choices=['type1', 'type2', 'type3'], required=True,
                          help='biomart type', dest='bm_type')
    parser_g.add_argument('--output', action='store', required=True,
                          help='output file path', metavar='<FILE>', dest='output')

    parser_d = subparsers.add_parser('database', help='parse outputs when algining with database')
    parser_d.add_argument('--m8_out', action='store', required=True,
                          help='blast m8 tabular', metavar='<FILE>', dest='m8_out')
    parser_d.add_argument('--db_type', action='store', choices=['transcript', 'cds', 'peptide'], required=True,
                          help='database type', dest='db_type')
    parser_d.add_argument('--bm_file', action='store', required=True,
                          help='biomart file', metavar='<FILE>', dest='bm_file')
    parser_d.add_argument('--bm_type', action='store', choices=['type1', 'type2', 'type3'], required=True,
                          help='biomart type', dest='bm_type')
    parser_d.add_argument('--g2t2p', action='store', required=True,
                          help='relation map', metavar='<FILE>', dest='g2t2p')
    parser_d.add_argument('--output', action='store', required=True,
                          help='output file path', metavar='<FILE>', dest='output')

    parser_u = subparsers.add_parser('upload', help='parse outputs when existing upload sequences')
    parser_u.add_argument('--m8_out', action='store', required=True,
                          help='blast m8 tabular', metavar='<FILE>', dest='m8_out')
    parser_u.add_argument('--output', action='store', required=True,
                          help='output file path', metavar='<FILE>', dest='output')

    if len(sys.argv) > 1 and sys.argv[1] not in ['-h', '--help']:
        if sys.argv[1] in ['genome', 'database', 'upload']:
            args = parser.parse_args()
        else:
            logging.error('unrecognized command ({})'.format(sys.argv[1]))
            sys.exit(-2)
    else:
        parser.print_help()
        sys.exit(-1)

    main(args)
