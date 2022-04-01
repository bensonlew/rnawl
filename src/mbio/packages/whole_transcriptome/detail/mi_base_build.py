# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import pickle
import re
import sqlite3
import os

import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

get_sequence_dict = lambda fasta: {record.id: str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}


def main(args):
    input_dict = json.load(open(args.json))

    transcript_detail_pk = input_dict['transcript_detail_pk']
    detail_dict, mirna_seq = get_detail_dict(input_dict['known_mirna_detail'], input_dict['novel_mirna_detail'])
    transcript_detail_data = get_transcript_detail_data(detail_dict)
    pickle.dump(transcript_detail_data, open(transcript_detail_pk, 'w'))
    logging.info('succeed in exporting {}'.format(transcript_detail_pk))
    seqdownload_mirna_seq = os.path.join(os.path.dirname(transcript_detail_pk), "seqdownload_mirna_seq")
    with open(seqdownload_mirna_seq, "wb") as f:
        pickle.dump(mirna_seq, f)


def get_detail_dict(known_detail, novel_detail):
    predict_dict = dict()
    mirna_seq_dict = dict()
    if os.path.exists(known_detail):
        with open(known_detail, 'r') as f:
            header = f.readline()

            for line in f:
                cols = line.strip("\n").split("\t")
                try:
                    chrom = cols[6].split("(")[0]
                    start = cols[6].split(":")[1].split("-")[0]
                    end = cols[6].split(":")[1].split("-")[1]
                    strand = cols[6].split("(")[1][0]
                except:
                    chrom = "-"
                    start = "-"
                    end = "-"
                    strand = "-"
                if cols[0] in predict_dict:
                    predict_dict[cols[0]]['pre_name'] = predict_dict[cols[0]]['pre_name'] + "," + "{}({})".format(cols[3], cols[5])
                    predict_dict[cols[0]]['pre_len'] = predict_dict[cols[0]]['pre_len'] + "," + "{}".format(cols[5])
                else:
                    predict_dict[cols[0]] = {'kind': "ref", 'transcript_id': cols[0],
                                             'mirna_seq': cols[1],
                                             'length': cols[2],
                                             'pre_name': "{}({})".format(cols[3], cols[5]),
                                             'pre_len': cols[5],
                                             'chrom': chrom,
                                             'start': start,
                                             'end': end,
                                             'strand': strand}
                    predict_dict[cols[0]]['url'] = 'http://www.mirbase.org/cgi-bin/query.pl?terms={}'.format(cols[0])
                mirna_seq_dict[cols[0]] = cols[1]
    if os.path.exists(novel_detail):
        with open(novel_detail, 'r') as f:
            header = f.readline()

            for line in f:
                cols = line.strip("\n").split("\t")
                chrom = cols[6].split("(")[0]
                start = cols[6].split(":")[1].split("-")[0]
                end = cols[6].split(":")[1].split("-")[1]
                strand = cols[6].split("(")[1][0]
                if cols[0] in predict_dict:
                    predict_dict[cols[0]]['pre_name'] = predict_dict[cols[0]]['pre_name'] + "," + "{}({})".format(cols[3], cols[5])
                    predict_dict[cols[0]]['pre_len'] = predict_dict[cols[0]]['pre_len'] + "," + "{}".format(cols[5])

                else:
                    predict_dict[cols[0]] = {'kind': "new", 'transcript_id': cols[0],
                                             'mirna_seq': cols[1],
                                             'length': cols[2],
                                             'pre_name': "{}({})".format(cols[3], cols[5]),
                                             'pre_len': cols[5],
                                             'chrom': chrom,
                                             'start': start,
                                             'end': end,
                                             'strand': strand}
                    predict_dict[cols[0]]['url'] = 'http://www.mirbase.org/cgi-bin/query.pl?terms={}'.format(cols[0])
                mirna_seq_dict[cols[0]] = cols[1]
    return predict_dict, mirna_seq_dict


def get_transcript_detail_data(detail_dict):
    data = list()
    for transcript_id in detail_dict:
        document = detail_dict[transcript_id]
        data.append(document)
    return data


def get_gene_detail_data(g2t_dict, gene_bed_dict, gene_seq_dict, transcript_bed_dict, transcript_seq_dict,
                         organism_name, source):
    data = list()
    for gene_id, transcript_ids in g2t_dict.items():
        transcripts = [get_transcript_info(transcript_id, transcript_bed_dict, transcript_seq_dict) for transcript_id in
                       transcript_ids]
        url = get_url(gene_id, organism_name=organism_name, source=source)
        document = {'gene_id': gene_id, 'chrom': gene_bed_dict[gene_id]['chrom'],
                    'start': gene_bed_dict[gene_id]['start'], 'end': gene_bed_dict[gene_id]['end'],
                    'strand': gene_bed_dict[gene_id]['strand'], 'length': len(gene_seq_dict[gene_id]),
                    'transcripts': transcripts, 'transcript_num': len(transcript_ids), 'url': url}
        data.append(document)
    return data


def get_transcript_info(transcript_id, transcript_bed_dict, transcript_seq_dict):
    return {'transcript_id': transcript_id,
            'start': transcript_bed_dict[transcript_id]['start'],
            'end': transcript_bed_dict[transcript_id]['end'],
            'strand': transcript_bed_dict[transcript_id]['strand'],
            'length': len(transcript_seq_dict[transcript_id])}


def get_url(seq_id, t2g_dict=dict(), level='G', organism_name=None, source=None):
    if source == 'ensembl':
        if level == 'G':
            url = 'asia.ensembl.org/{}/Gene/Summary?g={}'.format(organism_name, seq_id)
        elif level == 'T':
            url = 'asia.ensembl.org/{}/Transcript/Summary?g={};t={}'.format(organism_name, t2g_dict[seq_id], seq_id)
    elif source == 'ncbi':
        url = 'www.ncbi.nlm.nih.gov/search/all/?term={}'.format(seq_id)
    else:
        url = None
    return url


def get_bed_dict(bed):
    dct = dict()
    for line in open(bed):
        eles = line.strip().split('\t')
        if len(eles) > 5:
            dct[eles[3]] = {'chrom': eles[0], 'start': eles[1], 'end': eles[2], 'strand': eles[5]}
    return dct


def get_g2t_dict(t2g_dict):
    dct = dict()
    for transcript_id, gene_id in t2g_dict.items():
        if gene_id in dct:
            dct[gene_id].add(transcript_id)
        else:
            dct[gene_id] = {transcript_id}
    return dct


def get_sub_info_dict(t2g_dict, cds_info_dict, pep_info_dict):
    dct = dict()
    for transcript_id in t2g_dict:
        dct[transcript_id] = {
            'transcript_id': transcript_id,
            'cds_id': cds_info_dict.get(transcript_id, {'id': 'None'})['id'],
            'pep_id': pep_info_dict.get(transcript_id, {'id': 'None'})['id'],
            'cds_seq': cds_info_dict.get(transcript_id, {'sequence': 'None'})['sequence'],
            'pep_seq': pep_info_dict.get(transcript_id, {'sequence': 'None'})['sequence'],
            'orf_type': cds_info_dict.get(transcript_id, {'orf_type': 'complete'})['orf_type']}
    return dct


def get_t2g_dict(gtf):
    dct = dict()
    pg = re.compile(r'gene_id "(\S+)";')
    pt = re.compile(r'transcript_id "(\S+)";')
    for line in open(gtf):
        if line.strip() and line[0] != '#':
            eles = line.strip().split('\t')
            if len(eles) >= 8 and 'gene_id' in eles[8] and 'transcript_id' in eles[8]:
                mg = re.search(pg, eles[8])
                mt = re.search(pt, eles[8])
                if mg and mt:
                    g = mg.group(1)
                    t = mt.group(1)
                    dct[t] = g
    return dct


def get_p2t_dict(biomart_file, biomart_type, relation_file):
    df = pd.read_table(biomart_file, usecols=[1, {'type1': 6, 'type2': 4, 'type3': 2}[biomart_type]])
    df.columns = ('transcript_id', 'pep_id')
    df = df.dropna()
    dct = df.set_index('pep_id').to_dict()['transcript_id']
    for line in open(relation_file):
        eles = line.strip().split('\t')
        if len(eles) == 5:
            transcript_id = eles[0]
            pep_id = eles[4]
            dct[pep_id] = transcript_id
    return dct


def get_cds_info_dict(ref_cds_fa, new_cds_fa, t2g_dict):
    dct = dict()
    for record in SeqIO.parse(ref_cds_fa, 'fasta'):
        cds_id = record.id
        transcript_id = cds_id[:cds_id.rfind('.')] if '.' in cds_id else cds_id
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': cds_id, 'sequence': str(record.seq), 'orf_type': 'complete'}
    p = re.compile(r'(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*')
    for record in SeqIO.parse(new_cds_fa, 'fasta'):
        gene_id, transcript_id, orf_type, position = re.match(p, record.description).groups()
        cds_id = record.id
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': cds_id, 'sequence': str(record.seq), 'orf_type': orf_type}
    return dct


def get_pep_info_dict(ref_pep_fa, new_pep_fa, t2g_dict, p2t_dict):
    dct = dict()
    for record in SeqIO.parse(ref_pep_fa, 'fasta'):
        pep_id = record.id
        ret = p2t_dict.get(pep_id)
        if ret:
            transcript_id = ret
        else:
            ret = p2t_dict.get(pep_id[:pep_id.rfind('.')])
            if ret:
                transcript_id = ret
            else:
                continue
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': pep_id, 'sequence': str(record.seq)}
    p = re.compile(r'(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*')
    for record in SeqIO.parse(new_pep_fa, 'fasta'):
        gene_id, transcript_id, orf_type, position = re.match(p, record.description).groups()
        pep_id = record.id
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': pep_id, 'sequence': str(record.seq)}
    return dct


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Clean bedtools output FASTA file')
    parser.add_argument('--json', action='store', required=True,
                        help='setting JSON file', metavar='<FILE>', dest='json')
    parser.add_argument('--database', action='store', required=True,
                        help='SQLite database file', metavar='<FILE>', dest='database')

    args = parser.parse_args()

    main(args)
