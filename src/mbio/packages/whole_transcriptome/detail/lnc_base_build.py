# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import json
import logging
import pickle
import re
import sqlite3
import os
from biocluster.config import Config

import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

get_sequence_dict = lambda fasta: {record.id: str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}


def main(args):
    input_dict = json.load(open(args.json))
    lncrna_gtf = input_dict['lncrna_gtf']
    biomart_file = input_dict['biomart_file']
    biomart_type = input_dict['biomart_type']
    relation_file = input_dict['relation_file']


    transcript_fa = input_dict['transcript_fa']
    transcript_bed = input_dict['transcript_bed']
    organism_name = input_dict['organism_name']
    species_urls = input_dict['species_urls']
    source = input_dict['source']
    transcript_detail_pk = input_dict['transcript_detail_pk']


    transcript_seq_dict = get_sequence_dict(transcript_fa)
    logging.info('succeed in obtaining transcript sequence dict with {} items'.format(len(transcript_seq_dict)))

    t2g_dict = get_t2g_dict(lncrna_gtf)
    logging.info('succeed in obtaining transcript-gene id relation dict with {} items'.format(len(t2g_dict)))



    conn = sqlite3.connect(args.database)
    cursor = conn.cursor()

    for table_name, sequence_dict in (('transcript_seq', transcript_seq_dict), ):
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (seq_id TEXT, sequence TEXT)'.format(table_name))
        for seq_id, sequence in sequence_dict.items():
            cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, sequence))
        logging.info('succeed in creating {} TABLE in {}'.format(table_name, conn))
    cursor.execute('DROP TABLE IF EXISTS sub_info')
    cursor.execute(
        'CREATE TABLE sub_info (transcript_id TEXT, cds_id TEXT, pep_id TEXT, cds_seq TEXT, pep_seq TEXT, orf_type TEXT)')
    logging.info('succeed in creating sub_info TABLE in {}'.format(conn))
    conn.commit()
    conn.close()

    g2t_dict = get_g2t_dict(t2g_dict)
    logging.info('succeed in obtaining gene-transcripts id relation dict with {} items'.format(len(g2t_dict)))

    transcript_bed_dict = get_bed_dict(transcript_bed)
    logging.info('succeed in obtaining transcript bed dict with {} items'.format(len(transcript_bed_dict)))

    predict_dict  = get_detail_dict(input_dict['known_lncrna_detail'], input_dict['novel_lncrna_detail'], species_urls)

    transcript_detail_data = get_transcript_detail_data(t2g_dict, transcript_bed_dict, transcript_seq_dict,
                                                        organism_name, source, predict_dict)
    pickle.dump(transcript_detail_data, open(transcript_detail_pk, 'w'))
    logging.info('succeed in exporting {}'.format(transcript_detail_pk))
    seqdownload_lnc_seq = os.path.join(os.path.dirname(transcript_detail_pk), "seqdownload_lnc_seq")
    with open(seqdownload_lnc_seq, "wb") as f:
        pickle.dump(transcript_seq_dict, f)


def get_url_dict(species_urls):
    # 获取已知lncRNA数据库ID

    lnc_url = {
        "ENSEMBL" : '{}/Gene/Summary?g='.format(species_urls[:-11]) + '{}',
        "NCBI" : 'https://www.ncbi.nlm.nih.gov/search/all/?term={}',
        "NONCODE" : 'http://www.noncode.org/show_rna.php?id={}',
    }
    green_nc_ids = Config().SOFTWARE_DIR + '/database/lnc_rna/greenc/greenc_lncRNA_link.txt'
    with open(green_nc_ids, 'r') as f:
        green_dict = {line.strip().split()[0]: line.strip().split()[1] for line in f.readlines()}
    return lnc_url, green_dict


def get_detail_dict(known_detail, novel_detail, species_urls):
    lnc_url, green_dict = get_url_dict(species_urls)
    predict_dict = dict()
    if os.path.exists(known_detail):
        with open(known_detail, 'r') as f:
            header = f.readline()
            header_list = header.strip("\n").split("\t")
            type_col = header_list.index("type")
            database_col = header_list.index("database")
            for line in f:
                cols = line.strip("\n").split("\t")

                predict_dict[cols[0]] = {'kind': "ref", 'type': cols[type_col]}
                if cols[database_col].upper() in ["ENSEMBL", "NCBI", "NONCODE"]:
                    predict_dict[cols[0]]['url'] = lnc_url[cols[database_col].upper()].format(cols[0])
                if cols[database_col].upper() in ["NONCODE"]:
                    predict_dict[cols[0]]['url'] = lnc_url[cols[database_col].upper()].format(cols[0].split(".")[0])
                if cols[database_col].upper() == "GREENC":
                    predict_dict[cols[0]]['url'] = green_dict.get(cols[0], '')
    if os.path.exists(novel_detail):
        with open(novel_detail, 'r') as f:
            header = f.readline()
            header_list = header.strip("\n").split("\t")
            type_col = header_list.index("type")
            for line in f:
                cols = line.strip("\n").split("\t")
                predict_dict[cols[0]] = {'kind': "ref", 'type': cols[type_col]}
    return predict_dict

def get_transcript_detail_data(t2g_dict, transcript_bed_dict, transcript_seq_dict, organism_name, source, predict_dict):
    data = list()
    for transcript_id in t2g_dict:
        # url = get_url(transcript_id, t2g_dict, 'T', organism_name=organism_name, source=source)
        if transcript_id in transcript_bed_dict:
            document = {'transcript_id': transcript_id,
                        'chrom': transcript_bed_dict[transcript_id]['chrom'],
                        'start': transcript_bed_dict[transcript_id]['start'],
                        'end': transcript_bed_dict[transcript_id]['end'],
                        'strand': transcript_bed_dict[transcript_id]['strand'],
                        'length': len(transcript_seq_dict[transcript_id]),
                        'exon_num': transcript_bed_dict[transcript_id]['exon_num'],
                        'type': predict_dict[transcript_id]['type'],
                        'url': predict_dict[transcript_id].get('url', '')}

            data.append(document)
        else:
            print "transcript {} not in bed".format(transcript_id)
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
            dct[eles[3]] = {'chrom': eles[0], 'start': eles[1], 'end': eles[2], 'strand': eles[5], 'exon_num': eles[9]}
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




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Clean bedtools output FASTA file')
    parser.add_argument('--json', action='store', required=True,
                        help='setting JSON file', metavar='<FILE>', dest='json')
    parser.add_argument('--database', action='store', required=True,
                        help='SQLite database file', metavar='<FILE>', dest='database')

    args = parser.parse_args()

    main(args)
