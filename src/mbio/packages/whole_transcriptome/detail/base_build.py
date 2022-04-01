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
    mrna_gtf = input_dict['mrna_gtf']
    biomart_file = input_dict['biomart_file']
    biomart_type = input_dict['biomart_type']
    entrez = input_dict['entrez']
    species_urls = input_dict['species_urls']
    relation_file = input_dict['relation_file']
    ref_cds_fa = input_dict['ref_cds_fa']
    new_cds_fa = input_dict['new_cds_fa']
    ref_pep_fa = input_dict['ref_pep_fa']
    new_pep_fa = input_dict['new_pep_fa']
    gene_fa = input_dict['gene_fa']
    transcript_fa = input_dict['transcript_fa']
    gene_bed = input_dict['gene_bed']
    transcript_bed = input_dict['transcript_bed']
    organism_name = input_dict['organism_name']
    source = input_dict['source']
    gene_detail_pk = input_dict['gene_detail_pk']
    transcript_detail_pk = input_dict['transcript_detail_pk']
    trans_type = input_dict['trans_type']
    trans_type_dict = get_trans_type(trans_type)
    gene_enterz_dict = get_gene_enterz(entrez)
    gene_seq_dict = get_sequence_dict(gene_fa)
    logging.info('succeed in obtaining gene sequence dict with {} items'.format(len(gene_seq_dict)))
    transcript_seq_dict = get_sequence_dict(transcript_fa)
    logging.info('succeed in obtaining transcript sequence dict with {} items'.format(len(transcript_seq_dict)))
    t2g_dict = get_t2g_dict(mrna_gtf)
    logging.info('succeed in obtaining transcript-gene id relation dict with {} items'.format(len(t2g_dict)))
    p2t_dict = get_p2t_dict(biomart_file, biomart_type, relation_file)
    gene2name, gene2des, trans2des = get_des(biomart_file, biomart_type)
    logging.info('succeed in obtaining peptide-transcript id relation dict with {} items'.format(len(p2t_dict)))
    cds_info_dict, cds_seq = get_cds_info_dict(ref_cds_fa, new_cds_fa, t2g_dict)
    logging.info('succeed in obtaining cds information dict with {} items'.format(len(cds_info_dict)))
    pep_info_dict, pep_seq = get_pep_info_dict(ref_pep_fa, new_pep_fa, t2g_dict, p2t_dict)
    logging.info('succeed in obtaining peptide information dict with {} items'.format(len(pep_info_dict)))
    sub_info_dict = get_sub_info_dict(t2g_dict, cds_info_dict, pep_info_dict)
    logging.info('succeed in obtaining {} documents for inserting to database'.format(len(sub_info_dict)))

    conn = sqlite3.connect(args.database)
    cursor = conn.cursor()
    for table_name, sequence_dict in (('gene_seq', gene_seq_dict), ('transcript_seq', transcript_seq_dict)):
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (seq_id TEXT, sequence TEXT)'.format(table_name))
        for seq_id, sequence in sequence_dict.items():
            cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, sequence))
        logging.info('succeed in creating {} TABLE in {}'.format(table_name, conn))
    cursor.execute('DROP TABLE IF EXISTS sub_info')
    cursor.execute(
        'CREATE TABLE sub_info (transcript_id TEXT, cds_id TEXT, pep_id TEXT, cds_seq TEXT, pep_seq TEXT, orf_type TEXT)')
    for sub_info in sub_info_dict.values():
        cursor.execute(
            "INSERT INTO sub_info VALUES ('{transcript_id}', '{cds_id}', '{pep_id}', '{cds_seq}', '{pep_seq}', '{orf_type}')".format(
                **sub_info))
    logging.info('succeed in creating sub_info TABLE in {}'.format(conn))
    conn.commit()
    conn.close()

    g2t_dict = get_g2t_dict(t2g_dict)
    logging.info('succeed in obtaining gene-transcripts id relation dict with {} items'.format(len(g2t_dict)))
    gene_bed_dict = get_bed_dict(gene_bed)
    logging.info('succeed in obtaining gene bed dict with {} items'.format(len(gene_bed_dict)))
    transcript_bed_dict = get_bed_dict(transcript_bed)
    logging.info('succeed in obtaining transcript bed dict with {} items'.format(len(transcript_bed_dict)))
    gene_detail_data = get_gene_detail_data(g2t_dict, gene_bed_dict, gene_seq_dict, transcript_bed_dict,
                                            transcript_seq_dict, organism_name, source, gene_enterz_dict,
                                            species_urls, gene2name, gene2des, trans_type_dict)
    pickle.dump(gene_detail_data, open(gene_detail_pk, 'w'))
    logging.info('succeed in exporting {}'.format(gene_detail_pk))
    transcript_detail_data = get_transcript_detail_data(t2g_dict, transcript_bed_dict, transcript_seq_dict,
                                                        organism_name, source, species_urls, trans2des)
    pickle.dump(transcript_detail_data, open(transcript_detail_pk, 'w'))
    logging.info('succeed in exporting {}'.format(transcript_detail_pk))
    output_dir = os.path.dirname(transcript_detail_pk)
    gene_detail = dict()
    tran_detail = dict()
    for gene_id in g2t_dict:
        trans_id = g2t_dict[gene_id]
        gene_detail[gene_id] = {'tran_id': trans_id, 'cds_id': [cds_info_dict[i]['id'] if i in cds_info_dict else '-' for i in trans_id ], 'pep_id': [pep_info_dict[j]['id'] if j in pep_info_dict else '-' for j in trans_id ]}
        for tran_id in trans_id:
            try:
                if tran_id not in tran_detail:
                    tran_detail[tran_id]={
                        "cds_id":[],
                        "pep_id":[]
                    }
                if tran_id in cds_info_dict:
                    tran_detail[tran_id]["cds_id"].append(cds_info_dict[tran_id]['id'])
                    # tran_detail[tran_id] = {'cds_id': [cds_info_dict[tran_id]['id']]}
                else:
                    tran_detail[tran_id]["cds_id"].append('-')
                if tran_id in pep_info_dict:
                    tran_detail[tran_id]['pep_id'].append(pep_info_dict[tran_id]['id'])
                else:
                    tran_detail[tran_id]['pep_id'].append('-')
                #     # tran_detail[tran_id] = {'cds_id': [cds_info_dict[tran_id]['id']], 'pep_id': [pep_info_dict[tran_id]['id']]}
                # else:
                #     if tran_id in cds_info_dict:
                #         tran_detail[tran_id]['cds_id'].append(cds_info_dict[tran_id]['id'])
                #     else:
                #         tran_detail[tran_id]['cds_id'].append('-')
                #     if tran_id in pep_info_dict:
                #         tran_detail[tran_id]['pep_id'].append(pep_info_dict[tran_id]['id'])
                #     else:
                #         tran_detail[tran_id]['pep_id'].append('-')
            except:
                pass

    gene_id_list = list()
    tran_num = list()
    cds_num = list()
    pep_num = list()
    for i in gene_detail:
        gene_id_list.append(i)
        tran_id_list = gene_detail[i]['tran_id']
        tran_num.append(len(tran_id_list))
        cds_id_list = gene_detail[i]['cds_id']
        pep_id_list = gene_detail[i]['pep_id']
        cds_num.append(len([i for i in cds_id_list if i != '-']))
        pep_num.append(len([i for i in pep_id_list if i != '-']))
    gene_stat_df = pd.DataFrame({'gene_id': gene_id_list, 'tran_num': tran_num, 'cds_num': cds_num, 'pep_num': pep_num})
    seqdownload_gene_stat = os.path.join(output_dir, "seqdownload_gene_stat")
    gene_stat_df.to_csv(seqdownload_gene_stat, sep="\t", index=False, header=True)
    tran_id_list = list()
    cds_num_tran = list()
    pep_num_tran = list()
    for a in tran_detail:
        tran_id_list.append(a)
        if tran_detail[a] == '':
            continue
        if 'cds_id' in tran_detail[a]:
            cds_num_tran.append(len([i for i in tran_detail[a]['cds_id'] if i != '-']))
        else:
            cds_num_tran.append(0)
        if 'pep_id' in tran_detail[a]:
            pep_num_tran.append(len([i for i in tran_detail[a]['pep_id'] if i != '-']))
        else:
            pep_num_tran.append(0)
    tran_stat_df = pd.DataFrame({'tran_id': tran_id_list, 'cds_num': cds_num_tran, 'pep_num': pep_num_tran})
    seqdownload_tran_stat = os.path.join(output_dir, "seqdownload_tran_stat")
    tran_stat_df.to_csv(seqdownload_tran_stat, sep="\t", index=False, header=True)
    seqdownload_gene_seq = os.path.join(output_dir, "seqdownload_gene_seq")
    seqdownload_txpt_seq = os.path.join(output_dir, "seqdownload_txpt_seq")
    seqdownload_cds_seq = os.path.join(output_dir, "seqdownload_cds_seq")
    seqdownload_pep_seq = os.path.join(output_dir, "seqdownload_pep_seq")
    seqdownload_gene_detail = os.path.join(output_dir, "seqdownload_gene_detail")
    seqdownload_tran_detail = os.path.join(output_dir, "seqdownload_tran_detail")
    with open(seqdownload_gene_seq, "wb") as f:
        pickle.dump(gene_seq_dict, f)
    with open(seqdownload_txpt_seq, "wb") as f:
        pickle.dump(transcript_seq_dict, f)
    with open(seqdownload_cds_seq, "wb") as f:
        pickle.dump(cds_seq, f)
    with open(seqdownload_pep_seq, "wb") as f:
        pickle.dump(pep_seq, f)
    with open(seqdownload_gene_detail, "w") as f:
        pickle.dump(gene_detail, f)
    with open(seqdownload_tran_detail, "w") as f:
        pickle.dump(tran_detail, f)

def get_gene_enterz(enterz):
    gene2enterz = dict()
    with open(enterz, 'r') as f:
        for line in f:
            cols = line.strip("\n").split("\t")
            if len(cols) >= 3:
                if cols[2] != "":
                    gene2enterz[cols[0]] = cols[2]
    return gene2enterz


def get_des(des, des_type):
    """
    获取蛋白名或描述信息
    """
    gene2name = dict()
    gene2des = dict()
    trans2des = dict()

    with open(des, "rb") as f:
        # f.readline()
        for line in f.readlines():
            line = line.strip().split('\t')
            gene_id = line[0]
            tran_id = line[1]
            if des_type == "type3":
                symbol = ""
                des = line[3]
            elif des_type == "type2":
                symbol = line[2]
                des = line[5]
            elif des_type == "type1":
                symbol = line[2]
                des = line[7]
            else:
                pass

            gene2name[gene_id] = symbol
            gene2des[gene_id] = des
            trans2des[tran_id] = des
    return gene2name, gene2des, trans2des

def get_transcript_detail_data(t2g_dict, transcript_bed_dict, transcript_seq_dict, organism_name, source, species_urls, trans2des):
    data = list()
    for transcript_id in t2g_dict:
        url = get_url(transcript_id, t2g_dict, 'T', organism_name=organism_name, source=source, species_urls=species_urls)
        document = {'transcript_id': transcript_id,
                    'chrom': transcript_bed_dict[transcript_id]['chrom'],
                    'start': transcript_bed_dict[transcript_id]['start'],
                    'end': transcript_bed_dict[transcript_id]['end'],
                    'strand': transcript_bed_dict[transcript_id]['strand'],
                    'length': len(transcript_seq_dict[transcript_id]),
                    'url': url,
                    'description': trans2des.get(transcript_id, "")
        }
        data.append(document)
    return data


def get_gene_detail_data(g2t_dict, gene_bed_dict, gene_seq_dict, transcript_bed_dict, transcript_seq_dict,
                         organism_name, source, gene_enterz_dict, species_urls, gene2name, gene2des, trans_type_dict):
    data = list()
    for gene_id, transcript_ids in g2t_dict.items():
        transcripts = [get_transcript_info(transcript_id, transcript_bed_dict, transcript_seq_dict, trans_type_dict) for transcript_id in
                       transcript_ids]
        url = get_url(gene_id, organism_name=organism_name, source=source, species_urls=species_urls)
        entrez = gene_enterz_dict.get(gene_id, "")
        document = {'gene_id': gene_id, 'chrom': gene_bed_dict[gene_id]['chrom'],
                    'start': gene_bed_dict[gene_id]['start'], 'end': gene_bed_dict[gene_id]['end'],
                    'strand': gene_bed_dict[gene_id]['strand'], 'length': len(gene_seq_dict[gene_id]),
                    'gene_name': gene2name.get(gene_id, ""),
                    'description': gene2des.get(gene_id, ""),
                    'transcripts': transcripts, 'transcript_num': len(transcript_ids), 'url': url,  "entrez": entrez}
        data.append(document)
    return data


def get_transcript_info(transcript_id, transcript_bed_dict, transcript_seq_dict, trans_type_dict):
    return {'transcript_id': transcript_id,
            'category': trans_type_dict[transcript_id],
            'start': transcript_bed_dict[transcript_id]['start'],
            'end': transcript_bed_dict[transcript_id]['end'],
            'strand': transcript_bed_dict[transcript_id]['strand'],
            'length': len(transcript_seq_dict[transcript_id])}


def get_url(seq_id, t2g_dict=dict(), level='G', organism_name=None, source=None, species_urls=None):
    if source == 'ensembl':
        if "/Info" in species_urls:
            species_urls = species_urls.split("/Info")[0]
        if level == 'G':
            url = '{}/Gene/Summary?g={}'.format(species_urls, seq_id)
        elif level == 'T':
            url = '{}/Transcript/Summary?g={};t={}'.format(species_urls, t2g_dict[seq_id], seq_id)
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

def get_trans_type(trans_type):
    dct = dict()
    with open(trans_type, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            if items[0] not in dct:
                dct[items[0]] = items[2]
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
    dct_cds = dict()
    for record in SeqIO.parse(ref_cds_fa, 'fasta'):
        cds_id = record.id
        if cds_id in t2g_dict:
            dct[cds_id] = {'id': cds_id, 'sequence': str(record.seq), 'orf_type': 'complete'}
            dct_cds[cds_id] = str(record.seq)
        else:
            transcript_id = cds_id[:cds_id.rfind('.')] if '.' in cds_id else cds_id
            if transcript_id in t2g_dict:
                dct[transcript_id] = {'id': cds_id, 'sequence': str(record.seq), 'orf_type': 'complete'}
                dct_cds[transcript_id] = str(record.seq)
    p = re.compile(r'(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*')
    for record in SeqIO.parse(new_cds_fa, 'fasta'):
        gene_id, transcript_id, orf_type, position = re.match(p, record.description).groups()
        cds_id = record.id
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': cds_id, 'sequence': str(record.seq), 'orf_type': orf_type}
            dct_cds[cds_id] = str(record.seq)
    return dct, dct_cds


def get_pep_info_dict(ref_pep_fa, new_pep_fa, t2g_dict, p2t_dict):
    dct = dict()
    dct_pep = dict()
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
            dct_pep[pep_id] = str(record.seq)
    p = re.compile(r'(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*')
    for record in SeqIO.parse(new_pep_fa, 'fasta'):
        gene_id, transcript_id, orf_type, position = re.match(p, record.description).groups()
        pep_id = record.id
        if transcript_id in t2g_dict:
            dct[transcript_id] = {'id': pep_id, 'sequence': str(record.seq)}
            dct_pep[pep_id] = str(record.seq)
    return dct, dct_pep


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Clean bedtools output FASTA file')
    parser.add_argument('--json', action='store', required=True,
                        help='setting JSON file', metavar='<FILE>', dest='json')
    parser.add_argument('--database', action='store', required=True,
                        help='SQLite database file', metavar='<FILE>', dest='database')

    args = parser.parse_args()

    main(args)
