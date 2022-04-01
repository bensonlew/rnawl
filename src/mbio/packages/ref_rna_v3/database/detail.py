# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import re
import sqlite3
from collections import defaultdict
import pickle
from Bio import SeqIO
import pandas as pd
import os
import shutil
from mbio.packages.ref_rna_v3.functions import pkgsfuncdeco


@pkgsfuncdeco
def main(args):
    gene_seq = fasta_to_dict(args.gf)
    txpt_seq = fasta_to_dict(args.tf)
    biomart_info, p2t_dict = parse_biomart(args.bf, args.bt)
    t2g_dict = dict([x.strip().split()[:2] for x in open(args.t2g) if x.strip()])
    g2t_dict = getg2t(args.t2g)
    known_cds,t2c_dict = get_cds_seq(args.kcf)
    known_pep,t2p_dict,pep_dict = get_pep_seq(args.kpf, p2t_dict)
    new_cds_pep,new_t2cds_dict,new_t2pep_dict,new_cds_dict,new_pep_dict = get_predict_cds_pep(args.ncf, args.npf) if args.ncf and args.npf else dict(),dict(),dict(),dict(),dict()
    build_database(args.db, gene_seq, txpt_seq, t2g_dict, known_cds, known_pep, new_cds_pep)
    merge_results(out=args.detail,gene_seq = gene_seq ,txpt_seq = txpt_seq,g2t_dict=g2t_dict,t2c_dict=t2c_dict,t2p_dict =t2p_dict,
            known_cds_dict = known_cds,known_pep_dict = pep_dict,new_pep_dict=new_pep_dict,
            new_cds_dict=new_cds_dict,new_t2cds_dict= new_t2cds_dict,new_t2pep_dict=new_t2pep_dict)


    with open("gene_seq","wb") as f:
        pickle.dump(gene_seq, f)
        # f.write(str(gene_seq))
    with open("txpt_seq","wb") as f:
        pickle.dump(txpt_seq, f)
        # f.write(str(txpt_seq))
    with open("biomart_info","wb") as f:
        pickle.dump(biomart_info, f)
        # f.write(str(biomart_info))
    with open("p2t_dict","wb") as f:
        # f.write(str(p2t_dict))
        pickle.dump(p2t_dict, f)
    with open("t2g_dict", "wb") as f:
        pickle.dump(t2g_dict, f)
        # f.write(str(t2g_dict))
    with open("known_cds", "w") as f:
        pickle.dump(known_cds, f)
        # f.write(str(known_cds))
    with open("known_pep", "w") as f:
        pickle.dump(known_pep, f)
        # f.write(str(known_pep))
    with open("new_cds_pep", "w") as f:
        pickle.dump(new_cds_pep, f)
        # f.write(str(new_cds_pep))

@pkgsfuncdeco
def merge_results(out=None,gene_seq = None,txpt_seq = None,g2t_dict=None,t2c_dict=None,t2p_dict =None,
            known_cds_dict = None,known_pep_dict = None,new_pep_dict=None,
            new_cds_dict=None,new_t2cds_dict= None,new_t2pep_dict=None):
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    gene_detail_list = {}
    gene_stat_list = []
    trans_detail_list = {}
    trans_stat_list = []
    for gene_id in g2t_dict:
        gene_detail = dict()
        gene_stat = dict()
        gene_detail["gene_id"] = gene_id
        gene_stat["gene_id"] = gene_id
        gene_detail["trans_id"] = g2t_dict[gene_id]
        gene_stat["trans_num"] = len(g2t_dict[gene_id])
        gene_detail.setdefault("cds_id",[])
        gene_stat.setdefault("cds_num",0)
        gene_detail.setdefault("pep_id", [])
        gene_stat.setdefault("pep_num", 0)
        for trans_id in g2t_dict[gene_id]:
            trans_detail = dict()
            trans_stat = dict()
            trans_detail["trans_id"] = trans_id
            trans_stat["trans_id"] = trans_id
            if trans_id in t2c_dict:
                trans_detail["cds_id"] = t2c_dict[trans_id]
                trans_stat["cds_num"] = len(t2c_dict[trans_id])
                gene_detail["cds_id"].append(t2c_dict[trans_id])
                gene_stat["cds_num"] += len(t2c_dict[trans_id])
            elif trans_id in new_t2cds_dict:
                trans_detail["cds_id"] = new_t2cds_dict[trans_id]
                trans_stat["cds_num"] = len(new_t2cds_dict[trans_id])
                gene_detail["cds_id"].append(new_t2cds_dict[trans_id])
                gene_stat["cds_num"] += len(new_t2cds_dict[trans_id])
            else:
                trans_detail["cds_id"] = ""
                trans_stat["cds_num"] = 0
            if trans_id in t2p_dict:
                trans_detail["pep_id"] = t2p_dict[trans_id]
                trans_stat["pep_num"] = len(t2p_dict[trans_id])
                gene_detail["pep_id"].append(t2p_dict[trans_id])
                gene_stat["pep_num"] += len(t2p_dict[trans_id])
            elif trans_id in new_t2pep_dict:
                trans_detail["pep_id"] = new_t2pep_dict[trans_id]
                trans_stat["pep_num"] = len(new_t2pep_dict[trans_id])
                gene_detail["pep_id"].append(new_t2pep_dict[trans_id])
                gene_stat["pep_num"] += len(new_t2pep_dict[trans_id])
            else:
                trans_detail["pep_id"] = ""
                trans_stat["pep_num"] = 0
            trans_detail_list[trans_id]= trans_detail
            trans_stat_list.append(trans_stat)
        gene_detail_list[gene_id] = gene_detail
        gene_stat_list.append(gene_stat)

    # gene_detail_df = pd.DataFrame(gene_detail_list,columns=["gene_id","trans_id","cds_id","pep_id"])
    gene_stat_df = pd.DataFrame(gene_stat_list,columns=["gene_id","trans_num","cds_num","pep_num"])
    with open(os.path.join(out,"gene_detail"),"w") as f:
        pickle.dump(gene_detail_list, f)
    # trans_detail_df = pd.DataFrame(trans_detail_list,columns=["trans_id","cds_id","pep_id"])
    trans_stat_df = pd.DataFrame(trans_stat_list,columns=["trans_id","cds_num","pep_num"])
    # gene_detail_df.to_csv(os.path.join(out,"gene_detail"),sep="\t",index=False)
    gene_stat_df.to_csv(os.path.join(out, "gene_stat"),sep="\t",index=False)
    with open(os.path.join(out,"trans_detail"),"w") as f:
        pickle.dump(trans_detail_list, f)
    # trans_detail_df.to_csv(os.path.join(out, "trans_detail"),sep="\t",index=False)
    trans_stat_df.to_csv(os.path.join(out, "trans_stat"),sep="\t",index=False)

    #合并各序列字典
    with open(os.path.join(out, "gene_seq"),"w") as f:
        pickle.dump(gene_seq, f)
    with open(os.path.join(out, "txpt_seq"),"w") as f:
        pickle.dump(txpt_seq, f)
    with open(os.path.join(out, "cds_seq"), "w") as f:
        cds_dict = dict(known_cds_dict.items() + new_cds_dict.items())
        pickle.dump(cds_dict, f)
    with open(os.path.join(out, "pep_seq"), "w") as f:
        pep_dict = dict(known_pep_dict.items() + new_pep_dict.items())
        pickle.dump(pep_dict, f)



@pkgsfuncdeco
def fasta_to_dict(fasta):
    def get_id(record_id):
        if '|' in record_id:
            seq_id = record_id.split('|')[0]
        elif '(+)' in record_id:
            seq_id = record_id.split('(+)')[0]
        elif '(-)' in record_id:
            seq_id = record_id.split('(-)')[0]
        elif '(.)' in record_id:
            seq_id = record_id.split('(.)')[0]
        else:
            seq_id = record_id.strip()
        return seq_id

    return {get_id(record.id): str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}


def getg2t(t2g):
    g2t = defaultdict(list)
    with open(t2g, "r") as t2g:
        g2t = defaultdict(list)
        for line in t2g.readlines():
                line = line.strip().split()
                trans = line[0]
                gene = line[1]
                g2t[gene].append(trans)
    return g2t


@pkgsfuncdeco
def bed_to_loc(bed):
    info = dict()
    for line in open(bed):
        if not line.strip():
            continue
        items = line.strip('\n').split('\t')
        if info.get(items[3]) is None:
            info[items[3]] = {'chr': items[0], 'strand': items[5], 'start': int(items[1]) + 1, 'end': int(items[2])}
        else:
            raise Exception('find {} in multi-line at {}'.format(items[3], bed))
    return info


@pkgsfuncdeco
def parse_biomart(biomart_path, biomart_type):
    if biomart_type == 'type1':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 2
        chromosome_ind = 8
        gene_type_ind = 16
        desc_ind = 7
        strand_ind = 11
        start_ind = 9
        end_ind = 10
        pep_id_ind = 6
    elif biomart_type == 'type2':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 2
        chromosome_ind = 6
        gene_type_ind = 14
        desc_ind = 5
        strand_ind = 9
        start_ind = 7
        end_ind = 8
        pep_id_ind = 4
    elif biomart_type == 'type3':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 0
        chromosome_ind = 4
        gene_type_ind = 12
        desc_ind = 3
        strand_ind = 7
        start_ind = 5
        end_ind = 6
        pep_id_ind = 2
    else:
        raise Exception('biomart type should be among type1, type2 and type3')
    biomart_info = dict()
    pep2transcript = dict()
    with open(biomart_path) as f:
        for line in f:
            if not line.strip():
                continue
            line = line.replace('\t\t', '\t-\t')
            items = line.strip('\n').split('\t')
            gene_id = items[gene_id_ind]
            trans_id = items[trans_id_ind]
            gene_name = items[gene_name_ind]
            if biomart_type == 'type3':
                gene_name = '-'
            chromosome = items[chromosome_ind]
            gene_type = items[gene_type_ind]
            desc = items[desc_ind]
            strand_tmp = items[strand_ind]
            if strand_tmp == '1':
                strand = '+'
            elif strand_tmp == '-1':
                strand = '-'
            elif strand_tmp == '0':
                strand = '.'
            else:
                strand = strand_tmp
            start = items[start_ind]
            end = items[end_ind]
            pep_id = items[pep_id_ind]
            biomart_info.setdefault(gene_id, defaultdict(list))
            biomart_info[gene_id]['trans_id'].append(trans_id)
            biomart_info[gene_id]['gene_name'].append(gene_name)
            biomart_info[gene_id]['chromosome'].append(chromosome)
            biomart_info[gene_id]['gene_type'].append(gene_type)
            biomart_info[gene_id]['description'].append(desc)
            biomart_info[gene_id]['pep_id'].append(pep_id)
            biomart_info[gene_id]['strand'].append(strand)
            biomart_info[gene_id]['start'].append(start)
            biomart_info[gene_id]['end'].append(end)
            pep2transcript[pep_id] = trans_id
    if not biomart_info:
        raise Exception('biomart information is None')
    print 'information of {} genes was parsed from {}'.format(len(biomart_info), biomart_path)
    return biomart_info, pep2transcript


@pkgsfuncdeco
def get_cds_seq(cds_path):
    cds_dict = dict()
    trans2cds = defaultdict(list)
    cds_pattern_match = re.compile(r'>([^\s]+)').match
    with open(cds_path, 'r') as f:
        j = 0
        trans_id, cds_id, cds_sequence = '', '', ''
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    trans2cds[trans_id].append(cds_id)
                    trans2cds[cds_id].append(cds_id)
                    seq_len = len(cds_sequence)
                    cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                    cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                cds_id = cds_pattern_match(line).group(1)
                if '.' in cds_id:
                    trans_id = cds_id[:cds_id.rfind('.')]
                else:
                    trans_id = cds_id
                cds_sequence = ''
            else:
                cds_sequence += line.strip()
        else:
            trans2cds[trans_id].append(cds_id)
            trans2cds[cds_id].append(cds_id)
            seq_len = len(cds_sequence)
            cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
            cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
    if not cds_dict:
        print 'CDS information is None'
    print 'information of {} cds was parsed from {}'.format(len(cds_dict), cds_path)
    return cds_dict,trans2cds


@pkgsfuncdeco
def get_pep_seq(pep_path, p2t):
    pep_dict = dict()
    trans2pep = defaultdict(list)
    pep_dict_pd = dict()
    j, trans_id, trans_id_else, pep_sequence, pep_id = 0, '', '', '', ''
    pep_pattern = re.compile(r'>([^\s]+)')
    trans_pattern = re.compile(r'transcript:([^\s]+)')
    with open(pep_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    seq_len = len(pep_sequence)
                    if trans_id:
                        trans2pep[trans_id].append(pep_id)
                        trans2pep[trans_id_else].append(pep_id)
                        pep_dict_pd[pep_id] =  dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                        pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                        pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
                pep_id = pep_pattern.match(line).group(1)
                try:
                    trans_id = trans_pattern.search(line).group(1)
                except Exception:
                    if pep_id not in p2t:
                        if '.' in pep_id:
                            pep_id_else = pep_id[:pep_id.rfind('.')]
                            if pep_id_else not in p2t:
                                print 'transcript id -> protein {} failed in biomart'.format(pep_id_else)
                                trans_id = None
                                continue
                            else:
                                trans_id = p2t[pep_id_else]
                        else:
                            print 'transcript id -> protein {} failed in biomart'.format(pep_id)
                            trans_id = None
                            continue
                    else:
                        trans_id = p2t[pep_id]
                if '.' in trans_id:
                    trans_id_else = trans_id[:trans_id.rfind('.')]
                else:
                    trans_id_else = trans_id
                pep_sequence = ''
            else:
                pep_sequence += line.strip()
        else:
            seq_len = len(pep_sequence)
            pep_dict_pd[pep_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
            trans2pep[trans_id].append(pep_id)
            trans2pep[trans_id_else].append(pep_id)
            pep_dict[trans_id] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
            pep_dict[trans_id_else] = dict(name=pep_id, sequence=pep_sequence, sequence_length=seq_len)
    if not pep_dict:
        print 'PEP information is None'
    print 'information of {} pep was parsed from {}'.format(len(pep_dict), pep_path)
    return pep_dict,trans2pep,pep_dict_pd


@pkgsfuncdeco
def get_predict_cds_pep(cds_file, pep_file):
    def parse_seq_file(seq_file):
        with open(seq_file) as f:
            j = 0
            seq_tag = tuple()
            sequence = ""
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        yield seq_tag, sequence
                    seq_tag = line.strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            else:
                yield seq_tag, sequence

    pep_parser = parse_seq_file(pep_file)
    pep_dict = dict()
    for pep_desc, pep_seq in pep_parser:
        pep_dict[pep_desc] = pep_seq
    t2cds_pep = dict()
    new_t2cds_dict = defaultdict(list)
    new_t2pep_dict = defaultdict(list)
    new_cds_dict = defaultdict(str)
    new_pep_dict = defaultdict(str)
    match = re.compile(r'>(.*?)::(.*?)::.*type:(.*?)\s+len:\d+.*:(\d+-\d+\(.\)).*').match
    cds_parser = parse_seq_file(cds_file)
    for cds_desc, cds_seq in cds_parser:
        g_id, t_id, _type, cds_pos = match(cds_desc).groups()
        pep_seq = pep_dict[cds_desc]
        seq_id = re.match(r'>(.*?)\s.*', cds_desc).groups()[0]
        t2cds_pep.setdefault(t_id, list())
        t2cds_pep[t_id].append((seq_id, _type, cds_seq, pep_seq))  # pep_id = cds_id
        new_t2cds_dict[t_id].append(seq_id)
        new_t2pep_dict[t_id].append(seq_id)
        new_cds_dict[seq_id] = cds_seq
        new_pep_dict[seq_id] = pep_seq
    return t2cds_pep,new_t2cds_dict,new_t2pep_dict,new_cds_dict,new_pep_dict


@pkgsfuncdeco
def build_database(db_path, gene_seq, txpt_seq, t2g_dict, known_cds, known_pep, new_cds_pep):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    for table_name, seq_dict in dict(gene_seq=gene_seq, trans_seq=txpt_seq).items():
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} (seq_id text, sequence text)'.format(table_name))
        for seq_id, seq in seq_dict.items():
            cursor.execute("INSERT INTO {} VALUES ('{}', '{}')".format(table_name, seq_id, seq))
    table_name = "trans_annot"
    table_columns = 'transcript_id text, cds_id text, pep_id text, cds_seq text, pep_seq text, orf_type text'
    cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
    cursor.execute('CREATE TABLE {} ({})'.format(table_name, table_columns))
    for transcript in t2g_dict:
        if transcript in known_cds:
            cds_id = known_cds[transcript]['name']
            cds_seq = known_cds[transcript]['sequence']
            orf_type = 'complete'
            if transcript in known_pep:
                pep_id = known_pep[transcript]['name']
                pep_seq = known_pep[transcript]['sequence']
                orf_type = 'complete'
            else:
                pep_id = 'None'
                pep_seq = 'None'
            cursor.execute("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}', '{}', '{}')".format(
                table_name, transcript, cds_id, pep_id, cds_seq, pep_seq, orf_type))
        elif transcript in new_cds_pep:
            for seq_id, _type, cds_seq, pep_seq in new_cds_pep[transcript]:
                cursor.execute("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}', '{}', '{}')".format(
                    table_name, transcript, seq_id, seq_id, cds_seq, pep_seq, _type))
        else:
            print '{} has no CDS and PEP annotation'.format(transcript)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script for making sqlite3 database of sequences')
    parser.add_argument('--gf', action='store', required=True, dest='gf',
                        help='gene FASTA file')
    parser.add_argument('--tf', action='store', required=True, dest='tf',
                        help='transcript FASTA file')
    parser.add_argument('--bf', action='store', required=True, dest='bf',
                        help='biomart file')
    parser.add_argument('--bt', action='store', required=True, dest='bt',
                        help='biomart type')
    parser.add_argument('--t2g', action='store', required=True, dest='t2g',
                        help='relation map between transcript and gene file')
    parser.add_argument('--kcf', action='store', required=True, dest='kcf',
                        help='known cds FASTA file')
    parser.add_argument('--kpf', action='store', required=True, dest='kpf',
                        help='known pep FASTA file')
    parser.add_argument('--ncf', action='store', dest='ncf',
                        help='novel cds FASTA file')
    parser.add_argument('--npf', action='store', dest='npf',
                        help='novel pep FASTA file')
    parser.add_argument('--output', action='store', dest='db',
                        help='output sqlite3 database file')
    parser.add_argument('--output_detail', action='store', dest='detail',
                        help='output detail database file')

    args = parser.parse_args()

    main(args)
