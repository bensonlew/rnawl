# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import re

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

pg = re.compile(r'gene_id "(\S+)";')
pt = re.compile(r'transcript_id "(\S+)";')


def main(args):
    ref_gtf = args.ref
    merged_gtf = args.merge
    tmap = args.tmap
    output_dir = args.output

    ref_gene_id_set = set()
    ref_transcript_id_set = set()
    ref_record_list = list()
    new_transcript_id2class_code_dict = dict()
    new_record_list = list()
    new_gene_id_set = set()
    old_record_list = list()

    for line in open(ref_gtf):
        record = line.strip()
        if 'gene_id' in record:
            gene_id = get_gene_id(record)
            ref_gene_id_set.add(gene_id)
            if 'transcript_id' in record:
                transcript_id = get_transcript_id(record)
                ref_transcript_id_set.add(transcript_id)
                ref_record = add_class_code(record, '=')
            else:
                ref_record = record
            ref_record_list.append(ref_record)
    else:
        logging.info('succeed in finding {} ref genes in {}'.format(len(ref_gene_id_set), ref_gtf))
        logging.info('succeed in finding {} ref transcripts in {}'.format(len(ref_transcript_id_set), ref_gtf))

    tmap_df = pd.read_table(tmap)
    for index, row in tmap_df.iterrows():
        if row['qry_id'] not in ref_transcript_id_set:
            new_transcript_id2class_code_dict[row['qry_id']] = row['class_code']
    else:
        logging.info('succeed in finding {} new transcripts in {}'.format(len(new_transcript_id2class_code_dict), tmap))

    for line in open(merged_gtf):
        record = line.strip()
        if 'gene_id' in record and 'transcript_id' in record:
            transcript_id = get_transcript_id(record)
            if transcript_id in new_transcript_id2class_code_dict:
                class_code = new_transcript_id2class_code_dict[transcript_id]
                if class_code in list('ijoxu'):
                    new_record = add_class_code(record, class_code)
                    new_record_list.append(new_record)
                    gene_id = get_gene_id(new_record)
                    if gene_id not in ref_gene_id_set:
                        new_gene_id_set.add(gene_id)
                else:
                    if 'ref_gene_id' in record:
                        gene_id, ref_gene_id = re.findall("gene_id\s(.*?)\;", record)
                        record = record.replace(gene_id, ref_gene_id)
                        old_record = add_class_code(record, '=')
                        old_record_list.append(old_record)
            else:
                if 'ref_gene_id' in record:
                    gene_id, ref_gene_id = re.findall("gene_id\s(.*?)\;", record)
                    record = record.replace(gene_id, ref_gene_id)
                    old_record = add_class_code(record, '=')
                    old_record_list.append(old_record)
    else:
        logging.info(
            'succeed in finding {} records in {} by screening class code'.format(len(new_record_list), merged_gtf))

    all_gtf = os.path.join(output_dir, 'all.gtf')
    ref_gtf = os.path.join(output_dir, 'ref.gtf')
    new_gtf = os.path.join(output_dir, 'new.gtf')
    old_gtf = os.path.join(output_dir, 'old.gtf')
    with open(ref_gtf, 'w') as ref_gtf_handle, open(new_gtf, 'w') as new_gtf_handle, \
            open(all_gtf, 'w') as all_gtf_handle, open(old_gtf, 'w') as old_gtf_handle:
        for ref_record in ref_record_list:
            all_gtf_handle.write('{}\n'.format(ref_record))
            ref_gtf_handle.write('{}\n'.format(ref_record))
        else:
            logging.info('succeed in writing {} records to {}'.format(len(ref_record_list), all_gtf))
            logging.info('succeed in writing {} records to {}'.format(len(ref_record_list), ref_gtf))
        for new_record in new_record_list:
            all_gtf_handle.write('{}\n'.format(new_record))
            new_gtf_handle.write('{}\n'.format(new_record))
        else:
            logging.info('succeed in writing {} records to {}'.format(len(new_record_list), all_gtf))
            logging.info('succeed in writing {} records to {}'.format(len(new_record_list), new_gtf))
        for old_record in old_record_list:
            old_gtf_handle.write('{}\n'.format(old_record))
        else:
            logging.info('succeed in writing {} records to {}'.format(len(new_record_list), old_gtf))

    open(os.path.join(output_dir, 'new.gene_id.list'), 'w').writelines(
        '{}\n'.format(gene_id) for gene_id in new_gene_id_set)
    logging.info('succeed in writing {} records to {}'.format(
        len(new_gene_id_set), os.path.join(output_dir, 'new.gene_id.list')))

    transcript_id2gene_id_dict = get_transcript_id2gene_id_dict(all_gtf)
    open(os.path.join(output_dir, 't2g.txt'), 'w').writelines(
        '{}\t{}\n'.format(transcript_id, gene_id) for transcript_id, gene_id in transcript_id2gene_id_dict.items())
    logging.info('succeed in exporting {}'.format(os.path.join(output_dir, 't2g.txt')))
    open(os.path.join(output_dir, 'trans_type.xls'), 'w').writelines(
        '{}\t{}\tmRNA\t{}\n'.format(
            transcript_id, gene_id, 'known' if transcript_id in ref_transcript_id_set else 'novel')
        for transcript_id, gene_id in transcript_id2gene_id_dict.items())
    logging.info('succeed in exporting {}'.format(os.path.join(output_dir, 'trans_type.xls')))

    gene_id2transcript_id_set = get_gene_id2transcript_id_set_dict(transcript_id2gene_id_dict)
    open(os.path.join(output_dir, 'gene_type.xls'), 'w').writelines(
        '{}\t{}\tmRNA\t{}\n'.format(
            gene_id, ';'.join(transcript_id_set), 'known' if gene_id in ref_gene_id_set else 'novel')
        for gene_id, transcript_id_set in gene_id2transcript_id_set.items())
    logging.info('succeed in exporting {}'.format(os.path.join(output_dir, 'gene_type.xls')))


add_class_code = lambda record, class_code: '{} class_code "{}";'.format(record, class_code)


def get_gene_id(record):
    mg = re.search(pg, record)
    if mg:
        return mg.group(1)
    else:
        raise Exception('can not find gene_id in record -> ({})'.format(record))


def get_transcript_id(record):
    mt = re.search(pt, record)
    if mt:
        return mt.group(1)
    else:
        raise Exception('can not find transcript_id in record -> ({})'.format(record))


def get_transcript_id2gene_id_dict(all_gtf):
    transcript_id2gene_id_dict = dict()
    for line in open(all_gtf):
        if not line.strip():
            continue
        if line[0] == '#':
            continue
        eles = line.strip().split('\t')
        if len(eles) != 9:
            continue
        mg = re.search(pg, eles[8])
        mt = re.search(pt, eles[8])
        if mg and mt:
            gene_id = mg.group(1)
            transcript_id = mt.group(1)
            transcript_id2gene_id_dict[transcript_id] = gene_id
    return transcript_id2gene_id_dict


def get_gene_id2transcript_id_set_dict(transcript_id2gene_id_dict):
    gene_id2transcript_id_set_dict = dict()
    for transcript_id, gene_id in transcript_id2gene_id_dict.items():
        if gene_id in gene_id2transcript_id_set_dict:
            gene_id2transcript_id_set_dict[gene_id].add(transcript_id)
        else:
            gene_id2transcript_id_set_dict[gene_id] = {transcript_id}
    return gene_id2transcript_id_set_dict


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate requried GTF by majorbio rules')
    parser.add_argument('-r', action='store', required=True,
                        help='reference GTF file', metavar='<FILE>', dest='ref')
    parser.add_argument('-m', action='store', required=True,
                        help='merged GTF file', metavar='<FILE>', dest='merge')
    parser.add_argument('-t', action='store', required=True,
                        help='compare TMAP file', metavar='<FILE>', dest='tmap')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
