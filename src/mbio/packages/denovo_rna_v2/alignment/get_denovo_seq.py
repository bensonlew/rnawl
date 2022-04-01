# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

import pandas as pd
from Bio import SeqIO
from biocluster.config import Config
from biocluster.file import download
from bson.objectid import ObjectId

from mbio.packages.denovo_rna_v2.functions import watcher


@watcher
def main(args):
    {'query': get_query, 'subject': get_subject}[args.subparsers](args)


@watcher
def get_query(args):
    database = Config().get_mongo_client(mtype='denovo_rna_v2',db_version =args.version)[Config().get_mongo_dbname('denovo_rna_v2',db_version =args.version)]
    task_document = database['sg_task'].find_one({'task_id': args.task})
    if args.geneset != 'All':
        geneset_detail_document = database['sg_geneset_detail'].find_one({'geneset_id': ObjectId(args.geneset)})
        seq_list = geneset_detail_document['seq_list']
    output_dir = os.path.dirname(args.output)
    if args.unit == 'nucl':
        from_file = os.path.join(task_document['assembly_dir'],
                                 task_document['assembly_object'][{'T': 'filter', 'G': 'unigene'}[args.level]])
        to_file = os.path.join(output_dir, os.path.basename(from_file))
        logging.debug('start downing {} to {}'.format(from_file, to_file))
        download(from_file, to_file)
        if args.geneset == 'All':
            os.rename(to_file, args.output)
        else:
            logging.debug('start selecting required sequences from specified geneset')
            SeqIO.write([s for s in SeqIO.parse(to_file, 'fasta') if s.id in seq_list], args.output, 'fasta')
    elif args.unit == 'prot':
        s3_pep_file = os.path.join(task_document['assembly_dir'], task_document['assembly_object']['peptide'])
        my_pep_file = os.path.join(output_dir, os.path.basename(s3_pep_file))
        logging.debug('start downing {} to {}'.format(s3_pep_file, my_pep_file))
        download(s3_pep_file, my_pep_file)
        if args.geneset == 'All':
            os.rename(my_pep_file, args.output)
        else:
            s3_map_file = os.path.join(task_document['assembly_dir'], task_document['assembly_object']['map'])
            my_map_file = os.path.join(output_dir, os.path.basename(s3_map_file))
            logging.debug('start downing {} to {}'.format(s3_map_file, my_map_file))
            download(s3_map_file, my_map_file)
            seq_set = collect_required_peptide(my_map_file, args.level, seq_list)
            logging.debug('start selecting required sequences from specified geneset')
            SeqIO.write([s for s in SeqIO.parse(my_pep_file, 'fasta') if s.id in seq_set], args.output, 'fasta')


# @watcher
def collect_required_peptide(map_file, exp_level, seq_list):
    data = list()
    for line in open(map_file):
        items = line.strip().split('\t')
        if len(items) == 4:
            data.append(items + [str()])
        elif len(items) == 5:
            data.append(items)
    else:
        df = pd.DataFrame(data, columns=['transcript_id', 'gene_id', 'is_gene', 'length', 'peptide_ids'])
    expr = '{} in @seq_list'.format({'T': 'transcript_id', 'G': 'gene_id'}[exp_level])
    df = df.query(expr)
    peptide_id_set = set()
    peptide_id_set.update(*df['peptide_ids'].apply(lambda x: x.split(';')))
    peptide_id_set.remove(str())
    return peptide_id_set


@watcher
def get_subject(args):
    database = Config().get_mongo_client(mtype='denovo_rna_v2',db_version =args.version)[Config().get_mongo_dbname('denovo_rna_v2',db_version =args.version)]
    task_document = database['sg_task'].find_one({'task_id': args.task})
    output_dir = os.path.dirname(args.output)
    from_file = os.path.join(task_document['assembly_dir'], task_document['assembly_object'][args.source])
    to_file = os.path.join(output_dir, os.path.basename(from_file))
    logging.debug('start downing {} to {}'.format(from_file, to_file))
    download(from_file, to_file)
    os.rename(to_file, args.output)


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Obtain assembly fasta file')
    subparsers = parser.add_subparsers(dest='subparsers', help='sub-command help')

    parser_q = subparsers.add_parser('query', help='use file as query')
    parser_q.add_argument('--task', action='store', required=True,
                          help='task id', metavar='<STR>', dest='task')
    parser_q.add_argument('--unit', action='store', choices=['nucl', 'prot'], required=True,
                          help='sequence class', dest='unit')
    parser_q.add_argument('--level', action='store', choices=['T', 'G'], required=True,
                          help='separator type', dest='level')
    parser_q.add_argument('--geneset', action='store', required=True,
                          help='geneset id', metavar='<STR>', dest='geneset')
    parser_q.add_argument('--output', action='store', required=True,
                          help='output file path', metavar='<FILE>', dest='output')
    parser_q.add_argument('-v', '--version', dest='version', required=True, type=int, help='db_version')

    parser_s = subparsers.add_parser('subject', help='use file as subject')
    parser_s.add_argument('--task', action='store', required=True,
                          help='task id', metavar='<STR>', dest='task')
    parser_s.add_argument('--source', action='store', choices=['filter', 'raw'], required=True,
                          help='fasta file source', dest='source')
    parser_s.add_argument('--output', action='store', required=True,
                          help='output file path', metavar='<FILE>', dest='output')
    parser_s.add_argument('-v', '--version', dest='version', required=True, type=int, help='db_version')

    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    if len(sys.argv) > 1 and sys.argv[1] not in ['-h', '--help']:
        if sys.argv[1] in ['query', 'subject']:
            args = parser.parse_args()
        else:
            logging.error('unrecognized command ({})'.format(sys.argv[1]))
            sys.exit(-2)
    else:
        parser.print_help()
        sys.exit(-1)

    main(args)
