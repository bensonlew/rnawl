# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import sys

from biocluster.config import Config
from dateutil.parser import parse

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
base = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish')

collections = {'sg_genome_db': db['sg_genome_db'],
               'sg_task': db['sg_task'],
               'sg_species_information': db['sg_species_information'],
               'sg_species_information_detail': db['sg_species_information_detail']}


def fix_sg_species_information(file_path, species_id, task_id):
    insert_data_list = list()
    with open(file_path) as fr:
        fr.readline()
        for line in fr:
            tmp = line.strip().split('\t')
            chr = tmp[0]
            size = tmp[1]
            gc = tmp[2]
            gene = tmp[3]
            protein_coding = tmp[4]
            other_rna = tmp[5]
            pseudo = tmp[6]
            insert_data = {
                'species_id': species_id,
                'ref_seq': chr,
                'length': size,
                'qc_percent': gc,
                'gene': gene,
                'proteincoding': protein_coding,
                'other_rna': other_rna,
                'pseudogene': pseudo
            }
            insert_data_list.append(insert_data)
    collection = db['sg_species_information_detail']
    collection.delete_many({'species_id': species_id})
    collection.insert_many(insert_data_list)
    logging.info('succeed in inserting records to {} by task_id ({})'.format(collection, task_id))

def single_operation(document):
    task_id = document['task_id']
    document = collections['sg_task'].find_one({'task_id': task_id})
    logging.info('start processing {} ...'.format(task_id))
    gtf_suffix = document['ref_gtf'].split('Genome_DB_finish/')[1]
    gene_stat = collections['sg_genome_db'].find_one({'gtf': gtf_suffix})['gene_stat']
    file_path = os.path.join(base, gene_stat)
    species_id = collections['sg_species_information'].find_one({'task_id': task_id})['main_id']
    fix_sg_species_information(file_path, species_id, task_id)

if len(sys.argv) == 2:
    task_id = sys.argv[1]
    document = collections['sg_task'].find_one({'task_id': task_id})
    single_operation(document)
elif len(sys.argv) == 3:
    task_id = sys.argv[1]
    genome_id = sys.argv[2]
    document = collections['sg_task'].find_one({'task_id': task_id})
    logging.info('start processing {} ...'.format(task_id))
    gene_stat = collections['sg_genome_db'].find_one({'genome_id': genome_id})['gene_stat']
    file_path = os.path.join(base, gene_stat)
    species_id = collections['sg_species_information'].find_one({'task_id': task_id})['main_id']
    fix_sg_species_information(file_path, species_id, task_id)
else:
    cutoff_ts = parse('2019-05-13 00:00:00')
    for document in collections['sg_task'].find({}):
        if parse(document['created_ts']) > cutoff_ts and 'ref_gtf' in document:
            task_id = document['task_id']
            logging.info('start processing {} ...'.format(task_id))
            gtf_suffix = document['ref_gtf'].split('Genome_DB_finish/')[1]
            gene_stat = collections['sg_genome_db'].find_one({'gtf': gtf_suffix})['gene_stat']
            file_path = os.path.join(base, gene_stat)
            species_id = collections['sg_species_information'].find_one({'task_id': task_id})['main_id']
            fix_sg_species_information(file_path, species_id, task_id)
    else:
        logging.info('succeed in fix issues!')
