# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import logging
import os
import pickle

from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

mapped1_fq_fp = '/mnt/ilustre/users/isanger/test/mapped_data/O2_zt2_1.R1.fq'
raw_data_dir = '/mnt/ilustre/users/isanger/test/raw_data'
output_dir = '/mnt/ilustre/users/isanger/test/mapped_raw_data'

if __name__ == '__main__':
    if os.path.isfile('tmp.set.pkl'):
        mapped_id_set = pickle.load(open('tmp.set.pkl'))
    else:
        logger.debug('Parse {}'.format(mapped1_fq_fp))
        mapped_id_set = {r.id for r in SeqIO.parse(mapped1_fq_fp, 'fastq')}
        pickle.dump(mapped_id_set, open('tmp.set.pkl', 'w'))
    logger.debug('Mapped reads number: {}'.format(len(mapped_id_set)))
    records_dict = {'_R1': list(), '_R2': list()}
    for fp in glob.glob(os.path.join(raw_data_dir, '*.fastq')):
        logger.debug('Parse {}'.format(fp))
        for record in SeqIO.parse(fp, 'fastq'):
            if record.id in mapped_id_set:
                for k in records_dict:
                    if k in os.path.basename(fp):
                        records_dict[k].append(record)
    pickle.dump(records_dict, open('tmp.dict.pkl', 'w'))
    logger.debug('Sort records and export ...')
    for k, records in records_dict.items():
        records.sort(key=lambda r: r.id)
        SeqIO.write(records, os.path.join(output_dir, '1-zt_CAGATC_L003{}.fastq'.format(k)), 'fastq')
