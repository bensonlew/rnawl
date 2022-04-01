# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import shutil
import subprocess

import pandas as pd
from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
base = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish')
seqkit = os.path.join(Config().SOFTWARE_DIR, 'bioinfo/seq/seqkit')


def line2items(line):
    items = line.strip().split('\t')
    return items[0], round(float(items[-2]) / 10 ** 6, 2), items[-1]


collection = db['sg_genome_db']
for document in collection.find({}):
    if int(document['genome_id'][-3:]) >= 273:
        logging.info('processing {} ...'.format(document['genome_id']))
        fasta = os.path.join(base, document['dna_fa'])
        gtf = os.path.join(base, document['gtf'])
        fx2tab = '{}.fx2tab'.format(gtf)
        stat = os.path.join(base, document['gene_stat'])
        cmd = '{} fx2tab -g -l -n -i {} -o {}'.format(seqkit, fasta, fx2tab)
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = proc.wait()
        outs, errs = proc.communicate()
        if retcode:
            msg = '\n'.join(('fail to excecute {}'.format(cmd), 'STDOUT: {}'.format(outs), 'STDERR: {}'.format(errs)))
            raise Exception(msg)
        else:
            logging.info('succeed in excecuting command ({})'.format(cmd))
            logging.debug(outs)
            logging.error(errs)
        dfl = pd.DataFrame([line2items(line) for line in open(fx2tab)], columns=['Chr', 'Size(Mb)', 'GC%'])
        dfr = pd.read_table(stat, usecols=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
        df = pd.merge(dfl, dfr)
        bak = '{}.bak'.format(stat)
        os.rename(stat, bak)
        logging.info('succeed in moving {} to {}'.format(stat, bak))
        df.to_csv(stat, sep='\t', index=False)
        logging.info('succeed in exporting {}'.format(stat))
