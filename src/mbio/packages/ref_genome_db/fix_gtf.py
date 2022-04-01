# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import shutil

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
base = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish')


collection = db['sg_genome_db']
for document in collection.find({}):
    if int(document['genome_id'][-3:]) >= 273:
        logging.info('precessing {} ...'.format(document['genome_id']))
        gtf = os.path.join(base, document['gtf'])
        bak = '{}.bak'.format(gtf)
        shutil.copy(bak, gtf)
        logging.info('succeed in copying {} to {}'.format(bak, gtf))
