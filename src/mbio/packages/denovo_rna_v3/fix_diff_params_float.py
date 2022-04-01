# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main():
    database = Config().get_mongo_client(mtype='denovo_rna_v2')[Config().get_mongo_dbname('denovo_rna_v2')]
    collection = database['sg_diff']
    cursor = collection.find({'status': 'end'})
    count = 0
    for document in cursor:
        _id = document['_id']
        oldparams = document['params']
        try:
            dct = json.loads(oldparams)
            if 'fc' in dct:
                dct['fc'] = str(float(dct['fc']))
                if 'tpm_filter_threshold' in dct:
                    dct['tpm_filter_threshold'] = str(float(dct['tpm_filter_threshold']))
                newparams = json.dumps(dct, sort_keys=True, separators=(',', ':'))
                collection.update({'_id': _id}, {'$set': {'params': newparams}}, upsert=True)
                # logging.info('succeed in processing document by _id ({})'.format(_id))
                count += 1
        except:
            pass
        if not count % 10000:
            logging.info('succeed in fixing {} documents'.format(count))
    else:
        logging.info('succeed in fixing {} documents'.format(count))


if __name__ == '__main__':
    main()
