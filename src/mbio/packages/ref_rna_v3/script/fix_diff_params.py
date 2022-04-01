# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main():
    database = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
    collection = database['sg_diff']
    cursor = collection.find({'version': 'v3'})
    count = 0
    for document in cursor:
        _id = document['_id']
        oldparams = document['params']
        try:
            dct = json.loads(oldparams)
            if 'correct_method' not in dct:
                dct.update({'correct_method': 'BH'})
                newparams = json.dumps(dct, sort_keys=True, separators=(',', ':'))
                collection.update({'_id': _id}, {'$set': {'params': newparams}}, upsert=True)
                logging.info('succeed in processing document by _id ({})'.format(_id))
                count += 1
        except:
            pass
    else:
        logging.info('succeed in fixing {} documents'.format(count))


if __name__ == '__main__':
    main()
