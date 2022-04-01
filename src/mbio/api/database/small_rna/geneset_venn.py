# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from api_base import ApiBase
from bson.objectid import ObjectId

class GenesetVenn(ApiBase):
    def __init__(self, bind_object):
        super(GenesetVenn, self).__init__(bind_object)

    def insert_key_value_pair(self, main_id, key='workflow', value='used'):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        conn = self.db['sg_geneset_venn']
        conn.update({'_id': main_id}, {'$set': {key: value}})
