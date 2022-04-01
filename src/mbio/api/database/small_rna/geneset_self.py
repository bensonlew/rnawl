# -*- coding: utf-8 -*-
# __author__ = "fengyitong, qinjincheng"

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.small_rna.api_base import ApiBase
import os
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict
from biocluster.config import Config
from bson.objectid import ObjectId

class GenesetSelf(ApiBase):
    def __init__(self, bind_object):
        super(GenesetSelf, self).__init__(bind_object)
        self._project_type = 'small_rna'

    def add_geneset(self, geneset_table, name=None, gene_type="M", main_id=None):
        main_id = ObjectId(main_id)
        with open(geneset_table) as f:
            _ = f.readline()
            genes = [x.strip() for x in f]
        row_dict_list = [{'seq_list': genes}]
        tag_dict = dict(geneset_id=main_id)
        self.create_db_table('sg_geneset_detail', row_dict_list, tag_dict=tag_dict)
        gene_length=len(genes)
        if gene_length < 1:
            self.update_db_record('sg_geneset', record_id=main_id, status='failed', params=None)
            self.bind_object.set_error('预上传序列集中基因/转录本为空，不予上传，请核查ID书写是否规范或者从交互页面直接创建')
        else:
            self.update_db_record('sg_geneset', record_id=main_id, status='end', gene_length=gene_length)
            self.bind_object.logger.info('succeed in building geneset named {}'.format(
                self.db['sg_geneset'].find_one({'main_id': main_id})['name']
            ))
