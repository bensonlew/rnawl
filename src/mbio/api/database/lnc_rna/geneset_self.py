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
        self._project_type = 'lnc_rna'

    def add_geneset(self, geneset_table, task_id, name=None, gene_type="M", main_id=None):
        main_id = ObjectId(main_id)
        genes = self.db['sg_genes'].find_one({'task_id': task_id})
        if genes and "main_id" in genes:
            genes_main_id = genes['main_id']
        else:
            self.bind_object.set_error('can not find task id {}'.format(task_id))

        with open(geneset_table) as f:
            genes = [x.strip() for x in f]
        genes_choose = []
        genes_delete = []
        for gene in genes:
            if gene_type == "G":
                gene_find = self.db['sg_genes_detail'].find_one({'genes_id': genes_main_id, 'gene_id': gene, "level": "G", "gene_type": "mRNA"})
            if gene_type == "LG":
                gene_find = self.db['sg_genes_detail'].find_one({'genes_id': genes_main_id, 'gene_id': gene, "level": "G", "gene_type": "lncRNA"})
            if gene_type == "T":
                gene_find = self.db['sg_genes_detail'].find_one({'genes_id': genes_main_id, 'trans_id': gene, "level": "T", "gene_type": "mRNA"})
            if gene_type == "LT":
                gene_find = self.db['sg_genes_detail'].find_one({'genes_id': genes_main_id, 'trans_id': gene, "level": "T", "gene_type": "lncRNA"})
            if gene_find:
                genes_choose.append(gene)
            else:
                genes_delete.append(gene)

        row_dict_list = [{'seq_list': genes_choose}]
        tag_dict = dict(geneset_id=main_id)
        self.create_db_table('sg_geneset_detail', row_dict_list, tag_dict=tag_dict)
        gene_length=len(genes_choose)
        if gene_length < 1:
            self.update_db_record('sg_geneset', record_id=main_id, status='failed', gene_length=gene_length, params=None)
            # self.bind_object.set_error('no id in this task of uploaded seq list ')
            self.bind_object.set_error('预上传序列集中基因/转录本为空，不予上传，请核查ID书写是否规范或者从交互页面直接创建')
        elif len(genes_delete) != 0:
            self.update_db_record('sg_geneset', record_id=main_id, status='end', gene_length=gene_length, desc="列表中部分基因不满足条件，自动剔除")
            self.bind_object.logger.info('列表中部分基因不满足条件，自动剔除{}'.format(genes_delete))
        else:
            self.update_db_record('sg_geneset', record_id=main_id, status='end', gene_length=gene_length)
            self.bind_object.logger.info('succeed in building geneset named {}'.format(
                self.db['sg_geneset'].find_one({'main_id': main_id})['name']
            ))
