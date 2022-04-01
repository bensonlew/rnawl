# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class GeneFusionAnnot(ApiBase):
    def __init__(self, bind_object):
        super(GeneFusionAnnot, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_fusion_annot(self, main_id, result_file_path):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(result_file_path, header=0, sep='\t')
        df.rename(columns={'Source': 'source', 'Gene': 'gene', 'Genome': 'genome', "Disease": "disease",
                           "PUBMED_PMID": "pubmed_id"}, inplace=True)
        use_cols = ["source", "gene", "genome", "disease", "pubmed_id"]
        df = df[use_cols]
        columns_list = list()
        columns_list.append({'field': 'source', 'filter': False, 'sort': False, 'title': 'Source', 'type': 'string'})
        columns_list.append({'field': 'gene', 'filter': False, 'sort': False, 'title': 'Gene', 'type': 'string'})
        columns_list.append({'field': 'genome', 'filter': False, 'sort': False, 'title': 'Genome', 'type': 'string'})
        columns_list.append({'field': 'disease', 'filter': False, 'sort': False, 'title': 'Disease', 'type': 'string'})
        columns_list.append({'field': 'pubmed_id', 'filter': False, 'sort': False, 'title': 'PMID', 'type': 'string'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        df['gene_fusion_annot_id'] = main_id
        detail = df.to_dict('r')
        if len(detail) > 0:
            has_annot = ["yes"]
            self.create_db_table('sg_gene_fusion_annot_detail', detail)
        else:
            has_annot = ["no"]
        self.update_db_record('sg_gene_fusion_annot', main_id, main_id=main_id, column_data_detail=columns_data,
                              has_annot=has_annot)



