# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.medical_transcriptome.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import unittest

class Stem(ApiBase):
    '''
    last_modify: 2019.07.02
    '''
    def __init__(self, bind_object):
        super(Stem, self).__init__(bind_object)

    @report_check
    def add_stem(self, detail_table, cluster_table, main_id):
        main_id = ObjectId(main_id)
        idf = self.get_id_map_data_frame(main_id)
        ddf = pd.read_table(detail_table)
        time = list(ddf.columns)[2:]
        ddf['stem_id'] = main_id
        ddf = pd.merge(ddf, idf)
        for c in time:
            ddf[c] = ddf[c].apply(lambda x: str(x).replace(',', '')).astype(float)
        self.create_db_table('sg_stem_detail', ddf.to_dict('r'))
        cdf = pd.read_table(cluster_table)
        cluster = sorted(set(cdf['cluster']))
        cdf['stem_id'] = main_id
        self.create_db_table('sg_stem_cluster', cdf.to_dict('r'))
        self.update_db_record('sg_stem', main_id, main_id=main_id, time=time, cluster=cluster, status='end')

    def get_id_map_data_frame(self, main_id):
        exp_level = self.db['sg_stem'].find_one({'main_id': main_id})['level']
        task_id = self.db['sg_stem'].find_one({'main_id': main_id})['task_id']
        exp_id = self.db['sg_exp'].find_one({'task_id': task_id, 'level': exp_level, 'is_rmbe': False})['main_id']
        cursor = self.db['sg_exp_detail'].find({'exp_id': exp_id})
        # query_id = self.db['sg_annotation_query'].find_one({'task_id': task_id})['main_id']
        # cursor = self.db['sg_annotation_query_detail'].find({'query_id': query_id})
        data = list()
        for doc in cursor:
            if exp_level == 'T':
                data.append({'seq_id': doc['transcript_id'], 'gene_id': doc['gene_id'], 'gene_name': doc['gene_name']})
            elif exp_level == 'G':
                data.append({'seq_id': doc['gene_id'], 'gene_name': doc['gene_name']})
        else:
            df = pd.DataFrame(data).drop_duplicates()
            return df
