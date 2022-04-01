# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Stem(ApiBase):
    '''
    last_modify: 2020.01.07
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
        self.create_db_table('stem_detail', ddf.to_dict('r'))
        cdf = pd.read_table(cluster_table)
        cluster = sorted(set(cdf['cluster']))
        cdf['stem_id'] = main_id
        self.create_db_table('stem_cluster', cdf.to_dict('r'))
        self.update_db_record('stem', main_id, main_id=main_id, time=time, cluster=cluster, status='end')

    def get_id_map_data_frame(self, main_id):
        exp_level = self.db['stem'].find_one({'main_id': main_id})['level']
        task_id = self.db['stem'].find_one({'main_id': main_id})['task_id']
        query_id = self.db['exp'].find_one({'task_id': task_id, 'level': exp_level})['main_id']
        cursor = self.db['exp_detail'].find({'exp_id': query_id})
        data = list()
        for doc in cursor:
            if exp_level == 'T':
                try:
                    gene_id = doc['gene_id']
                except:
                    gene_id = ""
                try:
                    gene_name = doc['gene_name']
                except:
                    gene_name = ""
                try:
                    description = doc['description']
                except:
                    description = ""
                data.append({'seq_id': doc['transcript_id'], 'gene_id': gene_id, 'gene_name': gene_name,
                             'description': description})
            elif exp_level == 'G':
                try:
                    gene_name = doc['gene_name']
                except:
                    gene_name = ""
                try:
                    description = doc['description']
                except:
                    description = ""
                data.append({'seq_id': doc['gene_id'], 'gene_name': gene_name, 'description': description})
        else:
            df = pd.DataFrame(data).drop_duplicates()
            return df
