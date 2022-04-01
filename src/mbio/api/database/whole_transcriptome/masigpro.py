# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import pickle

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Masigpro(ApiBase):
    '''
    last_modify: 2019.07.19
    '''

    def __init__(self, bind_object):
        super(Masigpro, self).__init__(bind_object)

    @report_check
    def add_masigpro(self, result_table, main_id, s3_output):
        main_id = ObjectId(main_id)
        idf = self.get_id_map_data_frame(main_id)
        df = pd.read_table(result_table)
        df = df.rename(lambda x: x.replace('.', '_'), axis=1)
        left_df = df.set_index('seq_id')
        right_df = idf.set_index('seq_id')
        join_df = left_df.join(right_df)
        df = join_df.reset_index()
        # df = pd.merge(df, idf)
        rdf = df.copy()
        df = df.drop('description', axis=1)
        df['masigpro_id'] = main_id
        heatmaps = list({os.path.basename(i)[:-4] for i in glob.glob(
            os.path.join(self.bind_object.output_dir, 'heatmap.*')
        )})
        for data_pkl in glob.glob(os.path.join(self.bind_object.output_dir, '*.pkl')):
            seq_data = pickle.load(open(data_pkl))
            self.create_db_table('masigpro_seq', seq_data, {'masigpro_id': main_id})
        self.create_db_table('masigpro_detail', df.to_dict('r'))
        self.update_db_record('masigpro', main_id, main_id=main_id, s3_output=s3_output,
                              heatmaps=heatmaps, status='end')
        self.modify_output(rdf, result_table)

    def get_id_map_data_frame(self, main_id):
        exp_level = self.db['masigpro'].find_one({'main_id': main_id})['level']
        task_id = self.db['masigpro'].find_one({'main_id': main_id})['task_id']
        exp_id = self.db['exp'].find_one({'task_id': task_id, 'level': exp_level})['main_id']
        cursor = self.db['exp_detail'].find({'exp_id': exp_id})
        data = list()
        for doc in cursor:
            if exp_level == 'T':
                try:
                    gene_id = doc['gene_id']
                except:
                    gene_id = str()
                try:
                    gene_name = doc['gene_name']
                except:
                    gene_name = str()
                try:
                    description = doc['description']
                except:
                    description = str()
                data.append({'seq_id': doc['transcript_id'], 'gene_id': gene_id, 'gene_name': gene_name,
                             'description': description})
            elif exp_level == 'G':
                try:
                    gene_name = doc['gene_name']
                except:
                    gene_name = str()
                try:
                    description = doc['description']
                except:
                    description = str()
                data.append({'seq_id': doc['gene_id'], 'gene_name': gene_name, 'description': description})
        else:
            df = pd.DataFrame(data).drop_duplicates()
            return df

    def modify_output(self, df, result_table):
        dct = {'gene_name': 'Gene Name', 'description': 'gene description', 'p_value': 'p-value',
               'p_adjust': 'p-adjust', 'r_squared': 'r.squared', 'cluster': 'cluster'}
        if 'gene_id' in df.columns:
            dct.update({'seq_id': 'Transcript ID', 'gene_id': 'Gene ID'})
            index = ['Transcript ID']
        else:
            dct.update({'seq_id': 'Gene ID'})
            index = list()
        index.extend(['Gene ID', 'Gene Name', 'gene description', 'p-value', 'p-adjust', 'r.squared', 'cluster'])
        df = df.rename(dct, axis=1)
        df = df.reindex(index, axis=1)
        df.to_csv(result_table, sep='\t', index=False)
