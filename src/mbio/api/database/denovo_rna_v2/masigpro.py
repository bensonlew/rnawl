# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.denovo_rna_v2.api_base import ApiBase


class Masigpro(ApiBase):
    '''
    last_modify: 2019.08.13
    '''

    def __init__(self, bind_object):
        super(Masigpro, self).__init__(bind_object)

    @report_check
    def add_masigpro(self, result_table, main_id, s3_output):
        main_id = ObjectId(main_id)
        idf = self.get_id_map_data_frame(main_id)
        idf = idf.set_index('seq_id')
        df = pd.read_table(result_table)
        df = df.rename(lambda x: x.replace('.', '_'), axis=1)
        df = df.set_index('seq_id')
        df = df.join(idf)
        df = df.reset_index()
        rdf = df.copy()
        df = df.drop('description', axis=1)
        df['masigpro_id'] = main_id
        heatmaps = list({os.path.basename(i)[:-4] for i in glob.glob(
            os.path.join(self.bind_object.output_dir, 'heatmap.*')
        )})
        self.create_db_table('sg_masigpro_detail', df.to_dict('r'))
        self.update_db_record('sg_masigpro', main_id, main_id=main_id, s3_output=s3_output,
                              heatmaps=heatmaps, status='end')
        self.modify_output(rdf, result_table)

    def get_id_map_data_frame(self, main_id):
        exp_level = self.db['sg_masigpro'].find_one({'main_id': main_id})['exp_level']
        task_id = self.db['sg_masigpro'].find_one({'main_id': main_id})['task_id']
        nr_id = self.db['sg_annotation_nr'].find_one({'task_id': task_id})['main_id']
        data = list()
        if exp_level == 'T':
            cursor = self.db['sg_annotation_nr_detail'].find({'nr_id': nr_id})
            for doc in cursor:
                data.append(
                    {'seq_id': doc['transcript_id'], 'gene_id': doc['gene_id'], 'description': doc['description']})
        else:
            cursor = self.db['sg_annotation_nr_detail'].find({'nr_id': nr_id, 'is_gene': True})
            for doc in cursor:
                data.append({'seq_id': doc['gene_id'], 'description': doc['description']})
        df = pd.DataFrame(data).drop_duplicates()
        return df

        # cursor = self.db['sg_annotation_nr_detail'].find({'nr_id': nr_id, 'is_gene': True})
        # data = list()
        # for doc in cursor:
        #     if exp_level == 'T':
        #         data.append(
        #             {'seq_id': doc['transcript_id'], 'gene_id': doc['gene_id'], 'description': doc['description']})
        #     elif exp_level == 'G':
        #         data.append({'seq_id': doc['gene_id'], 'description': doc['description']})
        # else:
        #     df = pd.DataFrame(data).drop_duplicates()
        #     return df

    def modify_output(self, df, result_table):
        dct = {'description': 'NR description', 'p_value': 'p-value',
               'p_adjust': 'p-adjust', 'r_squared': 'r.squared', 'cluster': 'cluster'}
        if 'gene_id' in df.columns:
            dct.update({'seq_id': 'Transcript ID', 'gene_id': 'Gene ID'})
            index = ['Transcript ID']
        else:
            dct.update({'seq_id': 'Gene ID'})
            index = list()
        index.extend(['Gene ID', 'gene description', 'p-value', 'p-adjust', 'r.squared', 'cluster'])
        df = df.rename(dct, axis=1)
        df = df.reindex(index, axis=1)
        df.to_csv(result_table, sep='\t', index=False)
