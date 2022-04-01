# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import pickle

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId
import unittest
from mbio.api.database.medical_transcriptome.api_base import ApiBase


class Masigpro(ApiBase):
    '''
    last_modify: 2019.07.19
    '''

    def __init__(self, bind_object):
        super(Masigpro, self).__init__(bind_object)

    @report_check
    def add_masigpro(self, result_dir, main_id, s3_output):
        cmp_list=[]
        heatmaps={}
        main_id = ObjectId(main_id)
        cmp_list=os.listdir(result_dir)
        if len(cmp_list) == 1:
            s3_output=s3_output+"/"+cmp_list[0]
            multi = None
        else:
            multi="true"
        if multi:
            for cmp in cmp_list:
                idf = self.get_id_map_data_frame(main_id)
                df = pd.read_table(os.path.join(result_dir,cmp,"result.tsv"))
                df = df.rename(lambda x: x.replace('.', '_'), axis=1)
                df = pd.merge(df, idf)
                df = df.drop_duplicates()
                rdf = df.copy()
                df = df.drop('description', axis=1)
                df['masigpro_id'] = main_id
                df["cmp"]= cmp
                heatmaps[cmp] = list({os.path.basename(i)[:-4] for i in glob.glob(
                    os.path.join(os.path.join(result_dir,cmp), 'heatmap.*')
                )})
                for data_pkl in glob.glob(os.path.join(result_dir,cmp, '*.pkl')):
                    seq_data = pickle.load(open(data_pkl))
                    group= os.path.splitext(os.path.basename(data_pkl))[0]
                    self.create_db_table('sg_masigpro_seq', seq_data, {'masigpro_id': main_id,"cmp":cmp,"group":group})
                self.create_db_table('sg_masigpro_detail', df.to_dict('r'))
                self.modify_output(rdf, os.path.join(result_dir,cmp,"result.tsv"))
        else:
            idf = self.get_id_map_data_frame(main_id)
            df = pd.read_table(os.path.join(result_dir, cmp_list[0], "result.tsv"))
            df = df.rename(lambda x: x.replace('.', '_'), axis=1)
            df = pd.merge(df, idf)
            df = df.drop_duplicates()
            rdf = df.copy()
            df = df.drop('description', axis=1)
            df['masigpro_id'] = main_id
            heatmaps = list({os.path.basename(i)[:-4] for i in glob.glob(
                os.path.join(os.path.join(result_dir,cmp_list[0]), 'heatmap.*')
            )})
            for data_pkl in glob.glob(os.path.join(result_dir, cmp_list[0], '*.pkl')):
                seq_data = pickle.load(open(data_pkl))
                group = os.path.splitext(os.path.basename(data_pkl))[0]
                self.create_db_table('sg_masigpro_seq', seq_data, {'masigpro_id': main_id,  "group": group})
            self.create_db_table('sg_masigpro_detail', df.to_dict('r'))
            self.modify_output(rdf, os.path.join(result_dir, cmp_list[0], "result.tsv"))

        if multi:
            self.update_db_record('sg_masigpro', main_id, main_id=main_id, s3_output=s3_output,
                                  heatmaps=heatmaps,cmp_list=cmp_list, status='end',multi="true")
        else:
            self.update_db_record('sg_masigpro', main_id, main_id=main_id, s3_output=s3_output,
                                  heatmaps=heatmaps, cmp_list=cmp_list, status='end')


        # main_id = ObjectId(main_id)
        # idf = self.get_id_map_data_frame(main_id)
        # df = pd.read_table(result_table)
        # df = df.rename(lambda x: x.replace('.', '_'), axis=1)
        # df = pd.merge(df, idf)
        # df = df.drop_duplicates()
        # rdf = df.copy()
        # df = df.drop('description', axis=1)
        # df['masigpro_id'] = main_id
        # heatmaps = list({os.path.basename(i)[:-4] for i in glob.glob(
        #     os.path.join(self.bind_object.output_dir, 'heatmap.*')
        # )})
        # self.create_db_table('sg_masigpro_detail', df.to_dict('r'))
        # self.update_db_record('sg_masigpro', main_id, main_id=main_id, s3_output=s3_output,
        #                       heatmaps=heatmaps, status='end')
        # self.modify_output(rdf, result_table)

    def get_id_map_data_frame(self, main_id):
        level = self.db['sg_masigpro'].find_one({'main_id': main_id})['level']
        task_id = self.db['sg_masigpro'].find_one({'main_id': main_id})['task_id']
        exp_id = self.db['sg_exp'].find_one({'task_id': task_id, 'level': level, 'is_rmbe': False})['main_id']
        cursor = self.db['sg_exp_detail'].find({'exp_id': exp_id})

        # query_id = self.db['sg_annotation_query'].find_one({'task_id': task_id})['main_id']
        # cursor = self.db['sg_annotation_query_detail'].find({'query_id': query_id})
        data = list()
        for doc in cursor:
            if level == 'T':
                data.append({'seq_id': doc['transcript_id'], 'gene_id': doc['gene_id'], 'gene_name': doc['gene_name'],
                             'description': doc['description']})
            elif level == 'G':
                data.append(
                    {'seq_id': doc['gene_id'], 'gene_name': doc['gene_name'], 'description': doc['description']})
        else:
            df = pd.DataFrame(data).drop_duplicates()
            return df

    def modify_output(self, df, result_table):
        os.rename(result_table,os.path.join(os.path.dirname(result_table),"result.xls"))
        result_table_new = os.path.join(os.path.dirname(result_table),"result.xls")
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
        df.to_csv(result_table_new, sep='\t', index=False)

    def run1(self):
        result_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200518/Single_masigpro_7653_2484/Masigpro/output"
        main_id="5d4a4e6a17b2bf73c182822f"
        self.add_masigpro(result_dir,main_id,"s3://refrnav2/files/m_188/188_5d3e93b81d8b3/tsg_34958/interaction_results/Masigpro_20190807_120706")




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    # test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/api/database/' \
    #            'denovo_rna_v2/test_files'
    # toolbox = AllExp(None)
    # task_id = "tsg_36994"

    def test_diff(self):
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "tsg_34958",
            "project_sn": "tsg_34958",
            "type": "workflow",
            "name": "ref_rna_v2.refrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("ref_rna_v2.masigpro")
        wf.test_api.run1()

if __name__ == '__main__':
        unittest.main()