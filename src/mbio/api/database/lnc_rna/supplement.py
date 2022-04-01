# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import unittest

class Supplement(ApiBase):
    '''
    last_modify: 2019.05.27
    '''
    def __init__(self, bind_object):
        super(Supplement, self).__init__(bind_object)

    @report_check
    def perfect_annotation_query_detail(self, query_id, known_tsv, novel_tsv, mrna_list, t2g_file, entrez):
        query_id = ObjectId(query_id)
        kdf = pd.read_table(known_tsv, na_filter='')
        kdf = kdf.rename({'lncrna_id': 'transcript_id', 'gene_description': 'description', 'lncrna_length': 'length'}, axis=1)
        kdf['seq_type'] = 'ref'
        ndf = pd.read_table(novel_tsv, na_filter='')
        ndf = ndf.rename({'gene_description': 'description', 'lncrna_length': 'length'}, axis=1)
        ndf['seq_type'] = 'new'
        columns = ['transcript_id', 'gene_id', 'gene_name', 'description', 'seq_type', 'length']
        df = pd.concat([kdf.reindex(columns, axis=1), ndf.reindex(columns, axis=1)], ignore_index=True)
        max_length_lncrnas = set()
        g2mlt_dict = dict()
        for i in df.index:
            if df.loc[i, 'gene_id'] in g2mlt_dict:
                if df.loc[i, 'length'] > g2mlt_dict[df.loc[i, 'gene_id']][1]:
                    g2mlt_dict.update({df.loc[i, 'gene_id']: (df.loc[i, 'transcript_id'], int(df.loc[i, 'length']))})
            else:
                g2mlt_dict.update({df.loc[i, 'gene_id']: (df.loc[i, 'transcript_id'], int(df.loc[i, 'length']))})
        mrnas = {line.strip() for line in open(mrna_list)}
        t2g_dict = dict([line.strip().split('\t') for line in open(t2g_file)])
        genes = {t2g_dict[mrna] for mrna in mrnas if mrna in t2g_dict}
        for i in df.index:
            if df.loc[i, 'gene_id'] in genes:
                df.loc[i, 'is_gene'] = False
            else:
                if df.loc[i, 'transcript_id'] == g2mlt_dict[df.loc[i, 'gene_id']][0]:
                    df.loc[i, 'is_gene'] = True
                else:
                    df.loc[i, 'is_gene'] = False
        df['query_id'] = query_id
        df['rna_type'] = 'lncRNA'
        df = df.merge(pd.read_table(entrez, usecols=[0, 2], names=['gene_id', 'enterz'], na_filter=''), on='gene_id', how='left')
        df = df.drop_duplicates()
        columns.extend(['is_gene', 'query_id', 'cog', 'cog_description', 'length', 'enterz', 'rna_type', 'ko_id', 'ko_name', 'pathways', 'pfam', 'go', 'nr', 'swissprot'])
        df = df.reindex(columns, axis=1).fillna('')
        self.create_db_table('sg_annotation_query_detail', df.to_dict('r'))
        self.update_db_record('sg_annotation_query', query_id, has_lnc=True)

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'supplement_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'lnc_rna.lnc_rna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = LncRnaTestApiWorkflow(wheet)
        wf.sheet.id = 'lnc_rna'
        wf.sheet.project_sn = 'lnc_rna'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('lnc_rna.supplement')
        wf.test_api.perfect_annotation_query_detail(
            query_id='5caaad7d17b2bf2f3dfe91d3',
            known_tsv='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_identify/output/known_lncrna_detail.tsv',
            novel_tsv='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_identify/output/novel_lncrna_detail.tsv',
            mrna_list='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_identify/output/mrna.list',
            t2g_file='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_identify/output/t2g.pairs',
            entrez='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_identify/output/Homo_sapiens.GRCh38.biomart_enterz.txt',
        )

if __name__ == '__main__':
    unittest.main()