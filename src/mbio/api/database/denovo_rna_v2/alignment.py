# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.denovo_rna_v2.api_base import ApiBase

TITLE_DICT = {
    'genome': {
        'db_title': ['query_id', 'query_start', 'query_end', 'genomic_location', 'overlapping_gene_id',
                     'overlapping_gene_name', 'overlapping_gene_description', 'orientation', 'length', 'score',
                     'evalue', 'identity'],
        'pg_title': ['Query_id', 'Query_start', 'Query_end', 'Genomic_location', 'Overlapping_gene_id',
                     'Overlapping_gene_name', 'Overlapping_gene_description', 'Orientation', 'Length', 'Score',
                     'Evalue', 'Identity(%)']
    },
    'database': {
        'db_title': ['query_id', 'query_start', 'query_end', 'subject_id', 'gene_id', 'gene_name', 'gene_description',
                     'subject_start', 'subject_end', 'orientation', 'length', 'score', 'evalue', 'identity'],
        'pg_title': ['Query_id', 'Query_start', 'Query_end', 'Subject_id', 'Subject Gene_id', 'Subject Gene_name',
                     'Subject Gene_description', 'Subject_start', 'Subject_end', 'Orientation', 'Length', 'Score',
                     'Evalue', 'Identity(%)']
    },
    'upload': {
        'db_title': ['query_id', 'query_start', 'query_end', 'subject_id', 'subject_start', 'subject_end', 'length',
                     'score', 'evalue', 'identity'],
        'pg_title': ['Query_id', 'Query_start', 'Query_end', 'Subject_id', 'Subject_start', 'Subject_end', 'Length',
                     'Score', 'Evalue', 'Identity(%)']
    }
}


class Alignment(ApiBase):
    '''
    last_modify: 2019.08.12
    '''

    def __init__(self, bind_object):
        super(Alignment, self).__init__(bind_object)

    @report_check
    def add_blast_detail(self, result, about, main_id):
        main_id = ObjectId(main_id)
        df = pd.read_table(result)
        df['blast_id'] = main_id
        df = df.fillna('')
        self.create_db_table('sg_blast_detail', df.to_dict('r'))
        self.update_db_record('sg_blast', main_id, about=about, db_title=TITLE_DICT[about]['db_title'],
                              pg_title=TITLE_DICT[about]['pg_title'], main_id=main_id, status='end')
