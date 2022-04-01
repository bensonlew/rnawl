# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import os

class GoEnrich(ApiBase):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, bind_object):
        super(GoEnrich, self).__init__(bind_object)

    @report_check
    def add_go_enrich(self, result, main_id, s3_output):
        df = pd.read_table(result)
        categories = list(set(df['go_type']))
        result_dir = os.path.join(s3_output, os.path.basename(result))
        main_id = ObjectId(main_id)
        df['seq_str'] = df['seq_list'].apply(lambda x: x.split(';'))
        df['go_enrich_id'] = main_id
        self.create_db_table('sg_geneset_go_enrich_detail', df.to_dict('r'))
        self.update_db_record('sg_geneset_go_enrich', main_id, categories=categories,
                              result_dir=result_dir, main_id=main_id, status='end')

